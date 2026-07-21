# Copied from Results_ALL-FEM-main/reference solutions/fluid/11/fluid11.py for dataset Problem_27.
# Original benchmark: Fluid Mechanics Problem 11 (Medium).

# file: powerlaw_channel_steady.py
from __future__ import print_function
from dolfin import *
import numpy as np
import csv
import json

# -----------------------
# Geometry and mesh
# -----------------------
L, H = 2.0, 0.20
nx, ny = 240, 24
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), nx, ny, "crossed")  # triangles are robust

# -----------------------
# Parameters (SI)
# -----------------------
rho  = 1.0          # kg/m^3
mu0  = 0.5          # Pa*s^n (consistency index)
nexp = 0.5          # flow behaviour index (shear-thinning)
Ubar = 1.0          # mean inlet speed (m/s)

# Regularisation to avoid singular viscosity at zero shear
eps_reg = Constant(1.0e-6)

# -----------------------
# Boundary markers
# -----------------------
left_id, right_id, bot_id, top_id = 1, 2, 3, 4
facets = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)

tol = 1e-12
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0, tol)
class Right(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], L, tol)
class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.0, tol)
class Top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], H, tol)

Left().mark(facets, left_id)
Right().mark(facets, right_id)
Bottom().mark(facets, bot_id)
Top().mark(facets, top_id)

ds_ = Measure("ds", domain=mesh, subdomain_data=facets)

# -----------------------
# Function spaces (P2-P1)
# -----------------------
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)
mixed_el = MixedElement([V.ufl_element(), Q.ufl_element()])
W = FunctionSpace(mesh, mixed_el)

# Scalar CG1 for projections (viscosity, speed, etc.)
S = FunctionSpace(mesh, "CG", 1)

# Unknowns and tests
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

# Picard iterate fields
u_k = Function(V)                  # advecting velocity
mu_k = Function(S)                 # frozen effective viscosity

# -----------------------
# Inlet profile (parabolic with mean Ubar)
# u_x(y) = 6 Ubar y (H - y) / H^2, u_y = 0
# -----------------------
inlet_expr = Expression(
    ("6.0*Ubar*x[1]*(H - x[1])/(H*H)", "0.0"),
    degree=2, Ubar=Ubar, H=H
)

# -----------------------
# Boundary conditions
# -----------------------
bc_inlet_u = DirichletBC(W.sub(0), inlet_expr, facets, left_id)
bc_walls_u = DirichletBC(W.sub(0), Constant((0.0, 0.0)), facets, bot_id)
bc_wallt_u = DirichletBC(W.sub(0), Constant((0.0, 0.0)), facets, top_id)

# Pressure gauge at a single point on outlet to fix nullspace
class GaugePoint(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], L, 1e-10) and near(x[1], H/2.0, 1e-3)
gauge = GaugePoint()
bc_p0 = DirichletBC(W.sub(1), Constant(0.0), gauge, method="pointwise")

bcs = [bc_inlet_u, bc_walls_u, bc_wallt_u, bc_p0]

# -----------------------
# Helper functions
# -----------------------
def D(u):
    return sym(grad(u))

def shear_rate(u):
    # |D| = sqrt(2 D:D)
    return sqrt(2.0*inner(D(u), D(u)) + eps_reg**2)

def mu_eff_of(u):
    # mu_eff = mu0 * |D|^(n-1)
    return mu0*pow(shear_rate(u), nexp - 1.0)

# -----------------------
# Picard linearised variational problem
# Convection uses u_k; viscosity uses mu_k
# Outlet traction-free is natural in this weak form
# -----------------------
a = (
    rho*inner(grad(u)*u_k, v)*dx
  + 2.0*inner(mu_k*D(u), D(v))*dx
  - p*div(v)*dx
  + q*div(u)*dx
)
Lform = Constant(0.0)*q*dx  # no body force here

w = Function(W)

# Initial guess: Stokes with constant viscosity mu0 and zero advection
u_k.assign(Constant((0.0, 0.0)))
mu_k.assign(Constant(mu0))

# -----------------------
# Picard iteration
# -----------------------
max_it = 40
rtol = 1e-8
atol = 1e-10
conv = False

for it in range(1, max_it+1):
    solve(a == Lform, w, bcs, solver_parameters={
        "linear_solver": "mumps"
    })
    u_sol, p_sol = w.split(deepcopy=True)

    # Update viscosity from new velocity
    mu_new = project(mu_eff_of(u_sol), S)

    # Convergence check on velocity
    if u_k.vector().norm("l2") > 0:
        rel = (u_sol.vector() - u_k.vector()).norm("l2") / u_k.vector().norm("l2")
    else:
        rel = (u_sol.vector()).norm("l2")

    # Assign updates
    u_k.assign(u_sol)
    mu_k.assign(mu_new)

    print("Picard iter %d: rel_u = %.3e" % (it, rel))

    if rel < rtol or (u_sol.vector().norm("l2") < atol):
        conv = True
        break

if not conv:
    print("Warning: Picard did not reach rtol=%.1e in %d iterations." % (rtol, max_it))

# Final fields
u_sol, p_sol = u_k, p_sol
mu_eff = mu_k

# -----------------------
# Outputs
# -----------------------

# 1) Save to XDMF
xdmf = XDMFFile(mesh.mpi_comm(), "q12_solution.xdmf")
xdmf.parameters["flush_output"] = True
xdmf.parameters["functions_share_mesh"] = True
u_sol.rename("u", "velocity")
p_sol.rename("p", "pressure")
mu_eff.rename("mu_eff", "effective_viscosity")
xdmf.write(u_sol, 0.0)
xdmf.write(p_sol, 0.0)
xdmf.write(mu_eff, 0.0)
xdmf.close()
print("Saved (u, p, mu_eff) to q12_solution.xdmf")

# 2) Speed colour map |u|
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation

    V1 = FunctionSpace(mesh, "CG", 1)
    speed = project(sqrt(dot(u_sol, u_sol)), V1)
    coords = mesh.coordinates()
    cells = mesh.cells()

    if cells.shape[1] == 3:
        triang = Triangulation(coords[:,0], coords[:,1], cells)
        z = speed.compute_vertex_values(mesh)
        plt.figure(figsize=(10, 10*H/L))
        tpc = plt.tripcolor(triang, z, shading="gouraud")
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlim(0.0, L); plt.ylim(0.0, H)
        plt.xlabel("x (m)"); plt.ylabel("y (m)")
        cbar = plt.colorbar(tpc); cbar.set_label("|u| (m/s)")
        plt.title("Speed magnitude |u| (power-law n=%.2f)" % nexp)
    else:
        pplt = plot(speed); plt.colorbar(pplt)

    plt.tight_layout()
    plt.savefig("q12_speed.png", dpi=220)
    plt.close()
    print("Saved q12_speed.png")
except Exception as e:
    print("Plotting failed:", e)

# 3) Mid-channel profile u_x(y) at x = L/2
x_mid = 1.0
Ny = 201
ys = np.linspace(0.0, H, Ny)
ux_vals = []
# Slight epsilon shift inside the domain to avoid boundary lookup issues
xq = min(max(x_mid, 0.0 + 1e-12), L - 1e-12)

for y in ys:
    yq = min(max(float(y), 0.0 + 1e-12), H - 1e-12)
    ux_vals.append(u_sol(Point(xq, yq))[0])

# Save CSV
with open("q12_profile.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["y (m)", "u_x (m/s) @ x=%.3f" % x_mid])
    for y, ux in zip(ys, ux_vals):
        writer.writerow([y, ux])

# Report maximum on that line
ux_max = np.max(ux_vals)
y_at_max = ys[int(np.argmax(ux_vals))]
print("Max u_x on x = L/2: %.6f m/s at y = %.6f m" % (ux_max, y_at_max))
with open("max_velocity_profile.json", "w", encoding="utf-8") as f:
    json.dump({
        "line_x_m": float(x_mid),
        "max_u_x_m_per_s": float(ux_max),
        "y_at_max_m": float(y_at_max),
        "profile_csv": "q12_profile.csv",
    }, f, indent=2)
    f.write("\n")
# [[END_OF_CODE]]

# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
import json

mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

u_sol.rename("u", "velocity field")
p_sol.rename("p", "pressure field")
mu_eff.rename("mu_eff", "effective viscosity field")

solution_file = XDMFFile(mesh.mpi_comm(), "solution.xdmf")
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u_sol, 0.0)
solution_file.write(p_sol, 0.0)
solution_file.write(mu_eff, 0.0)
solution_file.close()

checkpoint_file = XDMFFile(mesh.mpi_comm(), "read_checkpoint.xdmf")
checkpoint_file.write_checkpoint(u_sol, "u", 0.0, XDMFFile.Encoding.HDF5, False)
checkpoint_file.write_checkpoint(p_sol, "p", 0.0, XDMFFile.Encoding.HDF5, True)
checkpoint_file.write_checkpoint(mu_eff, "mu_eff", 0.0, XDMFFile.Encoding.HDF5, True)
checkpoint_file.close()

read_checkpoint_metadata = {
    "format": "FEniCS XDMFFile.write_checkpoint",
    "time": 0.0,
    "fields": [
        {
            "symbol": "u_sol",
            "checkpoint_name": "u",
            "meaning": "velocity field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'u', -1)",
        },
        {
            "symbol": "p_sol",
            "checkpoint_name": "p",
            "meaning": "pressure field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'p', -1)",
        },
        {
            "symbol": "mu_eff",
            "checkpoint_name": "mu_eff",
            "meaning": "effective viscosity field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'mu_eff', -1)",
        },
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")
