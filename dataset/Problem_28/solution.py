# Copied from Results_ALL-FEM-main/reference solutions/fluid/12/fluid12.py for dataset Problem_28.
# Original benchmark: Fluid Mechanics Problem 12 (Hard).

# file: ns_tdep_viscosity_channel.py
from __future__ import print_function
from dolfin import *
import numpy as np
import csv
import json

# -----------------------
# Geometry and mesh
# -----------------------
L, H = 2.0, 0.20
nx, ny = 240, 24  # structured triangulation
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), nx, ny, "crossed")

# -----------------------
# Parameters (SI)
# -----------------------
rho    = 1.0          # kg m^-3
Ubar   = 1.0          # m s^-1
mu_ref = 0.02         # Pa s
beta   = 0.05         # K^-1
T_ref  = 300.0        # K
kappa  = 1.0e-3       # m^2 s^-1  (thermal diffusivity)

# Small regularisation to avoid zero-div etc. (not strictly required here)
eps_reg = Constant(1e-12)

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
# Mixed space: (u, p, T) = (P2^2, P1, P1)
# -----------------------
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)
S = FunctionSpace(mesh, "Lagrange", 1)
mixed_el = MixedElement([V.ufl_element(), Q.ufl_element(), S.ufl_element()])
W = FunctionSpace(mesh, mixed_el)

# Unknowns and tests
w = Function(W)
(u, p, T) = split(w)
(v, q, s) = TestFunctions(W)

# -----------------------
# Inlet profiles
# -----------------------
# u_x(y) = 6*Ubar*y*(H-y)/H^2, u_y = 0
uin = Expression(("6.0*Ubar*x[1]*(H - x[1])/(H*H)", "0.0"), degree=2, Ubar=Ubar, H=H)

# -----------------------
# Boundary conditions
# -----------------------
# Flow: inlet velocity; no-slip on y=0 and y=H; outlet traction-free (natural)
bc_u_in  = DirichletBC(W.sub(0), uin,     facets, left_id)
bc_u_bot = DirichletBC(W.sub(0), Constant((0.0, 0.0)), facets, bot_id)
bc_u_top = DirichletBC(W.sub(0), Constant((0.0, 0.0)), facets, top_id)

# Temperature: T = T_ref at inlet; T = T_ref+10 on bottom; top/outlet are natural (zero normal flux)
bc_T_in   = DirichletBC(W.sub(2), Constant(T_ref), facets, left_id)
bc_T_bot  = DirichletBC(W.sub(2), Constant(T_ref + 10.0), facets, bot_id)

# Pressure gauge at a point on outlet
class GaugePoint(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], L, 1e-10) and near(x[1], H/2.0, 1e-3)
gauge = GaugePoint()
bc_p0 = DirichletBC(W.sub(1), Constant(0.0), gauge, method="pointwise")

bcs = [bc_u_in, bc_u_bot, bc_u_top, bc_T_in, bc_T_bot, bc_p0]

# -----------------------
# Model definitions
# -----------------------
def eps(u):
    return sym(grad(u))

# Temperature-dependent viscosity
muT = mu_ref*exp(-beta*(T - T_ref))

I = Identity(2)

# Weak forms:
# Momentum: rho*(grad(u)*u) : v + 2*mu(T)*eps(u):eps(v) - p*div(v) = 0
# Continuity: q*div(u) = 0
# Energy: (u·grad T)*s + kappa*grad(T)·grad(s) = 0
F_mom = rho*inner(grad(u)*u, v)*dx + 2.0*inner(muT*eps(u), eps(v))*dx - p*div(v)*dx
F_cont = q*div(u)*dx
F_temp = inner(dot(u, grad(T)), s)*dx + kappa*inner(grad(T), grad(s))*dx

F = F_mom + F_cont + F_temp

# Jacobian for Newton
J = derivative(F, w, TrialFunction(W))

# -----------------------
# Solve nonlinear system (Newton)
# -----------------------
problem = NonlinearVariationalProblem(F, w, bcs, J)
solver = NonlinearVariationalSolver(problem)
prm = solver.parameters
prm["nonlinear_solver"] = "newton"
# The former 1e-8 relative tolerance stopped after two Newton updates while
# the coupled temperature field still undershot the 300 K boundary minimum.
# A tighter tolerance converges the monolithic (u, p, T) system in four updates.
prm["newton_solver"]["relative_tolerance"] = 1e-12
prm["newton_solver"]["absolute_tolerance"] = 1e-10
prm["newton_solver"]["maximum_iterations"] = 50
prm["newton_solver"]["linear_solver"] = "mumps"
solver.solve()

u_sol, p_sol, T_sol = w.split(deepcopy=True)

# Also build mu(x,y) from final T
S1 = FunctionSpace(mesh, "CG", 1)
mu_field = project(mu_ref*exp(-beta*(T_sol - T_ref)), S1)
mu_field.rename("mu", "viscosity")

# -----------------------
# Outputs
# -----------------------

# 1) XDMF fields
xdmf = XDMFFile(mesh.mpi_comm(), "q13_solution.xdmf")
xdmf.parameters["flush_output"] = True
xdmf.parameters["functions_share_mesh"] = True
u_sol.rename("u", "velocity")
p_sol.rename("p", "pressure")
T_sol.rename("T", "temperature")
xdmf.write(u_sol, 0.0)
xdmf.write(p_sol, 0.0)
xdmf.write(T_sol, 0.0)
xdmf.write(mu_field, 0.0)
xdmf.close()
print("Saved fields to q13_solution.xdmf")

# 2) Save μ(x,y) colour map
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation

    coords = mesh.coordinates()
    cells = mesh.cells()

    mu_v = mu_field.compute_vertex_values(mesh)
    if cells.shape[1] == 3:
        triang = Triangulation(coords[:,0], coords[:,1], cells)
        plt.figure(figsize=(10, 10*H/L))
        tpc = plt.tripcolor(triang, mu_v, shading="gouraud")
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlim(0.0, L); plt.ylim(0.0, H)
        plt.xlabel("x (m)"); plt.ylabel("y (m)")
        cbar = plt.colorbar(tpc); cbar.set_label("μ (Pa·s)")
        plt.title("Temperature-dependent viscosity μ(x,y)")
    else:
        pplt = plot(mu_field); plt.colorbar(pplt)
    plt.tight_layout()
    plt.savefig("q13_mu.png", dpi=220)
    plt.close()
    print("Saved q13_mu.png")
except Exception as e:
    print("Plotting failed:", e)

# 3) Mid-length profile u_x(y) at x = 1.0 -> CSV
x_mid = 1.0
Ny = 201
ys = np.linspace(0.0, H, Ny)
ux_vals = []
xq = min(max(x_mid, 0.0 + 1e-12), L - 1e-12)
for y in ys:
    yq = min(max(float(y), 0.0 + 1e-12), H - 1e-12)
    ux_vals.append(u_sol(Point(xq, yq))[0])

with open("q13_profile.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["y", "ux"])
    for y, ux in zip(ys, ux_vals):
        writer.writerow([y, ux])

ux_max = float(np.max(ux_vals))
y_at_max = float(ys[int(np.argmax(ux_vals))])
print("Max u_x at x = 1.0: %.6f m/s at y = %.6f m" % (ux_max, y_at_max))
with open("max_velocity_profile.json", "w", encoding="utf-8") as f:
    json.dump({
        "line_x_m": float(x_mid),
        "max_u_x_m_per_s": float(ux_max),
        "y_at_max_m": float(y_at_max),
        "profile_csv": "q13_profile.csv",
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
T_sol.rename("T", "temperature field")
mu_field.rename("mu", "viscosity field")

solution_file = XDMFFile(mesh.mpi_comm(), "solution.xdmf")
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u_sol, 0.0)
solution_file.write(p_sol, 0.0)
solution_file.write(T_sol, 0.0)
solution_file.write(mu_field, 0.0)
solution_file.close()

checkpoint_file = XDMFFile(mesh.mpi_comm(), "read_checkpoint.xdmf")
checkpoint_file.write_checkpoint(u_sol, "u", 0.0, XDMFFile.Encoding.HDF5, False)
checkpoint_file.write_checkpoint(p_sol, "p", 0.0, XDMFFile.Encoding.HDF5, True)
checkpoint_file.write_checkpoint(T_sol, "T", 0.0, XDMFFile.Encoding.HDF5, True)
checkpoint_file.write_checkpoint(mu_field, "mu", 0.0, XDMFFile.Encoding.HDF5, True)
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
            "symbol": "T_sol",
            "checkpoint_name": "T",
            "meaning": "temperature field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'T', -1)",
        },
        {
            "symbol": "mu_field",
            "checkpoint_name": "mu",
            "meaning": "viscosity field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'mu', -1)",
        },
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")
