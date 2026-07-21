# Copied from Results_ALL-FEM-main/reference solutions/fluid/13/fluid13.py for dataset Problem_29.
# Original benchmark: Fluid Mechanics Problem 13 (Medium).

# filename: q_advection_diffusion_transport.py
# legacy FEniCS (dolfin/mshr-style), steady advection–diffusion in a 2-D channel
from dolfin import *

# ----------------------------
# Parameters
# ----------------------------
L  = 1.0           # m
H  = 0.10          # m
Ux_max = 0.75      # m/s
D  = 1.0e-5        # m^2/s

nx, ny = 100, 10   # mesh cells (uniform 100 × 10)

# ----------------------------
# Mesh and function space
# ----------------------------
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), nx, ny)
V = FunctionSpace(mesh, "CG", 1)

# ----------------------------
# Velocity field: u = (u_x(y), 0),  u_x(y) = Umax * 4*y*(H - y)/H^2
# ----------------------------
u_expr = Expression(("4.0*Umax*x[1]*(H - x[1])/(H*H)", "0.0"), Umax=Ux_max, H=H, degree=2)
u = as_vector((u_expr[0], u_expr[1]))

# ----------------------------
# Boundary conditions
#   Γ_in  (x=0): c = 0
#   Γ_out (x=L): c = 1
#   Γ_w   (y=0 or y=H): zero diffusive normal flux (natural in weak form)
# ----------------------------
tol = 1e-12

class Inlet(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0, tol)

class Outlet(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], L, tol)

bc_in  = DirichletBC(V, Constant(0.0), Inlet())
bc_out = DirichletBC(V, Constant(1.0), Outlet())
bcs = [bc_in, bc_out]

# ----------------------------
# Peclet number with mesh size as length scale
#   Use h = min(L/nx, H/ny), U_char = Ux_max
#   Pe = U_char * h / (2*D)
# ----------------------------
h_ref = min(L/float(nx), H/float(ny))
Pe_ref = Ux_max * h_ref / (2.0*D)
print("Reference cell size h = %.6e m" % h_ref)
print("Reference Peclet number Pe = %.6e" % Pe_ref)

# ----------------------------
# Variational problem
#   Find c \in V:
#     u·grad(c) - D*Δc = 0 in Ω
#     c = 0 on Γ_in,  c = 1 on Γ_out,  (D ∇c·n = 0) on Γ_w
# ----------------------------
c  = TrialFunction(V)
v  = TestFunction(V)

a = ( dot(u, grad(c))*v + D*dot(grad(c), grad(v)) )*dx
Lform = Constant(0.0)*v*dx

# ----------------------------
# SUPG stabilisation (streamline-upwind) if Pe_ref > 1
#   Add: τ * (u·∇v) * [u·∇c - D Δc] dx
#   τ ≈ h_K / (2 |u|)   (guarded for |u| ~ 0)
# ----------------------------
use_supg = Pe_ref > 1.0
if use_supg:
    hK = CellDiameter(mesh)
    u_norm = sqrt(dot(u, u)) + DOLFIN_EPS
    tau = hK / (2.0*u_norm)
    residual_c = dot(u, grad(c)) - D*div(grad(c))
    a += tau * dot(u, grad(v)) * residual_c * dx
    # RHS remains zero since the PDE has zero source

print("SUPG enabled?" , use_supg)

# ----------------------------
# Solve
# ----------------------------
c_h = Function(V)
solve(a == Lform, c_h, bcs, solver_parameters={"linear_solver": "lu"})

# ----------------------------
# Output (XDMF)
# ----------------------------
cfile = XDMFFile(mesh.mpi_comm(), "q_transport_solution.xdmf")
cfile.parameters["flush_output"] = True
cfile.parameters["functions_share_mesh"] = True
cfile.write(c_h, 0.0)

# Also write mesh once (optional but handy for some viewers)
mfile = XDMFFile(mesh.mpi_comm(), "q_transport_mesh.xdmf")
mfile.write(mesh)
print("Saved: q_transport_solution.xdmf (and mesh).")

# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
import json

mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

c_h.rename("c", "concentration field")

solution_file = XDMFFile(mesh.mpi_comm(), "solution.xdmf")
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(c_h, 0.0)
solution_file.close()

checkpoint_file = XDMFFile(mesh.mpi_comm(), "read_checkpoint.xdmf")
checkpoint_file.write_checkpoint(c_h, "c", 0.0, XDMFFile.Encoding.HDF5, False)
checkpoint_file.close()

read_checkpoint_metadata = {
    "format": "FEniCS XDMFFile.write_checkpoint",
    "time": 0.0,
    "fields": [
        {
            "symbol": "c_h",
            "checkpoint_name": "c",
            "meaning": "concentration field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'c', -1)",
        },
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")
