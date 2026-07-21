# Copied from Results_ALL-FEM-main/reference solutions/solid/10/solid10.py for dataset Problem_10.
# Original benchmark: Solid Mechanics Problem 10 (Medium).

from dolfin import *
import matplotlib
matplotlib.use("Agg")  # headless-safe
import matplotlib.pyplot as plt

# ----------------------------
# Geometry and mesh (100 x 20)
# ----------------------------
Lx, Ly = 1.0, 0.20
nx, ny = 100, 20
mesh = RectangleMesh(Point(0.0, 0.0), Point(Lx, Ly), nx, ny, "right/left")

# ----------------------------
# Material (plane strain)
# ----------------------------
E  = 5.0e6            # Pa
nu = 0.49
mu     = E/(2.0*(1.0 + nu))
lmbda  = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))  # plane strain λ

# ----------------------------
# Function spaces (Taylor–Hood)
#   u in [P2]^2, p in P1
# ----------------------------
Ve = VectorElement("Lagrange", mesh.ufl_cell(), 2)
Qe = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
W  = FunctionSpace(mesh, MixedElement([Ve, Qe]))

# ----------------------------
# Boundary conditions
# ----------------------------
tol = 1e-14
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0, tol)
class Right(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], Lx, tol)

left  = Left()
right = Right()

ux_prescribed = 0.03  # metres
zero_vec = Constant((0.0, 0.0))
right_disp = Constant((ux_prescribed, 0.0))

bcs = [
    DirichletBC(W.sub(0), zero_vec, left),
    DirichletBC(W.sub(0), right_disp, right),
]
# Top/bottom are traction-free => natural BCs; nothing to do.

# ----------------------------
# Mixed formulation (Herrmann)
#   Find (u,p) s.t.
#       2μ ε(u):ε(v) - p div(v) - q div(u) - (1/λ) p q = (f, v)
#   with f = 0 here.
# ----------------------------
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

def eps(w):
    return sym(grad(w))

a = (2.0*mu*inner(eps(u), eps(v)) - p*div(v) - q*div(u) - (1.0/lmbda)*p*q) * dx
L = Constant(0.0) * q * dx  # no body force

# ----------------------------
# Solve
# ----------------------------
w = Function(W)
# Prefer a robust direct solver for the indefinite mixed system
solver_params = {"linear_solver": "mumps"}  # falls back if unavailable
try:
    solve(a == L, w, bcs, solver_parameters=solver_params)
except:
    solve(a == L, w, bcs)  # default fallback

u_sol, p_sol = w.split(deepcopy=True)

# ----------------------------
# Outputs
# ----------------------------
# 1) Displacement magnitude |u| as a colour map -> q10_disp.png
Vsig = FunctionSpace(mesh, "Lagrange", 1)
u_mag = project(sqrt(dot(u_sol, u_sol)), Vsig)

plt.figure()
c = plot(u_mag, title="Displacement magnitude |u| (m)")
plt.colorbar(c)
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.tight_layout()
plt.savefig("q10_disp.png", dpi=300)

# 2) Displacement field in XDMF
with XDMFFile(mesh.mpi_comm(), "q10_displacement.xdmf") as xdmf:
    xdmf.parameters["flush_output"] = True
    xdmf.parameters["functions_share_mesh"] = True
    xdmf.write(u_sol)

# (Optional) also write pressure if you wish to inspect incompressibility:
# with XDMFFile(mesh.mpi_comm(), "q10_pressure.xdmf") as xdmf_p:
#     xdmf_p.parameters["flush_output"] = True
#     xdmf_p.parameters["functions_share_mesh"] = True
#     xdmf_p.write(p_sol)

print("Done. Wrote: q10_disp.png and q10_displacement.xdmf")

# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
import json

mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

u_sol.rename("u", "displacement field")
p_sol.rename("p", "pressure-like Lagrange multiplier field")

solution_file = XDMFFile(mesh.mpi_comm(), "solution.xdmf")
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u_sol, 0.0)
solution_file.write(p_sol, 0.0)
solution_file.close()

checkpoint_file = XDMFFile(mesh.mpi_comm(), "read_checkpoint.xdmf")
checkpoint_file.write_checkpoint(u_sol, "u", 0.0, XDMFFile.Encoding.HDF5, False)
checkpoint_file.write_checkpoint(p_sol, "p", 0.0, XDMFFile.Encoding.HDF5, True)
checkpoint_file.close()

read_checkpoint_metadata = {
    "format": "FEniCS XDMFFile.write_checkpoint",
    "time": 0.0,
    "fields": [
        {
            "symbol": "u_sol",
            "checkpoint_name": "u",
            "meaning": "displacement field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'u', -1)",
        },
        {
            "symbol": "p_sol",
            "checkpoint_name": "p",
            "meaning": "pressure-like Lagrange multiplier field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'p', -1)",
        },
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")
