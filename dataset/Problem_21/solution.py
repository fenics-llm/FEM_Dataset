# Copied from Results_ALL-FEM-main/reference solutions/fluid/5/fluid5.py for dataset Problem_21.
# Original benchmark: Fluid Mechanics Problem 5 (Easy).

# q5_lid_driven_navier_stokes.py
# Legacy FEniCS (dolfin): steady incompressible Navier–Stokes (Taylor–Hood P2–P1)
from dolfin import *
import matplotlib.pyplot as plt

# --- Mesh ---
nx = ny = 128
mesh = UnitSquareMesh(nx, ny)

# --- Parameters ---
rho = Constant(1.0)     # kg m^-3
mu  = Constant(0.01)    # Pa·s

# --- FE spaces: Taylor–Hood ---
Ve = VectorElement("Lagrange", mesh.ufl_cell(), 2)
Pe = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
W  = FunctionSpace(mesh, MixedElement([Ve, Pe]))

# --- Boundary subsets ---
tol = DOLFIN_EPS

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 1.0, tol)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.0, tol)

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0, tol)

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 1.0, tol)

top    = Top()
bottom = Bottom()
left   = Left()
right  = Right()

# --- Unknowns/tests ---
w  = Function(W)
u, p = split(w)
v, q = TestFunctions(W)

def eps(u):
    return sym(grad(u))

# Nonlinear residual for steady Navier–Stokes
F = (rho*inner(grad(u)*u, v) + 2.0*mu*inner(eps(u), eps(v)) - p*div(v) + q*div(u)) * dx

# Jacobian
dw = TrialFunction(W)
J  = derivative(F, w, dw)

# --- Pressure pin for uniqueness: p(0,0) = 0 ---
class PPoint(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0, tol) and near(x[1], 0.0, tol)

pin_bc = DirichletBC(W.sub(1), Constant(0.0), PPoint(), method="pointwise")

# --- Continuation on lid speed to improve robustness at Re~100 ---
# We solve for increasing lid speeds: 0.2 -> 0.4 -> ... -> 1.0
lid_speeds = [0.2, 0.4, 0.6, 0.8, 1.0]

# Initial guess
assign(w.sub(0), project(Constant((0.0, 0.0)), W.sub(0).collapse()))
assign(w.sub(1), project(Constant(0.0),            W.sub(1).collapse()))

# Apply boundary conditions and solve for each lid speed
for Utop in lid_speeds:
    u_top = Constant((Utop, 0.0))
    bcs = [
        DirichletBC(W.sub(0), u_top,  top),
        DirichletBC(W.sub(0), Constant((0.0, 0.0)), bottom),
        DirichletBC(W.sub(0), Constant((0.0, 0.0)), left),
        DirichletBC(W.sub(0), Constant((0.0, 0.0)), right),
        pin_bc
    ]

    problem = NonlinearVariationalProblem(F, w, bcs, J)
    solver = NonlinearVariationalSolver(problem)
    prm = solver.parameters
    prm["nonlinear_solver"] = "newton"
    prm["newton_solver"]["maximum_iterations"] = 30
    prm["newton_solver"]["relative_tolerance"] = 1e-8
    prm["newton_solver"]["absolute_tolerance"] = 1e-10
    prm["newton_solver"]["report"] = True
    prm["newton_solver"]["linear_solver"] = "mumps"
    solver.solve()

# Split final solution
u_sol, p_sol = w.split(True)
u_sol.rename("u", "velocity")
p_sol.rename("p", "pressure")

# --- Post-processing: speed and colour map ---
V1 = FunctionSpace(mesh, "Lagrange", 1)
speed = project(sqrt(dot(u_sol, u_sol)), V1)
speed.rename("speed", "|u|")

plt.figure()
pm = plot(speed, title="Speed |u|", mode="color")
plt.colorbar(pm)
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.savefig("q5_speed.png", dpi=200)
plt.close()

# --- Write fields to XDMF ---
xdmf = XDMFFile(mesh.mpi_comm(), "q5_soln.xdmf")
xdmf.parameters["flush_output"] = True
xdmf.parameters["functions_share_mesh"] = True
xdmf.write(u_sol, 0.0)
xdmf.write(p_sol, 0.0)
xdmf.close()

print("Saved colour map to q5_speed.png and solution fields to q5_soln.xdmf.")

# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
import json

mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

u_sol.rename("u", "velocity field")
p_sol.rename("p", "pressure field")

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
            "meaning": "velocity field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'u', -1)",
        },
        {
            "symbol": "p_sol",
            "checkpoint_name": "p",
            "meaning": "pressure field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'p', -1)",
        },
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")
