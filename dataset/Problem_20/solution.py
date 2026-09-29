# Copied from Results_ALL-FEM-main/reference solutions/fluid/4/fluid4.py for dataset Problem_20.
# Original benchmark: Fluid Mechanics Problem 4 (Easy).

# q4_navier_stokes_channel.py
# Legacy FEniCS (dolfin) steady incompressible Navier–Stokes with Taylor–Hood (P2–P1)

from dolfin import *
import matplotlib.pyplot as plt

# --- Geometry & mesh ---
L = 2.0
H = 0.20
nx, ny = 160, 16
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), nx, ny)

# --- Parameters ---
rho = Constant(1.0)     # kg m^-3
mu  = Constant(0.01)    # Pa·s
Ubar = 2.5              # m s^-1 (mean inflow speed)

# --- Function spaces: Taylor–Hood P2–P1 ---
Ve = VectorElement("Lagrange", mesh.ufl_cell(), 2)  # velocity
Pe = FiniteElement("Lagrange", mesh.ufl_cell(), 1)  # pressure
W  = FunctionSpace(mesh, MixedElement([Ve, Pe]))

# --- Boundary markers ---
tol = DOLFIN_EPS

class Inlet(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0, tol)

class Outlet(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], L, tol)

class Walls(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[1], 0.0, tol) or near(x[1], H, tol))

inlet  = Inlet()
outlet = Outlet()
walls  = Walls()

# --- Inlet parabolic profile: ux(y) = 6 Ubar (y/H)(1 - y/H), uy = 0 ---
inlet_expr = Expression(("6.0*Ubar*(x[1]/H)*(1.0 - x[1]/H)", "0.0"),
                        degree=2, Ubar=Ubar, H=H)

# --- Dirichlet BCs on velocity ---
bcs = [
    DirichletBC(W.sub(0), inlet_expr, inlet),                 # inflow velocity
    DirichletBC(W.sub(0), Constant((0.0, 0.0)), walls)        # no-slip walls
]
# No Dirichlet at outlet (do-nothing traction-free is natural)

# The prescribed outlet traction fixes the pressure level, so no pressure
# pin or mean-zero constraint is required.

# --- Unknowns, tests ---
w  = Function(W)               # current iterate (u, p)
u, p = split(w)
v, q = TestFunctions(W)

# --- Nonlinear form for steady Navier–Stokes ---
# F = rho (u · grad) u · v + 2 mu sym(grad u) : sym(grad v) - p div v + q div u
def eps(u):
    return sym(grad(u))

F = (rho*inner(grad(u)*u, v)                             # convection
     + 2.0*mu*inner(eps(u), eps(v))                      # viscous
     - p*div(v) + q*div(u)) * dx                         # pressure/incompressibility

# Note: Traction-free at outlet is the natural boundary condition for this weak form.

# --- Solve with Newton ---
dw = TrialFunction(W)
J = derivative(F, w, dw)

problem = NonlinearVariationalProblem(F, w, bcs, J)
solver  = NonlinearVariationalSolver(problem)

prm = solver.parameters
prm["nonlinear_solver"] = "newton"
prm["newton_solver"]["maximum_iterations"] = 30
prm["newton_solver"]["relative_tolerance"] = 1e-8
prm["newton_solver"]["absolute_tolerance"] = 1e-10
prm["newton_solver"]["report"] = True
prm["newton_solver"]["linear_solver"] = "mumps"   # robust direct solver (PETSc/MUMPS)

# Provide a reasonable initial guess: zero velocity, zero pressure is fine here
assign(w.sub(0), project(Constant((0.0, 0.0)), W.sub(0).collapse()))
assign(w.sub(1), project(Constant(0.0), W.sub(1).collapse()))

solver.solve()

u_sol, p_sol = w.split(True)
u_sol.rename("u", "velocity")
p_sol.rename("p", "pressure")

# --- Post-process: colour map of u_x ---
V1 = FunctionSpace(mesh, "Lagrange", 1)
ux = project(u_sol.sub(0), V1)
ux.rename("ux", "u_x")

plt.figure()
pm = plot(ux, title="Streamwise velocity u_x", mode="color")
plt.colorbar(pm)
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.savefig("q4_ux.png", dpi=200)
plt.close()

# --- Write solution to XDMF ---
xdmf = XDMFFile(mesh.mpi_comm(), "q4_soln.xdmf")
xdmf.parameters["flush_output"] = True
xdmf.parameters["functions_share_mesh"] = True
xdmf.write(u_sol, 0.0)
xdmf.write(p_sol, 0.0)
xdmf.close()

print("Saved colour map to q4_ux.png and fields to q4_soln.xdmf.")

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
