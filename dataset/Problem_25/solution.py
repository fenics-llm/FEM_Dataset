# Copied from Results_ALL-FEM-main/reference solutions/fluid/9/fluid9.py for dataset Problem_25.
# Original benchmark: Fluid Mechanics Problem 9 (Medium).

# filename: q10_navier_stokes_advection.py
from dolfin import *
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# 1. Mesh and geometry
# ------------------------------------------------------------
L, H = 2.0, 0.20
nx, ny = 200, 20
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), nx, ny, "crossed")

# ------------------------------------------------------------
# 2. Physical parameters
# ------------------------------------------------------------
rho = 1.0          # kg/m^3
mu  = 0.01         # Pa·s
Ubar = 0.1         # m/s
kappa = 1e-3       # m^2/s

# ------------------------------------------------------------
# 3. Boundary definitions
# ------------------------------------------------------------
tol = 1E-10
class Inlet(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0, tol)
class Outlet(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], L, tol)
class Walls(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[1], 0.0, tol) or near(x[1], H, tol))

inlet = Inlet()
outlet = Outlet()
walls  = Walls()

# ------------------------------------------------------------
# 4. Taylor–Hood (P2-P1) mixed space for Navier–Stokes
# ------------------------------------------------------------
Ve = VectorElement("Lagrange", mesh.ufl_cell(), 2)
Pe = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = MixedElement([Ve, Pe])
W  = FunctionSpace(mesh, TH)

# ------------------------------------------------------------
# 5. Boundary conditions for velocity & pressure
# ------------------------------------------------------------
# Inlet parabolic profile
u_inlet_expr = Expression(("6*Ubar*(x[1]/H)*(1.0 - x[1]/H)", "0.0"),
                          degree=2, Ubar=Ubar, H=H)

bcu_inlet = DirichletBC(W.sub(0), u_inlet_expr, inlet)
bcu_walls = DirichletBC(W.sub(0), Constant((0.0, 0.0)), walls)

# Pressure gauge at a point (0,0) to fix nullspace
class PointGauge(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0, tol) and near(x[1], 0.0, tol)
point_gauge = PointGauge()
bcp_point = DirichletBC(W.sub(1), Constant(0.0), point_gauge, method="pointwise")

bcs_NS = [bcu_inlet, bcu_walls, bcp_point]

# ------------------------------------------------------------
# 6. Navier–Stokes variational problem (steady)
# ------------------------------------------------------------
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

# Nonlinear term: (u·∇)u
U = Function(W)               # current iterate (u,p)
(u_, p_) = split(U)

F_NS = (2*mu*inner(sym(grad(u_)), sym(grad(v))) * dx
        + rho*inner(dot(u_, nabla_grad(u_)), v) * dx
        - div(v)*p_ * dx
        + q*div(u_) * dx)

J = derivative(F_NS, U, TrialFunction(W))

# Solve nonlinear problem
solve(F_NS == 0, U, bcs_NS, J=J,
      solver_parameters={"newton_solver":
                         {"relative_tolerance": 1e-6,
                          "absolute_tolerance": 1e-8,
                          "maximum_iterations": 25,
                          "linear_solver": "mumps"}})

# Extract velocity and pressure
u_sol, p_sol = U.split(deepcopy=True)
u, p = u_sol, p_sol

# ------------------------------------------------------------
# 7. Advection–Diffusion for concentration c
# ------------------------------------------------------------
Vc = FunctionSpace(mesh, "Lagrange", 1)
c = TrialFunction(Vc)
v_c = TestFunction(Vc)

# Dirichlet BCs for concentration
bcs_c = [DirichletBC(Vc, Constant(0.0), inlet),
         DirichletBC(Vc, Constant(1.0), outlet)]

a_c = (kappa*dot(grad(c), grad(v_c)) + dot(u, grad(c))*v_c) * dx
L_c = Constant(0.0) * v_c * dx

c_sol = Function(Vc)
solve(a_c == L_c, c_sol, bcs_c,
      solver_parameters={"linear_solver": "mumps"})

# ------------------------------------------------------------
# 8. Output to XDMF
# ------------------------------------------------------------
xdmf = XDMFFile(mesh.mpi_comm(), "q10_solution.xdmf")
xdmf.parameters["functions_share_mesh"] = True
xdmf.parameters["flush_output"] = True
xdmf.write(u, 0.0)
xdmf.write(p, 0.0)
xdmf.write(c_sol, 0.0)
xdmf.close()

# ------------------------------------------------------------
# 9. Plot concentration and save as PNG
# ------------------------------------------------------------
plt.figure(figsize=(8, 2))
p = plot(c_sol, title="Concentration", cmap="viridis")
plt.colorbar(p)
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.tight_layout()
plt.savefig("q10_conc.png", dpi=300)
plt.close()

# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
import json

mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

u_sol.rename("u", "velocity field")
p_sol.rename("p", "pressure field")
c_sol.rename("c", "concentration field")

solution_file = XDMFFile(mesh.mpi_comm(), "solution.xdmf")
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u_sol, 0.0)
solution_file.write(p_sol, 0.0)
solution_file.write(c_sol, 0.0)
solution_file.close()

checkpoint_file = XDMFFile(mesh.mpi_comm(), "read_checkpoint.xdmf")
checkpoint_file.write_checkpoint(u_sol, "u", 0.0, XDMFFile.Encoding.HDF5, False)
checkpoint_file.write_checkpoint(p_sol, "p", 0.0, XDMFFile.Encoding.HDF5, True)
checkpoint_file.write_checkpoint(c_sol, "c", 0.0, XDMFFile.Encoding.HDF5, True)
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
            "symbol": "c_sol",
            "checkpoint_name": "c",
            "meaning": "concentration field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'c', -1)",
        },
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")
