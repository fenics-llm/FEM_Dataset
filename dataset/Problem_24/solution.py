# Copied from Results_ALL-FEM-main/reference solutions/fluid/8/fluid8.py for dataset Problem_24.
# Original benchmark: Fluid Mechanics Problem 8 (Easy).

# file: periodic_channel_bodyforce.py
from __future__ import print_function
from dolfin import *

# -----------------------
# Geometry and mesh
# -----------------------
Lx, Ly = 1.0, 0.20
nx, ny = 128, 32
mesh = RectangleMesh(Point(0.0, 0.0), Point(Lx, Ly), nx, ny, "crossed")  # triangles are robust

# -----------------------
# Periodic boundary in x
# -----------------------
class PeriodicBoundary(SubDomain):
    # Left boundary is "target"
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0) and (not near(x[1], Ly))
    # Map right boundary onto left
    def map(self, x, y):
        y[0] = x[0] - Lx
        y[1] = x[1]

pbc = PeriodicBoundary()

# -----------------------
# Walls (y=0, y=Ly) markers
# -----------------------
walls_id = 1
facets = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)

tol = 1e-10
class WallBottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.0, tol)

class WallTop(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], Ly, tol)

WallBottom().mark(facets, walls_id)
WallTop().mark(facets, walls_id)
ds_ = Measure("ds", domain=mesh, subdomain_data=facets)

# -----------------------
# Parameters
# -----------------------
rho = 1.0          # kg m^-3
mu  = 0.01         # Pa s
G   = 1.0          # N m^-3 (body force in x)

# -----------------------
# Spaces (P2-P1) with periodic constraint
# -----------------------
V = VectorFunctionSpace(mesh, "Lagrange", 2, constrained_domain=pbc)
Q = FunctionSpace(mesh, "Lagrange", 1, constrained_domain=pbc)
mixed_el = MixedElement([V.ufl_element(), Q.ufl_element()])
W = FunctionSpace(mesh, mixed_el, constrained_domain=pbc)

# -----------------------
# Boundary conditions
# -----------------------
# No-slip on walls for velocity
bc_walls = DirichletBC(W.sub(0), Constant((0.0, 0.0)), facets, walls_id)

bcs = [bc_walls]

# -----------------------
# Variational formulation (steady Navier–Stokes + body force)
# -----------------------
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

w = Function(W)
u_, p_ = split(w)

def eps(u):
    return sym(grad(u))

f = Constant((G, 0.0))
n = FacetNormal(mesh)
I = Identity(2)

# Nonlinear residual
F = (
    rho*inner(grad(u_)*u_, v)*dx        # convection
  + 2.0*mu*inner(eps(u_), eps(v))*dx    # viscous
  - p_*div(v)*dx                        # pressure-velocity coupling
  + q*div(u_)*dx                        # incompressibility
  - inner(f, v)*dx                      # body force
)

# Jacobian for Newton
J = derivative(F, w, TrialFunction(W))

# -----------------------
# Solve (Newton)
# -----------------------
problem = NonlinearVariationalProblem(F, w, bcs, J)
solver = NonlinearVariationalSolver(problem)
prm = solver.parameters
prm["nonlinear_solver"] = "newton"
prm["newton_solver"]["relative_tolerance"] = 1e-9
prm["newton_solver"]["absolute_tolerance"] = 1e-11
prm["newton_solver"]["maximum_iterations"] = 50
prm["newton_solver"]["linear_solver"] = "mumps"
solver.solve()

u_sol, p_sol = w.split(deepcopy=True)

# Pressure is determined up to an additive constant. Normalize the stored
# primitive pressure field to zero domain mean instead of pinning a point value.
domain_area = assemble(Constant(1.0)*dx(domain=mesh))
p_mean = assemble(p_sol*dx(domain=mesh)) / domain_area
p_sol.vector()[:] -= p_mean
p_sol.vector().apply("insert")

# -----------------------
# Output: XDMF
# -----------------------
xdmf = XDMFFile(mesh.mpi_comm(), "q9_soln.xdmf")
xdmf.parameters["flush_output"] = True
xdmf.parameters["functions_share_mesh"] = True
u_sol.rename("u", "velocity")
p_sol.rename("p", "pressure")
xdmf.write(u_sol, 0.0)
xdmf.write(p_sol, 0.0)
xdmf.close()

# Quick sanity prints
Um = assemble(u_sol[0]*dx)/domain_area  # domain-average u_x
mean_pressure = assemble(p_sol*dx(domain=mesh)) / domain_area
print("Mean streamwise velocity <u_x> =", Um)
print("Mean pressure <p> after normalization =", mean_pressure)
print("Saved (u,p) to q9_soln.xdmf")
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
