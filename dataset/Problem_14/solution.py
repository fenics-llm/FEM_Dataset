# Copied from Results_ALL-FEM-main/reference solutions/solid/14/q14_neo_hookean.py for dataset Problem_14.
# Original benchmark: Solid Mechanics Problem 14 (Hard).

# filename: q14_neo_hookean.py
from dolfin import *
from mshr import Rectangle, Circle, generate_mesh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --------------------------------------------------------------
# Geometry and mesh
# --------------------------------------------------------------
L, H = 1.0, 0.20
a = 0.04
hole1 = Circle(Point(0.40, 0.10), a, 64)
hole2 = Circle(Point(0.60, 0.10), a, 64)
domain = Rectangle(Point(0.0, 0.0), Point(L, H)) - hole1 - hole2
mesh = generate_mesh(domain, 80)   # increase for finer resolution

# --------------------------------------------------------------
# Function spaces (Taylor–Hood)
# --------------------------------------------------------------
V_el = VectorElement("Lagrange", mesh.ufl_cell(), 2)
P_el = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH   = MixedElement([V_el, P_el])
W    = FunctionSpace(mesh, TH)

# --------------------------------------------------------------
# Boundary markers
# --------------------------------------------------------------
tol = 1E-6
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)

class Left(SubDomain):
    def inside(self, x, on):
        return on and near(x[0], 0.0, tol)
class Right(SubDomain):
    def inside(self, x, on):
        return on and near(x[0], L, tol)
class Top(SubDomain):
    def inside(self, x, on):
        return on and near(x[1], H, tol)
class Bottom(SubDomain):
    def inside(self, x, on):
        return on and near(x[1], 0.0, tol)
class Hole1(SubDomain):
    def inside(self, x, on):
        return on and ( (x[0]-0.40)**2 + (x[1]-0.10)**2 < (a+tol)**2 )
class Hole2(SubDomain):
    def inside(self, x, on):
        return on and ( (x[0]-0.60)**2 + (x[1]-0.10)**2 < (a+tol)**2 )

Left().mark(boundaries, 1)
Right().mark(boundaries, 2)
Top().mark(boundaries, 3)
Bottom().mark(boundaries, 4)
Hole1().mark(boundaries, 5)
Hole2().mark(boundaries, 6)

ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# --------------------------------------------------------------
# Material parameters (Neo‑Hookean, quasi‑incompressible)
# --------------------------------------------------------------
E  = 5.0e6          # Pa
nu = 0.49
mu    = E/(2.0*(1.0+nu))
# kappa = E/(3.0*(1.0-2.0*nu))   # not needed explicitly

# --------------------------------------------------------------
# Mixed variational problem
# --------------------------------------------------------------
w = Function(W)               # (u,p)
(u, p) = split(w)
(v, q) = TestFunctions(W)

I = Identity(2)
F = I + grad(u)               # deformation gradient
J = det(F)

# First Piola‑Kirchhoff stress for incompressible Neo‑Hookean
P = -p*J*inv(F).T + mu*F

# External follower pressure on holes
P_hole = 0.10e6                # Pa
n0 = FacetNormal(mesh)        # reference outward normal

# Traction vector (follower pressure)
traction = -P_hole * J * inv(F).T * n0

# Weak form
R = inner(P, grad(v))*dx + q*(J - 1.0)*dx \
    + dot(traction, v)*ds(5) \
    + dot(traction, v)*ds(6)

# --------------------------------------------------------------
# Boundary conditions
# --------------------------------------------------------------
# Left edge: u = (0,0)
bc_left = DirichletBC(W.sub(0), Constant((0.0, 0.0)), boundaries, 1)

# Right edge: prescribed displacement (+0.060,0)
bc_right = DirichletBC(W.sub(0), Expression(("0.060","0.0"), degree=1), boundaries, 2)

# Pressure gauge (fix p at a point to avoid nullspace)
class PointGauge(SubDomain):
    def inside(self, x, on):
        return near(x[0], 0.0, tol) and near(x[1], 0.0, tol)
gauge = PointGauge()
bc_pressure = DirichletBC(W.sub(1), Constant(0.0), gauge, method='pointwise')

bcs = [bc_left, bc_right, bc_pressure]

# --------------------------------------------------------------
# Newton solver
# --------------------------------------------------------------
J_form = derivative(R, w)
problem = NonlinearVariationalProblem(R, w, bcs, J_form)
solver  = NonlinearVariationalSolver(problem)
prm = solver.parameters
prm['newton_solver']['relative_tolerance'] = 1e-6
prm['newton_solver']['absolute_tolerance'] = 1e-8
prm['newton_solver']['maximum_iterations'] = 25
prm['newton_solver']['linear_solver'] = 'mumps'
solver.solve()

# --------------------------------------------------------------
# Post‑processing
# --------------------------------------------------------------
(u_h, p_h) = w.split(deepcopy=True)

# Cauchy stress
sigma = (1.0/J)*(P*F.T)

# von Mises stress
s = sigma - (1.0/3)*tr(sigma)*I
vonMises = sqrt(3.0/2.0*inner(s, s))

Vsig = FunctionSpace(mesh, "P", 1)
vonMises_h = project(vonMises, Vsig)

# Save displacement (XDMF)
xdmf_disp = XDMFFile(mesh.mpi_comm(), "q14_disp.xdmf")
xdmf_disp.write(u_h)

# Save pressure image
plt.figure(figsize=(6,4))
p_plot = plot(p_h, title="Pressure (Pa)", cmap="viridis")
plt.colorbar(p_plot)
plt.axis('off')
plt.tight_layout()
plt.savefig("q14_p.png", dpi=300)

# Save von Mises image
plt.figure(figsize=(6,4))
vm_plot = plot(vonMises_h, title="von Mises Stress (Pa)", cmap="viridis")
plt.colorbar(vm_plot)
plt.axis('off')
plt.tight_layout()
plt.savefig("q14_vm.png", dpi=300)

print("Simulation completed.")
print("Displacement saved to q14_disp.xdmf")
print("Pressure image saved to q14_p.png")
print("von Mises image saved to q14_vm.png")

# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
import json

mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

u_h.rename("u", "displacement field")
p_h.rename("p", "pressure-like Lagrange multiplier field")

solution_file = XDMFFile(mesh.mpi_comm(), "solution.xdmf")
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u_h, 0.0)
solution_file.write(p_h, 0.0)
solution_file.close()

checkpoint_file = XDMFFile(mesh.mpi_comm(), "read_checkpoint.xdmf")
checkpoint_file.write_checkpoint(u_h, "u", 0.0, XDMFFile.Encoding.HDF5, False)
checkpoint_file.write_checkpoint(p_h, "p", 0.0, XDMFFile.Encoding.HDF5, True)
checkpoint_file.close()

read_checkpoint_metadata = {
    "format": "FEniCS XDMFFile.write_checkpoint",
    "time": 0.0,
    "fields": [
        {
            "symbol": "u_h",
            "checkpoint_name": "u",
            "meaning": "displacement field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'u', -1)",
        },
        {
            "symbol": "p_h",
            "checkpoint_name": "p",
            "meaning": "pressure-like Lagrange multiplier field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'p', -1)",
        },
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")
