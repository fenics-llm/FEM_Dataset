# Copied from Results_ALL-FEM-main/reference solutions/solid/9/solid9.py for dataset Problem_9.
# Original benchmark: Solid Mechanics Problem 9 (Medium).

# q9_plane_stress_grading.py
from dolfin import *
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Mesh (structured 100 x 20 over (0,1) x (0,0.20))
# -----------------------------------------------------------------------------
Lx, Ly = 1.0, 0.20
nx, ny = 100, 20
mesh = RectangleMesh(Point(0.0, 0.0), Point(Lx, Ly), nx, ny, "right")

# -----------------------------------------------------------------------------
# Function spaces
# -----------------------------------------------------------------------------
V = VectorFunctionSpace(mesh, "CG", 1)   # for displacement
Q = FunctionSpace(mesh, "CG", 1)         # for scalar post-processing

# -----------------------------------------------------------------------------
# Material: E(y) graded, plane-stress with constant nu
# -----------------------------------------------------------------------------
nu = Constant(0.30)

# Young's modulus varying with y: E = 100e9 + 100e9*(y/0.20)
E = Expression("1.0e11 + 1.0e11*(x[1]/0.20)", degree=1)

# 3D Lamé parameters (functions of E, nu)
mu = E / (2.0*(1.0 + nu))
lmbda = E*nu / ((1.0 + nu)*(1.0 - 2.0*nu))

# Effective plane-stress lambda (2D constitutive reduction)
lambda_ps = 2.0*mu*lmbda / (lmbda + 2.0*mu)

def eps(u):
    return sym(grad(u))

def sigma(u):
    return 2.0*mu*eps(u) + lambda_ps*tr(eps(u))*Identity(2)

# -----------------------------------------------------------------------------
# Boundary markers
# -----------------------------------------------------------------------------
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0, DOLFIN_EPS)

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], Lx, DOLFIN_EPS)

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
Left().mark(boundaries, 1)
Right().mark(boundaries, 2)

ds = Measure("ds", domain=mesh, subdomain_data=boundaries)

# -----------------------------------------------------------------------------
# Boundary conditions and traction
# -----------------------------------------------------------------------------
zero = Constant((0.0, 0.0))
bc_left = DirichletBC(V, zero, boundaries, 1)

# Uniform traction on the right edge (per unit thickness): (2e6, 0) N/m = Pa
t_right = Constant((2.0e6, 0.0))

# -----------------------------------------------------------------------------
# Variational problem (linear elasticity, plane stress)
# -----------------------------------------------------------------------------
u = TrialFunction(V)
v = TestFunction(V)

a = inner(sigma(u), eps(v))*dx
L = dot(t_right, v)*ds(2)   # Neumann load only on the right boundary

# Solve
uh = Function(V, name="displacement")
solve(a == L, uh, bc_left)

# -----------------------------------------------------------------------------
# Post-processing: displacement magnitude and PNG colour map
# -----------------------------------------------------------------------------
umag = project(sqrt(dot(uh, uh)), Q, solver_type="cg", preconditioner_type="ilu")

plt.figure()
c = plot(umag, title="|u| (displacement magnitude)")
plt.colorbar(c)
plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()
plt.savefig("q9_disp.png", dpi=300)
plt.close()

# -----------------------------------------------------------------------------
# Output: XDMF (field) + also write magnitude if you wish to inspect it
# -----------------------------------------------------------------------------
# In legacy FEniCS, pass the communicator explicitly:
xf = XDMFFile(MPI.comm_world, "displacement.xdmf")
xf.parameters["flush_output"] = True
xf.parameters["functions_share_mesh"] = True
xf.write(uh)

# Optional: write magnitude as a separate scalar field
xf_mag = XDMFFile(MPI.comm_world, "displacement_magnitude.xdmf")
xf_mag.parameters["flush_output"] = True
xf_mag.parameters["functions_share_mesh"] = True
umag.rename("umag", "umag")
xf_mag.write(umag)

print("Saved: q9_disp.png, displacement.xdmf (and displacement_magnitude.xdmf)")

# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
import json

mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

uh.rename("u", "displacement field")

solution_file = XDMFFile(mesh.mpi_comm(), "solution.xdmf")
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(uh, 0.0)
solution_file.close()

checkpoint_file = XDMFFile(mesh.mpi_comm(), "read_checkpoint.xdmf")
checkpoint_file.write_checkpoint(uh, "u", 0.0, XDMFFile.Encoding.HDF5, False)
checkpoint_file.close()

read_checkpoint_metadata = {
    "format": "FEniCS XDMFFile.write_checkpoint",
    "time": 0.0,
    "fields": [
        {
            "symbol": "uh",
            "checkpoint_name": "u",
            "meaning": "displacement field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'u', -1)",
        },
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")
