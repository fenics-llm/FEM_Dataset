# Copied from Results_ALL-FEM-main/reference solutions/solid/2/solid2.py for dataset Problem_2.
# Original benchmark: Solid Mechanics Problem 2 (Easy).

# q2_plane_stress_top_traction.py
# Legacy FEniCS (dolfin) script: plane-stress elasticity with top-edge traction
from dolfin import *
import numpy as np

# ----------------------------------------------------------------------
# Mesh
# ----------------------------------------------------------------------
Lx, Ly = 1.0, 0.20
nx, ny = 40, 8
mesh = RectangleMesh(Point(0.0, 0.0), Point(Lx, Ly), nx, ny)

# ----------------------------------------------------------------------
# Function space
# ----------------------------------------------------------------------
V = VectorFunctionSpace(mesh, "CG", 1)

# ----------------------------------------------------------------------
# Material (plane stress)
# E = Young's modulus, nu = Poisson's ratio
# mu = E/(2*(1+nu))
# lambda_ps = 2*mu*nu/(1 - nu)   (effective Lamé parameter for plane stress)
# ----------------------------------------------------------------------
E  = 200e9
nu = 0.30
mu = E/(2.0*(1.0 + nu))
lambda_ps = 2.0*mu*nu/(1.0 - nu)

def eps(u):
    return sym(grad(u))

def sigma(u):
    return 2.0*mu*eps(u) + lambda_ps*tr(eps(u))*Identity(2)

# ----------------------------------------------------------------------
# Boundary conditions
# Left edge fixed: u = (0, 0)
# Top edge traction: t = (0, -2000) N/m
# Right and bottom edges: traction-free (natural)
# ----------------------------------------------------------------------
tol = 1e-14

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0, tol)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], Ly, tol)

left = Left()
top = Top()

bcs = [DirichletBC(V, Constant((0.0, 0.0)), left)]

# Mark facets to apply traction only on the top edge
facets = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
top.mark(facets, 1)
ds_ = Measure("ds", domain=mesh, subdomain_data=facets)

# ----------------------------------------------------------------------
# Variational problem
# ----------------------------------------------------------------------
u = TrialFunction(V)
v = TestFunction(V)

f = Constant((0.0, 0.0))                 # no body force
t_top = Constant((0.0, -2000.0))         # N/m

a = inner(sigma(u), eps(v)) * dx
L = dot(f, v) * dx + dot(t_top, v) * ds_(1)   # only on marked Top edge

# ----------------------------------------------------------------------
# Solve
# ----------------------------------------------------------------------
uh = Function(V, name="displacement")
solve(a == L, uh, bcs)

# ----------------------------------------------------------------------
# Save results: XDMF of displacement field
# ----------------------------------------------------------------------
with XDMFFile(MPI.comm_world, "q2_displacement.xdmf") as xdmf:
    xdmf.write(uh)

# ----------------------------------------------------------------------
# Save colour map of vertical displacement u_y to q2_uy.png
# ----------------------------------------------------------------------
try:
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    # Extract u_y (component 1)
    uh_split = uh.split(deepcopy=True)
    uy = uh_split[1]
    vals = uy.compute_vertex_values(mesh)

    coords = mesh.coordinates()
    cells = mesh.cells()

    triang = mtri.Triangulation(coords[:, 0], coords[:, 1], cells)
    plt.figure()
    im = plt.tripcolor(triang, vals, shading="gouraud")
    plt.colorbar(im, label="u_y [m]")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title("Vertical displacement u_y")
    plt.tight_layout()
    plt.savefig("q2_uy.png", dpi=300)
    plt.close()
except Exception as e:
    # Matplotlib may not be available in some environments; continue silently
    print("Plotting skipped:", e)

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
