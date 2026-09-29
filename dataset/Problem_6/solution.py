# Copied from Results_ALL-FEM-main/reference solutions/solid/6/oss-multi-solid6.py for dataset Problem_6.
# Original benchmark: Solid Mechanics Problem 6 (Easy).

# -*- coding: utf-8 -*-
from dolfin import *
from mshr import Rectangle, Circle, generate_mesh
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# 1. Geometry (rectangle with a semicircular notch)
# ----------------------------------------------------------------------
a = 0.05
domain = Rectangle(Point(0.0, 0.0), Point(1.0, 0.20)) \
         - Circle(Point(0.5, 0.20), a, 64)
mesh = generate_mesh(domain, 128)

# ----------------------------------------------------------------------
# 2. Function space (plane-stress displacement)
# ----------------------------------------------------------------------
V = VectorFunctionSpace(mesh, "Lagrange", 1)

# ----------------------------------------------------------------------
# 3. Plane-stress material parameters
# ----------------------------------------------------------------------
E  = 200e9
nu = 0.30
mu    = E / (2.0 * (1.0 + nu))
lmbda = 2.0 * mu * nu / (1.0 - nu)

def epsilon(u):
    return sym(grad(u))

def sigma(u):
    return lmbda * tr(epsilon(u)) * Identity(2) + 2.0 * mu * epsilon(u)

# ----------------------------------------------------------------------
# 4. Bottom Dirichlet (fixed) boundary condition
# ----------------------------------------------------------------------
tol = 1e-6
class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.0, tol)

bottom = Bottom()
bc_bottom = DirichletBC(V, Constant((0.0, 0.0)), bottom)

# ----------------------------------------------------------------------
# 5. Mark top edge (excluding the notch) for traction
# ----------------------------------------------------------------------
facets = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)

class TopNoNotch(SubDomain):
    def inside(self, x, on_boundary):
        # The semicircular notch arc lies below y = 0.20.  Marking the
        # straight top boundary by height includes the facets adjacent to the
        # two line--arc junctions without applying traction to the arc.
        return on_boundary and near(x[1], 0.20, tol)

top_no_notch = TopNoNotch()
top_no_notch.mark(facets, 1)

ds = Measure("ds", domain=mesh, subdomain_data=facets)

# ----------------------------------------------------------------------
# 6. Variational problem (plane-stress elasticity)
# ----------------------------------------------------------------------
u = TrialFunction(V)
v = TestFunction(V)

traction = Constant((0.0, -10e6))
a = inner(sigma(u), epsilon(v)) * dx
L = dot(traction, v) * ds(1)

# ----------------------------------------------------------------------
# 7. Solve
# ----------------------------------------------------------------------
u_sol = Function(V)
solve(a == L, u_sol, bc_bottom)

# ----------------------------------------------------------------------
# 8. Post-processing: von Mises stress (plane-stress)
# ----------------------------------------------------------------------
W = TensorFunctionSpace(mesh, "Lagrange", 1)
sigma_h = project(sigma(u_sol), W)

sigma_zz = nu * (sigma_h[0, 0] + sigma_h[1, 1])

S = as_tensor(((sigma_h[0, 0], sigma_h[0, 1], 0.0),
               (sigma_h[1, 0], sigma_h[1, 1], 0.0),
               (0.0,           0.0,           sigma_zz)))  # <-- closed

s_dev = S - (1.0/3.0) * tr(S) * Identity(3)
von_mises = sqrt(3.0/2.0 * inner(s_dev, s_dev))

Vv = FunctionSpace(mesh, "Lagrange", 1)
von_mises_proj = project(von_mises, Vv)

# ----------------------------------------------------------------------
# 9. Output: PNG (von Mises) and XDMF (displacement)
# ----------------------------------------------------------------------
plt.figure(figsize=(8, 3))
p = plot(von_mises_proj, cmap="viridis")
plt.colorbar(p, label="von Mises stress (Pa)")
plt.title("Plane-stress von Mises stress")
plt.xlabel("x (m)"); plt.ylabel("y (m)")
plt.tight_layout()
plt.savefig("q6_vm.png", dpi=300)

xdmf = XDMFFile(mesh.mpi_comm(), "q6_disp.xdmf")
xdmf.write(u_sol)
xdmf.close()

# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
import json

mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

u_sol.rename("u", "displacement field")

solution_file = XDMFFile(mesh.mpi_comm(), "solution.xdmf")
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u_sol, 0.0)
solution_file.close()

checkpoint_file = XDMFFile(mesh.mpi_comm(), "read_checkpoint.xdmf")
checkpoint_file.write_checkpoint(u_sol, "u", 0.0, XDMFFile.Encoding.HDF5, False)
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
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")
