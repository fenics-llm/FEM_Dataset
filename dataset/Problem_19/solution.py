# Copied from Results_ALL-FEM-main/reference solutions/fluid/3/ref_fluid3.py for dataset Problem_19.
# Original benchmark: Fluid Mechanics Problem 3 (Easy).

# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dolfin import *

# ------------------------------------------------------------
# 0. Settings
# ------------------------------------------------------------
set_log_level(LogLevel.WARNING)

# ------------------------------------------------------------
# 1. Mesh and geometry: unit square with 96 x 96 structured cells
# ------------------------------------------------------------
nx, ny = 96, 96
mesh = UnitSquareMesh(nx, ny)

# ------------------------------------------------------------
# 2. Taylor–Hood mixed finite element (P2 velocity, P1 pressure)
# ------------------------------------------------------------
Ve = VectorElement("Lagrange", mesh.ufl_cell(), 2)  # P2
Pe = FiniteElement("Lagrange", mesh.ufl_cell(), 1)  # P1
W  = FunctionSpace(mesh, MixedElement([Ve, Pe]))

# ------------------------------------------------------------
# 3. Boundary conditions
#    Lid (top): u = (1, 0)
#    Other walls: u = (0, 0)
#    Pressure: fix nullspace by pinning p = 0 at one point
# ------------------------------------------------------------
tol = 1e-14

class Lid(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 1.0, tol)

class Walls(SubDomain):
    def inside(self, x, on_boundary):
        # left OR right OR bottom, but explicitly NOT the top (avoids corner conflict)
        return (on_boundary
                and (near(x[0], 0.0, tol) or near(x[0], 1.0, tol) or near(x[1], 0.0, tol))
                and (not near(x[1], 1.0, tol)))

lid = Lid()
walls = Walls()

u_lid  = Constant((1.0, 0.0))
u_wall = Constant((0.0, 0.0))

bc_u_lid   = DirichletBC(W.sub(0), u_lid,  lid)
bc_u_walls = DirichletBC(W.sub(0), u_wall, walls)

# Pressure pinning at (0,0) to remove constant-nullspace in p
class PinPoint(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0, tol) and near(x[1], 0.0, tol)

bc_p = DirichletBC(W.sub(1), Constant(0.0), PinPoint(), method="pointwise")

bcs = [bc_u_walls, bc_u_lid, bc_p]

# ------------------------------------------------------------
# 4. Variational problem (steady Stokes)
#    -mu * Δu + ∇p = f
#     div(u) = 0
# ------------------------------------------------------------
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

mu = Constant(1.0)          # dynamic viscosity (Pa·s); density not used in steady Stokes
f  = Constant((0.0, 0.0))   # body force

# Standard mixed weak form (as in FEniCS Stokes demos)
a = mu * inner(grad(u), grad(v)) * dx - div(v) * p * dx + q * div(u) * dx
L = inner(f, v) * dx

# ------------------------------------------------------------
# 5. Assemble and solve (direct solver)
# ------------------------------------------------------------
w = Function(W)
solve(a == L, w, bcs, solver_parameters={"linear_solver": "mumps"})

# ------------------------------------------------------------
# 6. Split solution
# ------------------------------------------------------------
u_h, p_h = w.split(deepcopy=True)
u_h.rename("u", "velocity")
p_h.rename("p", "pressure")

# ------------------------------------------------------------
# 7. Post-processing: speed |u| and PNG output
# ------------------------------------------------------------
speed = sqrt(dot(u_h, u_h))
V1 = FunctionSpace(mesh, "CG", 1)
speed_cg = project(speed, V1)

plt.figure(figsize=(6, 5))
h = plot(speed_cg, title="Speed |u|", cmap="viridis")
plt.colorbar(h)
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.savefig("q3_speed.png", dpi=300)
plt.close()

# ------------------------------------------------------------
# 8. XDMF export of velocity and pressure
# ------------------------------------------------------------
xdmf = XDMFFile(mesh.mpi_comm(), "q3_soln.xdmf")
xdmf.parameters["functions_share_mesh"] = True
xdmf.parameters["rewrite_function_mesh"] = False
xdmf.parameters["flush_output"] = True

xdmf.write(u_h, 0.0)
xdmf.write(p_h, 0.0)
xdmf.close()

# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
import json

mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

u_h.rename("u", "velocity field")
p_h.rename("p", "pressure field")

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
            "meaning": "velocity field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'u', -1)",
        },
        {
            "symbol": "p_h",
            "checkpoint_name": "p",
            "meaning": "pressure field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'p', -1)",
        },
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")
