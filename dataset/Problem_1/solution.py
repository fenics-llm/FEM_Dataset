# Copied from Results_ALL-FEM-main/reference solutions/solid/1/ref_solid1.py for dataset Problem_1.
# Original benchmark: Solid Mechanics Problem 1 (Easy).

# -*- coding: utf-8 -*-
"""
Plane‑stress linear elasticity on a rectangular plate.
Domain: (0,1) × (0,0.20) m
Mesh: 20 × 4 Q1 cells
Material: E = 200 GPa, ν = 0.30 (plane stress)
BCs: left fixed, right u_x = 0.001 m, top/bottom traction‑free
Outputs: ux colormap (q1_ux.png) and full displacement (displacement.xdmf)
"""

import matplotlib
matplotlib.use('Agg')               # headless backend
import matplotlib.pyplot as plt

from dolfin import (Point, RectangleMesh, VectorFunctionSpace,
                    Constant, DirichletBC, TestFunction,
                    TrialFunction, Function, inner, sym,
                    grad, dx, solve, XDMFFile, plot,
                    SubDomain, near,
                    tr, Identity, DOLFIN_EPS, dot)   # added dot

# --------------------------------------------------------------
# 1. Mesh
mesh = RectangleMesh(Point(0.0, 0.0), Point(1.0, 0.20), 40, 8, "crossed")

# --------------------------------------------------------------
# 2. Function space (Q1 vector)
V = VectorFunctionSpace(mesh, "Lagrange", 1)

# --------------------------------------------------------------
# 3. Material parameters (plane stress)
E  = 200e9          # Pa
nu = 0.30
mu = E / (2.0 * (1.0 + nu))
lmbda = E * nu / (1.0 - nu**2)   # plane‑stress λ

# --------------------------------------------------------------
# 4. Boundary definitions
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0, DOLFIN_EPS)

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 1.0, DOLFIN_EPS)

left  = Left()
right = Right()

# --------------------------------------------------------------
# 5. Dirichlet BCs
zero_disp = Constant((0.0, 0.0))
right_disp = Constant((0.001, 0.0))

bc_left  = DirichletBC(V, zero_disp, left)
bc_right = DirichletBC(V, right_disp, right)

bcs = [bc_left, bc_right]

# --------------------------------------------------------------
# 6. Variational problem (plane stress)
u = TrialFunction(V)
v = TestFunction(V)

def epsilon(w):
    return sym(grad(w))

def sigma(w):
    eps = epsilon(w)
    return lmbda * tr(eps) * Identity(2) + 2.0 * mu * eps

a = inner(sigma(u), epsilon(v)) * dx
L = dot(Constant((0.0, 0.0)), v) * dx   # zero right‑hand side

# --------------------------------------------------------------
# 7. Solve
u_sol = Function(V)
solve(a == L, u_sol, bcs)

# --------------------------------------------------------------
# 8. Post‑processing
# 8.1 Save ux colormap
plt.figure(figsize=(6, 2.5))
p = plot(u_sol.sub(0), title=r"$u_x$ (m)", mode="color")
plt.colorbar(p)
plt.tight_layout()
plt.savefig("q1_ux.png", dpi=300)

# 8.2 Save full displacement in XDMF
with XDMFFile(mesh.mpi_comm(), "displacement.xdmf") as xdmf:
    xdmf.write(u_sol)

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
