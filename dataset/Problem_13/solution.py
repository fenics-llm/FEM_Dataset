# Copied from Results_ALL-FEM-main/reference solutions/solid/13/solid13.py for dataset Problem_13.
# Original benchmark: Solid Mechanics Problem 13 (Medium).

# q13_neo_hooke_incompressible_plane_strain.py
# Legacy FEniCS (dolfin + mshr)
#
# Geometry: Ω = (0,1.0)×(0,0.20) m minus a circular hole of radius a=0.04 m centred at (0.50, 0.10).
# BCs: left edge fixed; hole has follower pressure P=0.10 MPa; other outer edges traction-free.
# Model: large-deformation, isotropic incompressible Neo-Hookean (mixed (u,p) formulation), plane strain.
# Output: q13_def.png (deformed ×5), q13_vm.png (von Mises), q13_disp.xdmf (displacement).

from __future__ import print_function
import sys
import numpy as np
import dolfin as df
from mshr import Rectangle, Circle, generate_mesh

# -------------------------
# Parameters
# -------------------------
Lx, Ly = 1.0, 0.20
cx, cy, a = 0.50, 0.10, 0.04
P_hole = 0.10e6      # Pa
E = 5.0e6            # Pa
nu = 0.5
mu = E/(2.0*(1.0+nu))  # shear modulus (Pa)
magnify = 5.0
mesh_res = 96         # increase for accuracy if needed

# -------------------------
# Geometry and mesh (mshr CSG)
# -------------------------
domain = Rectangle(df.Point(0.0, 0.0), df.Point(Lx, Ly)) - Circle(df.Point(cx, cy), a, 64)
mesh = generate_mesh(domain, mesh_res)

# -------------------------
# Boundary markers
# -------------------------
left_id, right_id, top_id, bot_id, hole_id = 1, 2, 3, 4, 5
boundaries = df.MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)

class Left(df.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and df.near(x[0], 0.0, 1e-8)
class Right(df.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and df.near(x[0], Lx, 1e-8)
class Top(df.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and df.near(x[1], Ly, 1e-8)
class Bottom(df.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and df.near(x[1], 0.0, 1e-8)
class Hole(df.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and ( (x[0]-cx)**2 + (x[1]-cy)**2 ) <= (a + 1e-6)**2

Left().mark(boundaries, left_id)
Right().mark(boundaries, right_id)
Top().mark(boundaries, top_id)
Bottom().mark(boundaries, bot_id)
Hole().mark(boundaries, hole_id)

ds = df.Measure('ds', domain=mesh, subdomain_data=boundaries)

# -------------------------
# Mixed space (Taylor–Hood)
# -------------------------
V = df.VectorElement("Lagrange", mesh.ufl_cell(), 2)
Q = df.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
W = df.FunctionSpace(mesh, df.MixedElement([V, Q]))

w = df.Function(W)          # unknowns (u, p)
u, p = df.split(w)
v, q = df.TestFunctions(W)

I = df.Identity(2)
F = I + df.grad(u)          # 2D deformation gradient
J2 = df.det(F)              # 2D Jacobian (area ratio)
b = F*F.T                   # left Cauchy-Green

# Incompressible NH (Kirchhoff): τ = -p I + μ b, with plane strain handled separately in σ_vm
sigma = -p*I + mu*b

# Internal virtual work + incompressibility constraint
res_int = df.inner(df.grad(v), sigma)*df.dx + q*(J2 - 1.0)*df.dx

# Follower pressure on hole: −∫ P v · (Cof(F) N) ds  (reference surface)
Cof = J2*df.inv(F).T
N = df.FacetNormal(mesh)
res_ext = - P_hole*df.dot(v, df.dot(Cof, N))*ds(hole_id)

Ftotal = res_int + res_ext
Jtotal = df.derivative(Ftotal, w, df.TrialFunction(W))

# -------------------------
# Boundary conditions: left edge fixed
# -------------------------
bc_left = df.DirichletBC(W.sub(0), df.Constant((0.0, 0.0)), boundaries, left_id)
bcs = [bc_left]

# -------------------------
# Solve nonlinear system
# -------------------------
problem = df.NonlinearVariationalProblem(Ftotal, w, bcs, Jtotal)
solver = df.NonlinearVariationalSolver(problem)
prm = solver.parameters['newton_solver']
prm['absolute_tolerance'] = 1e-9
prm['relative_tolerance'] = 1e-8
prm['maximum_iterations'] = 40
prm['linear_solver'] = 'mumps'
prm['report'] = True
solver.solve()

u_sol, p_sol = w.split(deepcopy=True)

# -------------------------
# Save displacement field (XDMF)
# -------------------------
xdmf = df.XDMFFile(df.MPI.comm_world, "q13_disp.xdmf")
xdmf.parameters["flush_output"] = True
xdmf.parameters["functions_share_mesh"] = True
xdmf.write(u_sol, 0.0)
xdmf.close()

# -------------------------
# Von Mises (plane strain)
# For plane strain with incompressibility: F33 = 1/J2 ⇒ b_zz = F33^2 = 1/J2^2, J3D=1
# σ = -p I + μ b in 3D ⇒ σ_zz = -p + μ b_zz
# σ_vm = sqrt( ( (σx-σy)^2 + (σy-σz)^2 + (σz-σx)^2 + 6 τxy^2 )/2 )
# -------------------------
V1 = df.FunctionSpace(mesh, "Lagrange", 1)
sx = df.project(( -p + mu*b[0,0] ), V1)
sy = df.project(( -p + mu*b[1,1] ), V1)
txy = df.project(( mu*b[0,1] ), V1)  # Cauchy shear
J2_fun = df.project(J2, V1)
bzz = df.project(1.0/(J2_fun*J2_fun), V1)
sz = df.project( -p + mu*bzz, V1 )

vm = df.Function(V1)
vm_expr = df.sqrt(0.5*((sx - sy)**2 + (sy - sz)**2 + (sz - sx)**2 + 6.0*txy**2))
vm.assign(df.project(vm_expr, V1))

# -------------------------
# Plots: deformed mesh (×5) and von Mises on deformed config
# -------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

# Deformed nodes (× magnify)
coords = mesh.coordinates()
u_vert = u_sol.compute_vertex_values(mesh).reshape(2, -1).T
def_coords = coords + magnify*u_vert
cells = mesh.cells()
tri = Triangulation(def_coords[:,0], def_coords[:,1], cells)

# Deformed wireframe
plt.figure()
plt.triplot(tri)
plt.gca().set_aspect('equal', 'box')
plt.title("Deformed configuration (×{:.0f})".format(magnify))
plt.xlabel("x [m]"); plt.ylabel("y [m]")
plt.tight_layout()
plt.savefig("q13_def.png", dpi=200); plt.close()

# Von Mises colour map (on deformed)
vm_vert = vm.compute_vertex_values(mesh)
plt.figure()
pc = plt.tripcolor(tri, vm_vert, shading='gouraud')
plt.gca().set_aspect('equal', 'box')
plt.title("Von Mises equivalent stress [Pa]")
plt.xlabel("x [m]"); plt.ylabel("y [m]")
cb = plt.colorbar(pc); cb.set_label("σ_vm [Pa]")
plt.tight_layout()
plt.savefig("q13_vm.png", dpi=200); plt.close()

print("Done. Files written: q13_def.png, q13_vm.png, q13_disp.xdmf")

# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
import json

mesh_file = df.XDMFFile(df.MPI.comm_world, "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

u_sol.rename("u", "displacement field")
p_sol.rename("p", "pressure-like Lagrange multiplier field")

solution_file = df.XDMFFile(df.MPI.comm_world, "solution.xdmf")
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u_sol, 0.0)
solution_file.write(p_sol, 0.0)
solution_file.close()

checkpoint_file = df.XDMFFile(df.MPI.comm_world, "read_checkpoint.xdmf")
checkpoint_file.write_checkpoint(u_sol, "u", 0.0, df.XDMFFile.Encoding.HDF5, False)
checkpoint_file.write_checkpoint(p_sol, "p", 0.0, df.XDMFFile.Encoding.HDF5, True)
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
        {
            "symbol": "p_sol",
            "checkpoint_name": "p",
            "meaning": "pressure-like Lagrange multiplier field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'p', -1)",
        },
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")
