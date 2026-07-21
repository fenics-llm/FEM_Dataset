# Copied from Results_ALL-FEM-main/reference solutions/solid/11/solid11.py for dataset Problem_11.
# Original benchmark: Solid Mechanics Problem 11 (Medium).

# q11_nearly_incompressible_plane_strain.py
from __future__ import print_function
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Geometry: rectangle (0,1) x (0,0.20) with a circular hole (r=0.04) at (0.50,0.10)
# -----------------------------------------------------------------------------
rect = Rectangle(Point(0.0, 0.0), Point(1.0, 0.20))
hole = Circle(Point(0.50, 0.10), 0.04)
domain = rect - hole

# Mesh (tune resolution if needed)
mesh = generate_mesh(domain, 140)  # increase for finer results

# -----------------------------------------------------------------------------
# Material (plane strain): E = 5 MPa, nu = 0.49
# Use Herrmann mixed formulation with bulk modulus kappa
# -----------------------------------------------------------------------------
E = 5.0e6          # Pa
nu = 0.49
mu = E/(2.0*(1.0 + nu))                # shear modulus
kappa = E/(3.0*(1.0 - 2.0*nu))         # bulk modulus (3D), used in Herrmann
# (lambda_ = kappa - 2*mu/3 also available if needed)

# -----------------------------------------------------------------------------
# Function spaces: Taylor–Hood (P2 for u, P1 for p)
# -----------------------------------------------------------------------------
Ve = VectorElement("Lagrange", mesh.ufl_cell(), 2)  # displacement
Pe = FiniteElement("Lagrange", mesh.ufl_cell(), 1)  # pressure-like variable
W = FunctionSpace(mesh, MixedElement([Ve, Pe]))

# -----------------------------------------------------------------------------
# Boundary conditions
# Left edge fixed: u = (0,0)
# Right edge: u = (0.001, 0)
# Top/bottom/hole are traction-free (natural)
# -----------------------------------------------------------------------------
tol = 1e-14

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0, tol)

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 1.0, tol)

left = Left()
right = Right()

U0 = Constant((0.0, 0.0))
U1 = Constant((0.001, 0.0))  # prescribed horizontal displacement

bc_left = DirichletBC(W.sub(0), U0, left)   # on displacement subspace
bc_right = DirichletBC(W.sub(0), U1, right)

bcs = [bc_left, bc_right]

# -----------------------------------------------------------------------------
# Variational formulation (Herrmann mixed method)
# Unknowns: (u, p)
# For plane strain we still use 3D mu/kappa; the formulation is robust.
# -----------------------------------------------------------------------------
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

def eps(w):
    return sym(grad(w))

# Bilinear and linear forms:
a = (2.0*mu*inner(eps(u), eps(v)) - p*div(v) - q*div(u) - (1.0/kappa)*p*q)*dx
L = Constant(0.0)*q*dx  # no body forces, tractions are natural BCs

# -----------------------------------------------------------------------------
# Solve
# -----------------------------------------------------------------------------
w = Function(W, name="solution")
solve(a == L, w, bcs,
      solver_parameters={
          "linear_solver": "mumps" if has_linear_algebra_backend("PETSc") else "default",
          "preconditioner": "default"
      })

u_sol, p_sol = w.split(deepcopy=True)
u_sol.rename("u", "displacement")
p_sol.rename("p", "pressure_like")

# -----------------------------------------------------------------------------
# Stress tensor and von Mises (plane strain via Herrmann: sigma = 2 mu eps(u) - p I)
# For plane strain, sigma_zz = -p, tau_xz = tau_yz = 0.
# -----------------------------------------------------------------------------
I = Identity(2)
sigma_2D = 2.0*mu*eps(u_sol) - p_sol*I   # 2x2 in-plane stress

# Components for von Mises in 3D (plane strain extension)
sx = sigma_2D[0, 0]
sy = sigma_2D[1, 1]
txy = sigma_2D[0, 1]
# Out-of-plane normal stress (plane strain): szz = -p
szz = -p_sol

# von Mises (3D) formula:
# vm = sqrt(0.5*((sx - sy)^2 + (sy - szz)^2 + (szz - sx)^2) + 3*(txy^2))
vm_expr = sqrt(0.5*((sx - sy)**2 + (sy - szz)**2 + (szz - sx)**2) + 3.0*(txy**2))

# Project to a discontinuous space for sharper colour maps
DG0 = FunctionSpace(mesh, "DG", 0)
vm = project(vm_expr, DG0)
vm.rename("von_Mises", "von_Mises")

# Also project u_x for plotting
CG1 = FunctionSpace(mesh, "CG", 1)
ux = project(u_sol[0], CG1)
ux.rename("u_x", "u_x")

# -----------------------------------------------------------------------------
# Outputs
#   1) PNG colour maps: von Mises (q11_vm.png) and u_x (q11_ux.png)
#   2) XDMF: displacement field (and also pressure + vm for convenience)
# -----------------------------------------------------------------------------
# 1) PNGs
plt.figure()
p1 = plot(vm, title="von Mises stress (Pa)")
plt.colorbar(p1)
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()
plt.savefig("q11_vm.png", dpi=300)
plt.close()

plt.figure()
p2 = plot(ux, title="Horizontal displacement u_x (m)")
plt.colorbar(p2)
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()
plt.savefig("q11_ux.png", dpi=300)
plt.close()

# 2) XDMF (displacement; plus extras for post-processing)
with XDMFFile(mesh.mpi_comm(), "q11_results.xdmf") as xdmf:
    xdmf.parameters["flush_output"] = True
    xdmf.parameters["functions_share_mesh"] = True
    xdmf.write(mesh)
    xdmf.write(u_sol, 0.0)   # displacement field
    # Optional: also write pressure-like field and von Mises
    xdmf.write(p_sol, 0.0)
    xdmf.write(vm, 0.0)

print("Saved: q11_vm.png, q11_ux.png, q11_results.xdmf")

# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
import json

mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

u_sol.rename("u", "displacement field")
p_sol.rename("p", "pressure-like Lagrange multiplier field")

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
