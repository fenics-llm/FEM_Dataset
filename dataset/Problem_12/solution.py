# Copied from Results_ALL-FEM-main/reference solutions/solid/12/solid12.py for dataset Problem_12.
# Original benchmark: Solid Mechanics Problem 12 (Medium).

# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dolfin import *
from mshr import Rectangle, Circle, generate_mesh

# ----------------------------------------------------------------------
# 1. Geometry and mesh
# ----------------------------------------------------------------------
Lx, Ly = 1.0, 0.20
a = 0.04
hole_center = Point(0.50, 0.10)

domain = Rectangle(Point(0.0, 0.0), Point(Lx, Ly)) - Circle(hole_center, a)
mesh = generate_mesh(domain, 80)   # 80 ≈ mesh resolution (adjust if needed)

# ----------------------------------------------------------------------
# 2. Boundary markers
# ----------------------------------------------------------------------
class Left(SubDomain):
    def inside(self, x, on):
        return near(x[0], 0.0) and on

class Right(SubDomain):
    def inside(self, x, on):
        return near(x[0], Lx) and on

class Top(SubDomain):
    def inside(self, x, on):
        return near(x[1], Ly) and on

class Bottom(SubDomain):
    def inside(self, x, on):
        return near(x[1], 0.0) and on

class Hole(SubDomain):
    def inside(self, x, on):
        return ( (x[0]-hole_center.x())**2 + (x[1]-hole_center.y())**2 < (a+1e-8)**2 ) and on

boundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
Left().mark(boundaries, 1)
Right().mark(boundaries, 2)
Top().mark(boundaries, 3)
Bottom().mark(boundaries, 4)
Hole().mark(boundaries, 5)

ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# ----------------------------------------------------------------------
# 3. Function spaces (Taylor–Hood)
# ----------------------------------------------------------------------
Ve = VectorElement('Lagrange', mesh.ufl_cell(), 2)
Pe = FiniteElement('Lagrange', mesh.ufl_cell(), 1)
TH = MixedElement([Ve, Pe])
W = FunctionSpace(mesh, TH)

# ----------------------------------------------------------------------
# 4. Trial / test functions, solution
# ----------------------------------------------------------------------
w = Function(W)               # (u,p)
(u, p) = split(w)
(v, q) = TestFunctions(W)

# ----------------------------------------------------------------------
# 5. Material parameters (plane strain)
# ----------------------------------------------------------------------
E = 5.0e6
nu = 0.5
mu = E / (2.0 * (1.0 + nu))

# ----------------------------------------------------------------------
# 6. Kinematics
# ----------------------------------------------------------------------
I = Identity(2)
F = I + grad(u)
J = det(F)

# ----------------------------------------------------------------------
# 7. Strain energy and Cauchy stress (incompressible neo-Hookean)
# ----------------------------------------------------------------------
psi = mu/2.0 * (inner(F, F) - 2.0) - p * (J - 1.0)
sigma = mu * F * F.T - p * I

# ----------------------------------------------------------------------
# 8. Weak form (residual and Jacobian)
# ----------------------------------------------------------------------
P = mu*F - p*inv(F).T
R = inner(P, grad(v)) * dx - q * (J - 1.0) * dx
Jac = derivative(R, w, TrialFunction(W))

# ----------------------------------------------------------------------
# 9. Boundary conditions
# ----------------------------------------------------------------------
# Left edge: u = (0,0)
bc_left = DirichletBC(W.sub(0), Constant((0.0, 0.0)), boundaries, 1)

# Right edge: u = (0.060, 0.0)
bc_right = DirichletBC(W.sub(0), Constant((0.060, 0.0)), boundaries, 2)

bcs = [bc_left, bc_right]

# ----------------------------------------------------------------------
# 10. Solve nonlinear problem (Newton)
# ----------------------------------------------------------------------
solve(R == 0, w, bcs, J=Jac,
      solver_parameters={'newton_solver':
                         {'relative_tolerance': 1e-6,
                          'absolute_tolerance': 1e-8,
                          'maximum_iterations': 30,
                          'linear_solver': 'mumps'}})

# ----------------------------------------------------------------------
# 11. Extract fields
# ----------------------------------------------------------------------
(u_h, p_h) = w.split(deepcopy=True)

# ----------------------------------------------------------------------
# 12. Von Mises stress (plane strain)
# ----------------------------------------------------------------------
F_h = I + grad(u_h)
J_h = det(F_h)
F3 = as_matrix([[F_h[0,0], F_h[0,1], 0.0],
                [F_h[1,0], F_h[1,1], 0.0],
                [0.0,       0.0,    1.0/J_h]])
I3 = Identity(3)
B3 = F3*F3.T
sigma3 = mu*B3 - p_h*I3
s3 = sigma3 - (tr(sigma3)/3.0)*I3
von_mises = sqrt(3.0/2.0*inner(s3, s3))

Vsig = FunctionSpace(mesh, 'P', 1)
von_mises_h = project(von_mises, Vsig)

# ----------------------------------------------------------------------
# 13. Save results
# ----------------------------------------------------------------------
# XDMF for displacement
with XDMFFile(mesh.mpi_comm(), "q12_displacement.xdmf") as xdmf:
    xdmf.write(u_h, 0.0)

# PNG for pressure
plt.figure()
p_plot = plot(p_h, title='Hydrostatic pressure (Pa)', cmap='viridis')
plt.colorbar(p_plot)
plt.savefig("q12_p.png", dpi=300, bbox_inches='tight')
plt.close()

# PNG for von Mises stress
plt.figure()
vm_plot = plot(von_mises_h, title='Von Mises stress (Pa)', cmap='viridis')
plt.colorbar(vm_plot)
plt.savefig("q12_vm.png", dpi=300, bbox_inches='tight')
plt.close()

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
