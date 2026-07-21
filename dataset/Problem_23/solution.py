# Copied from Results_ALL-FEM-main/reference solutions/fluid/7/ref_fluid7.py for dataset Problem_23.
# Original benchmark: Fluid Mechanics Problem 7 (Medium).

# -*- coding: utf-8 -*-
#
#  Navier–Stokes flow in a 2-D channel with a circular obstacle
#  (legacy FEniCS – dolfin)
#
from dolfin import *
from mshr import Rectangle, Circle, generate_mesh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json
import numpy as np

# ----------------------------------------------------------------------
# 1. Geometry & mesh
# ----------------------------------------------------------------------
L, H = 2.2, 0.41
R, cx, cy = 0.05, 0.20, 0.20

channel = Rectangle(Point(0.0, 0.0), Point(L, H))
obstacle = Circle(Point(cx, cy), R, 160)
domain = channel - obstacle

# mesh resolution: finer near the obstacle
mesh = generate_mesh(domain, 320)   # 80 ≈ global cell size, adjust if needed

# ----------------------------------------------------------------------
# 2. Boundary markers
# ----------------------------------------------------------------------
class Inlet(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0) and on_boundary

class Outlet(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], L) and on_boundary

class Walls(SubDomain):
    def inside(self, x, on_boundary):
        return (near(x[1], 0.0) or near(x[1], H)) and on_boundary

class CircleBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and ( (x[0]-cx)**2 + (x[1]-cy)**2 < (R+1e-8)**2 )

boundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
Inlet().mark(boundaries, 1)
Outlet().mark(boundaries, 2)
Walls().mark(boundaries, 3)
CircleBoundary().mark(boundaries, 4)

ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# ----------------------------------------------------------------------
# 3. Function spaces (Taylor–Hood)
# ----------------------------------------------------------------------
Ve = VectorElement('Lagrange', mesh.ufl_cell(), 2)
Pe = FiniteElement('Lagrange', mesh.ufl_cell(), 1)
TH = MixedElement([Ve, Pe])
W  = FunctionSpace(mesh, TH)

# ----------------------------------------------------------------------
# 4. Boundary conditions
# ----------------------------------------------------------------------
Umean = 0.2
mu    = 1e-3
rho   = 1.0

# Inlet parabolic profile
inlet_expr = ('6.0*Umean*x[1]*(H - x[1]) / (H*H)', '0.0')
inlet_bc = DirichletBC(W.sub(0), Expression(inlet_expr, Umean=Umean, H=H, degree=2), boundaries, 1)

# No-slip on walls and obstacle
noslip_bc = DirichletBC(W.sub(0), Constant((0.0, 0.0)), boundaries, 3)
circle_bc = DirichletBC(W.sub(0), Constant((0.0, 0.0)), boundaries, 4)

# Pressure reference at a single outlet point (pointwise BC)
# Pressure reference at a single point (pointwise BC)  --- FIXED
class PressurePoint(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], L, DOLFIN_EPS) and near(x[1], H/2.0, DOLFIN_EPS)

p0_bc = DirichletBC(W.sub(1), Constant(0.0), PressurePoint(), method="pointwise")


bcs = [inlet_bc, noslip_bc, circle_bc, p0_bc]

# ----------------------------------------------------------------------
# 5. Variational formulation (steady Navier–Stokes)
# ----------------------------------------------------------------------
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

U = Function(W)          # current solution (u,p)
(u_, p_) = split(U)      # for the nonlinear term

F = ( rho*dot(dot(u_, nabla_grad(u_)), v)*dx
      + 2*mu*inner(sym(grad(u_)), sym(grad(v)))*dx
      - div(v)*p_*dx
      + q*div(u_)*dx )

# Jacobian
J = derivative(F, U, TrialFunction(W))

# ----------------------------------------------------------------------
# 6. Solve (Newton)
# ----------------------------------------------------------------------
solve(F == 0, U, bcs, J=J,
      solver_parameters={'newton_solver':
                         {'relative_tolerance': 1e-6,
                          'absolute_tolerance': 1e-8,
                          'maximum_iterations': 30,
                          'linear_solver': 'mumps'}})

# Split solution
(u_sol, p_sol) = U.split(deepcopy=True)

# ----------------------------------------------------------------------
# 7. Drag force on the circular obstacle
# ----------------------------------------------------------------------
n = FacetNormal(mesh)
sigma = -p_sol*Identity(2) + mu*(grad(u_sol) + grad(u_sol).T)
drag = -assemble(dot(dot(sigma, n), Constant((1.0, 0.0)))*ds(4))

D = 2.0*R
C_D = 2.0*drag/(rho*Umean**2*D)

print('Drag force   = {:.6e} N'.format(drag))
print('Drag coeff.  = {:.6e}'.format(C_D))
with open('drag_coefficient.json', 'w', encoding='utf-8') as f:
    json.dump({
        'drag_force_N': float(drag),
        'drag_coefficient': float(C_D),
        'diameter_m': float(D),
        'mean_inlet_velocity_m_per_s': float(Umean),
    }, f, indent=2)
    f.write('\n')

# ----------------------------------------------------------------------
# 8. Post-processing: speed field PNG
# ----------------------------------------------------------------------
speed = sqrt(dot(u_sol, u_sol))
V_cg = FunctionSpace(mesh, 'CG', 1)
speed_cg = project(speed, V_cg)

plt.figure(figsize=(8, 3))
p = plot(speed_cg, cmap='viridis')
plt.colorbar(p, label='|u| (m/s)')
plt.title('Speed magnitude')
plt.axis('off')
plt.savefig('q7_speed_2.png', dpi=300, bbox_inches='tight')
plt.close()

# ----------------------------------------------------------------------
# 9. Save solution to XDMF
# ----------------------------------------------------------------------
with XDMFFile(mesh.mpi_comm(), 'q7_soln_2.xdmf') as xdmf:
    xdmf.write(u_sol, 0.0)
    xdmf.write(p_sol, 0.0)

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
