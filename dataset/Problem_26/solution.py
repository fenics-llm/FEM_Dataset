# Copied from Results_ALL-FEM-main/reference solutions/fluid/10/ref_fluid10.py for dataset Problem_26.
# Original benchmark: Fluid Mechanics Problem 10 (Medium).

# -*- coding: utf-8 -*-
#
# 2-D steady Boussinesq cavity (legacy FEniCS)
#
#  Ω = [0,1]×[0,1]  (unit square)
#  ρ = 1.0 kg/m³
#  μ = 1.5e-5 Pa·s
#  α = 2.1e-5 m²/s
#  gβ = 3.15e-5 m·s⁻²·K⁻¹
#  T_ref = 0.5 K
#
#  Left wall  : T = 1, u = (0,0)
#  Right wall : T = 0, u = (0,0)
#  Top/Bottom : ∂T/∂n = 0, u = (0,0)
#
#  Output:
#    - q11_T.png          : temperature colormap
#    - q11_solution.xdmf  : (u,p,T)
#    - average Nusselt number on left wall printed to stdout
#
# --------------------------------------------------------------
from dolfin import *
import matplotlib.pyplot as plt
import json
from dolfin import Point   # for pressure pinning

# --------------------------------------------------------------
# 1. Mesh
# --------------------------------------------------------------
mesh = UnitSquareMesh(128, 128)          # refine as needed

# --------------------------------------------------------------
# 2. Function spaces (Taylor–Hood for (u,p) + CG1 for T)
# --------------------------------------------------------------
V_el = VectorElement("Lagrange", mesh.ufl_cell(), 2)   # velocity
P_el = FiniteElement("Lagrange", mesh.ufl_cell(), 1)   # pressure
U_el = FiniteElement("Lagrange", mesh.ufl_cell(), 1)   # temperature
TH_el = MixedElement([V_el, P_el, U_el])
W = FunctionSpace(mesh, TH_el)

# --------------------------------------------------------------
# 3. Boundary markers
# --------------------------------------------------------------
tol = 1E-14
class Left(SubDomain):
    def inside(self, x, on_boundary): return on_boundary and near(x[0], 0.0, tol)
class Right(SubDomain):
    def inside(self, x, on_boundary): return on_boundary and near(x[0], 1.0, tol)
class Bottom(SubDomain):
    def inside(self, x, on_boundary): return on_boundary and near(x[1], 0.0, tol)
class Top(SubDomain):
    def inside(self, x, on_boundary): return on_boundary and near(x[1], 1.0, tol)

left   = Left()
right  = Right()
bottom = Bottom()
top    = Top()

# --------------------------------------------------------------
# 4. Dirichlet BCs
# --------------------------------------------------------------
noslip = Constant((0.0, 0.0))

# velocity BCs (no-slip everywhere)
bcu = [DirichletBC(W.sub(0), noslip, left),
       DirichletBC(W.sub(0), noslip, right),
       DirichletBC(W.sub(0), noslip, bottom),
       DirichletBC(W.sub(0), noslip, top)]

# temperature BCs (hot left, cold right)
T_left  = Constant(1.0)
T_right = Constant(0.0)
bct = [DirichletBC(W.sub(2), T_left, left),
        DirichletBC(W.sub(2), T_right, right)]

# pin pressure at a corner using the collapsed pressure space
class Corner(SubDomain):
    def inside(self, x, on_boundary): return near(x[0], 0.0, tol) and near(x[1], 0.0, tol)
p_bc = DirichletBC(W.sub(1), Constant(0.0), Corner(), method="pointwise")

# combine all BCs
bcs = bcu + bct + [p_bc]

# --------------------------------------------------------------
# 5. Physical parameters
# --------------------------------------------------------------
rho   = Constant(1.0)
mu    = Constant(1.5e-5)
alpha = Constant(2.1e-5)
gbeta = Constant(3.15e-5)
Tref  = Constant(0.5)

# --------------------------------------------------------------
# 6. Trial / Test functions
# --------------------------------------------------------------
w = Function(W)               # (u,p,T)
(u, p, T) = split(w)
(v, q, theta) = TestFunctions(W)

# --------------------------------------------------------------
# 7. Variational forms (steady Boussinesq)
# --------------------------------------------------------------
F_mom = rho*dot(dot(u, nabla_grad(u)), v)*dx \
        + 2*mu*inner(sym(grad(u)), sym(grad(v)))*dx \
        - p*div(v)*dx \
        - rho*gbeta*(T - Tref)*v[1]*dx

F_cont = q*div(u)*dx

F_en = dot(u, nabla_grad(T))*theta*dx + alpha*dot(grad(T), grad(theta))*dx

F = F_mom + F_cont + F_en

# --------------------------------------------------------------
# 8. Non-linear solve (Newton)
# --------------------------------------------------------------
solver = NewtonSolver()
solver.parameters["relative_tolerance"] = 1e-6
solver.parameters["absolute_tolerance"] = 1e-8
solver.parameters["maximum_iterations"] = 30
solver.parameters["linear_solver"] = "mumps"

solve(F == 0, w, bcs, J=derivative(F, w, TrialFunction(W)),
      solver_parameters={"newton_solver": {
          "relative_tolerance": 1e-6,
          "absolute_tolerance": 1e-8,
          "maximum_iterations": 30,
          "linear_solver": "mumps"}})

# --------------------------------------------------------------
# 9. Extract components
# --------------------------------------------------------------
(u_sol, p_sol, T_sol) = w.split(deepcopy=True)

# --------------------------------------------------------------
# 10. Average Nusselt number on left wall
# --------------------------------------------------------------
n = FacetNormal(mesh)
Nu = -dot(grad(T_sol), n)   # -∂T/∂n

# mark left boundary for ds measure
left_mark = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)
left.mark(left_mark, 1)
ds_left = Measure("ds", domain=mesh, subdomain_data=left_mark)

Nu_avg = assemble(Nu*ds_left(1))
print("Average Nusselt number on left wall: {:.6f}".format(Nu_avg))
with open("nusselt_number.json", "w", encoding="utf-8") as f:
    json.dump({
        "average_nusselt_number_left_wall": float(Nu_avg),
    }, f, indent=2)
    f.write("\n")

# --------------------------------------------------------------
# 11. Save results
# --------------------------------------------------------------
xdmf = XDMFFile(mesh.mpi_comm(), "q11_solution.xdmf")
xdmf.write(u_sol, 0.0)
xdmf.write(p_sol, 0.0)
xdmf.write(T_sol, 0.0)
xdmf.close()

plt.figure(figsize=(6,5))
p = plot(T_sol, title="Temperature")
plt.colorbar(p)
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig("q11_T.png", dpi=300)
plt.close()

# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
import json

mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

u_sol.rename("u", "velocity field")
p_sol.rename("p", "pressure field")
T_sol.rename("T", "temperature field")

solution_file = XDMFFile(mesh.mpi_comm(), "solution.xdmf")
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u_sol, 0.0)
solution_file.write(p_sol, 0.0)
solution_file.write(T_sol, 0.0)
solution_file.close()

checkpoint_file = XDMFFile(mesh.mpi_comm(), "read_checkpoint.xdmf")
checkpoint_file.write_checkpoint(u_sol, "u", 0.0, XDMFFile.Encoding.HDF5, False)
checkpoint_file.write_checkpoint(p_sol, "p", 0.0, XDMFFile.Encoding.HDF5, True)
checkpoint_file.write_checkpoint(T_sol, "T", 0.0, XDMFFile.Encoding.HDF5, True)
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
        {
            "symbol": "T_sol",
            "checkpoint_name": "T",
            "meaning": "temperature field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'T', -1)",
        },
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")
