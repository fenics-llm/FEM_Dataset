# Copied from Results_ALL-FEM-main/reference solutions/multiphysics - Copy/7/multi-7.py for dataset Problem_38.
# Original benchmark: Multiphysics Problem 7 (Hard).

# file: q14_corrected.py
from __future__ import print_function

from dolfin import *
import numpy as np


# -----------------------------------------------------------------------------
# Geometry and mesh
# -----------------------------------------------------------------------------
L = 1.0
H = 0.20
xa, xb = 0.4, 0.6                 # porous block x-extent

nx, ny = 200, 40                  # aligned with x = 0.4 and x = 0.6
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), nx, ny, "crossed")


# -----------------------------------------------------------------------------
# Physical parameters, SI units
# -----------------------------------------------------------------------------
rho = 1.0                         # kg/m^3
mu = 0.01                         # Pa*s
K = 1.0e-6                        # m^2
Ubar = 1.0                        # m/s

rho_c = Constant(rho)
mu_c = Constant(mu)
drag_c = Constant(mu / K)


# -----------------------------------------------------------------------------
# Cell markers: 0 = free fluid Omega_f, 1 = porous block Pi
# -----------------------------------------------------------------------------
cell_markers = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)


class PorousBlock(SubDomain):
    def inside(self, x, on_boundary):
        # For cell marking this is evaluated at cell points/midpoints.
        return (x[0] >= xa - DOLFIN_EPS) and (x[0] <= xb + DOLFIN_EPS)


PorousBlock().mark(cell_markers, 1)
dxm = Measure("dx", domain=mesh, subdomain_data=cell_markers)


# -----------------------------------------------------------------------------
# Boundary markers
# -----------------------------------------------------------------------------
INLET, OUTLET, BOTTOM, TOP = 1, 2, 3, 4
facet_markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)
btol = 1.0e-12


class Inlet(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0, btol)


class Outlet(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], L, btol)


class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.0, btol)


class Top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], H, btol)


Inlet().mark(facet_markers, INLET)
Outlet().mark(facet_markers, OUTLET)
Bottom().mark(facet_markers, BOTTOM)
Top().mark(facet_markers, TOP)

# ds is defined for completeness; the traction-free outlet is imposed naturally.
ds = Measure("ds", domain=mesh, subdomain_data=facet_markers)


# -----------------------------------------------------------------------------
# Function space: Taylor-Hood P2/P1 mixed element
# -----------------------------------------------------------------------------
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, MixedElement([P2, P1]))

w = Function(W)
(u, p) = split(w)
(v, q) = TestFunctions(W)


# -----------------------------------------------------------------------------
# Boundary conditions
# -----------------------------------------------------------------------------
inlet_velocity = Expression(
    ("6.0*Ubar*x[1]*(H - x[1])/(H*H)", "0.0"),
    degree=2,
    Ubar=Ubar,
    H=H,
)

zero_velocity = Constant((0.0, 0.0))

bc_inlet = DirichletBC(W.sub(0), inlet_velocity, facet_markers, INLET)
bc_bottom = DirichletBC(W.sub(0), zero_velocity, facet_markers, BOTTOM)
bc_top = DirichletBC(W.sub(0), zero_velocity, facet_markers, TOP)

# No pressure Dirichlet condition is applied: the outlet traction-free condition
# is the natural boundary condition of the stress weak form.
bcs = [bc_inlet, bc_bottom, bc_top]


# -----------------------------------------------------------------------------
# Operators
# -----------------------------------------------------------------------------
def eps(a):
    return sym(grad(a))


# -----------------------------------------------------------------------------
# Variational problem
#
# In Omega_f:
#   rho (u . grad)u = -grad(p) + div(2 mu eps(u))
#
# In Pi:
#   0 = -grad(p) + div(2 mu eps(u)) - (mu/K)u
#
# Because div(u)=0 and mu is constant, div(2 mu eps(u)) = mu Laplacian(u).
# The conforming mixed space plus one stress weak form gives velocity continuity
# and weak traction continuity across the internal interface.
# -----------------------------------------------------------------------------
F = (
    rho_c * inner(grad(u) * u, v) * dxm(0)       # convection only in free fluid
    + 2.0 * mu_c * inner(eps(u), eps(v)) * dxm   # viscous stress everywhere
    - p * div(v) * dxm                           # pressure term everywhere
    + q * div(u) * dxm                           # incompressibility everywhere
    + drag_c * inner(u, v) * dxm(1)              # Brinkman drag only in Pi
)

J = derivative(F, w, TrialFunction(W))


# -----------------------------------------------------------------------------
# Optional Stokes-Brinkman initial guess, improves Newton robustness
# -----------------------------------------------------------------------------
try:
    linear_solver = "mumps" if has_lu_solver_method("mumps") else "lu"
except Exception:
    linear_solver = "lu"

w0 = TrialFunction(W)
(u0, p0) = split(w0)

A0 = (
    2.0 * mu_c * inner(eps(u0), eps(v)) * dxm
    - p0 * div(v) * dxm
    + q * div(u0) * dxm
    + drag_c * inner(u0, v) * dxm(1)
)
L0 = inner(Constant((0.0, 0.0)), v) * dxm

print("Solving Stokes-Brinkman initial guess...")
solve(A0 == L0, w, bcs, solver_parameters={"linear_solver": linear_solver})


# -----------------------------------------------------------------------------
# Nonlinear Navier-Stokes/Brinkman solve
# -----------------------------------------------------------------------------
problem = NonlinearVariationalProblem(F, w, bcs, J)
solver = NonlinearVariationalSolver(problem)

prm = solver.parameters
prm["nonlinear_solver"] = "newton"
prm["newton_solver"]["linear_solver"] = linear_solver
prm["newton_solver"]["absolute_tolerance"] = 1.0e-10
prm["newton_solver"]["relative_tolerance"] = 1.0e-8
prm["newton_solver"]["maximum_iterations"] = 40
prm["newton_solver"]["relaxation_parameter"] = 1.0
try:
    prm["newton_solver"]["error_on_nonconvergence"] = True
except Exception:
    pass

print("Solving nonlinear Navier-Stokes/Brinkman system...")
solver.solve()

u_sol, p_sol = w.split(deepcopy=True)
u_sol.rename("u", "velocity")
p_sol.rename("p", "pressure")


# -----------------------------------------------------------------------------
# Output 1: velocity and pressure field to XDMF
# -----------------------------------------------------------------------------
xdmf = XDMFFile(mesh.mpi_comm(), "q14_solution.xdmf")
xdmf.parameters["flush_output"] = True
xdmf.parameters["functions_share_mesh"] = True
xdmf.write(u_sol, 0.0)
xdmf.write(p_sol, 0.0)
xdmf.close()
print("Saved q14_solution.xdmf")


# -----------------------------------------------------------------------------
# Output 2: pressure drop across porous block at centerline
# -----------------------------------------------------------------------------
y_c = 0.5 * H
h_x = L / float(nx)
x_before = max(xa - 0.5 * h_x, DOLFIN_EPS)
x_after = min(xb + 0.5 * h_x, L - DOLFIN_EPS)

p_before = float(p_sol(Point(x_before, y_c)))
p_after = float(p_sol(Point(x_after, y_c)))
delta_p = p_before - p_after

# Save only the scalar value, which is safest for automated graders.
with open("q14_dp.txt", "w") as f:
    f.write("%.12e\n" % delta_p)

print("p_before = %.8e Pa at x = %.6f, y = %.6f" % (p_before, x_before, y_c))
print("p_after  = %.8e Pa at x = %.6f, y = %.6f" % (p_after, x_after, y_c))
print("Delta p  = %.8e Pa" % delta_p)
print("Saved q14_dp.txt")


# -----------------------------------------------------------------------------
# Output 3: color map of velocity magnitude |u|
# -----------------------------------------------------------------------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation
    from matplotlib.patches import Rectangle

    Vmag = FunctionSpace(mesh, "CG", 1)
    speed = project(sqrt(dot(u_sol, u_sol)), Vmag)

    coordinates = mesh.coordinates()
    triangles = mesh.cells()
    vertex_values = speed.compute_vertex_values(mesh)

    triangulation = Triangulation(coordinates[:, 0], coordinates[:, 1], triangles)

    fig, ax = plt.subplots(figsize=(10.0, 2.4))
    image = ax.tripcolor(triangulation, vertex_values, shading="gouraud")

    ax.add_patch(Rectangle((xa, 0.0), xb - xa, H, fill=False, linewidth=1.2))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(0.0, L)
    ax.set_ylim(0.0, H)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("Velocity magnitude |u|; porous block outlined")

    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("|u| (m/s)")

    fig.tight_layout()
    fig.savefig("q14_speed.png", dpi=220)
    plt.close(fig)
    print("Saved q14_speed.png")
except Exception as exc:
    print("Failed to save q14_speed.png: %s" % exc)

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
