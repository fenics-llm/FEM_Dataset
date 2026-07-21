# Copied from Results_ALL-FEM-main/reference solutions/fluid/6/fluid6.py for dataset Problem_22.
# Original benchmark: Fluid Mechanics Problem 6 (Medium).

# Backward-facing step – steady incompressible Navier–Stokes (legacy dolfin)

from dolfin import *
from mshr import Rectangle, generate_mesh
import matplotlib.pyplot as plt
import json

# ------ Geometry & mesh (single conforming mesh) ------
H      = 1.0                     # step height (m)
L_up   = 3.0 * H                 # upstream length
L_down = 20.0 * H                # downstream length

# The fluid domain is the L-shaped union of the upstream channel and the
# expanded downstream channel. A full enclosing rectangle would incorrectly
# include the upstream upper solid block and turn the step wall into an
# interior line with no no-slip boundary condition.
upstream_channel = Rectangle(Point(-L_up, 0.0), Point(0.0, H))
downstream_channel = Rectangle(Point(0.0, 0.0), Point(L_down, 2.0*H))
mesh_resolution = 128
mesh = generate_mesh(upstream_channel + downstream_channel, mesh_resolution)

# ------ Boundary markers (single mesh) ------
tol = 1.0e-8
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)

class Inlet(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], -L_up, tol) and (0.0 - tol <= x[1] <= H + tol)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.0, tol)

class TopUp(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], H, tol) and (x[0] <= 0.0 + tol)

class TopDown(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 2.0*H, tol) and (x[0] >= 0.0 - tol)

class StepWall(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0, tol) and (x[1] >= H - tol)

class Outlet(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], L_down, tol) and (0.0 - tol <= x[1] <= 2.0*H + tol)

Inlet().mark(boundaries, 1)
Bottom().mark(boundaries, 2)
TopUp().mark(boundaries, 3)
TopDown().mark(boundaries, 4)
StepWall().mark(boundaries, 5)
Outlet().mark(boundaries, 6)

ds = Measure("ds", domain=mesh, subdomain_data=boundaries)

boundary_lengths = {
    "inlet": assemble(Constant(1.0)*ds(1)),
    "bottom": assemble(Constant(1.0)*ds(2)),
    "upstream_top": assemble(Constant(1.0)*ds(3)),
    "downstream_top": assemble(Constant(1.0)*ds(4)),
    "step_wall": assemble(Constant(1.0)*ds(5)),
    "outlet": assemble(Constant(1.0)*ds(6)),
}
missing_boundaries = [
    name for name, length in boundary_lengths.items() if length <= 1.0e-10
]
if missing_boundaries:
    raise RuntimeError(
        "Missing boundary markers on the backward-facing-step mesh: "
        + ", ".join(missing_boundaries)
    )

# ------ Function spaces (Taylor–Hood) ------
Ve = VectorElement("Lagrange", mesh.ufl_cell(), 2)
Pe = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = MixedElement([Ve, Pe])
W  = FunctionSpace(mesh, TH)

# ------ Boundary conditions ------
U_bar = 1.0
mu    = 1.0e-2
rho   = 1.0

inlet_expr = Expression(("6.0*U_bar*x[1]/H*(1.0 - x[1]/H)", "0.0"),
                        degree=2, U_bar=U_bar, H=H)

noslip = Constant((0.0, 0.0))
bcu_inlet   = DirichletBC(W.sub(0), inlet_expr, boundaries, 1)
bcu_bottom  = DirichletBC(W.sub(0), noslip, boundaries, 2)
bcu_topup   = DirichletBC(W.sub(0), noslip, boundaries, 3)
bcu_topdown = DirichletBC(W.sub(0), noslip, boundaries, 4)
bcu_step    = DirichletBC(W.sub(0), noslip, boundaries, 5)

class Corner(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], -L_up, tol) and near(x[1], 0.0, tol)

bcp_corner = DirichletBC(W.sub(1), Constant(0.0), Corner(), method="pointwise")
bcs = [bcu_inlet, bcu_bottom, bcu_topup, bcu_topdown, bcu_step, bcp_corner]

# ------ Variational formulation (steady Navier–Stokes) ------
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

U = Function(W)          # current iterate
(u_k, p_k) = split(U)    # for Newton linearisation

def epsilon(w):
    return sym(grad(w))

# ---- corrected Navier–Stokes residual ----
F = ( rho*dot(dot(u_k, nabla_grad(u_k)), v)*dx
      + 2*mu*inner(epsilon(u_k), epsilon(v))*dx
      - div(v)*p_k*dx
      + q*div(u_k)*dx )
J = derivative(F, U, TrialFunction(W))

# ------ Stokes pre-solve (good initial guess) ------
U0 = Function(W)                     # zero initial guess
(u0, p0) = split(U0)

# ---- corrected Stokes residual ----
F_stokes = ( 2*mu*inner(epsilon(u0), epsilon(v))*dx
             - div(v)*p0*dx
             + q*div(u0)*dx )
J_stokes = derivative(F_stokes, U0, TrialFunction(W))

solve(F_stokes == 0, U0, bcs, J=J_stokes,
      solver_parameters={"newton_solver": {"linear_solver": "mumps"}})   # legacy key

U.assign(U0)   # use Stokes solution as Newton start

# ------ Newton solve ------
solve(F == 0, U, bcs, J=J,
      solver_parameters={"newton_solver":
                         {"linear_solver": "mumps",
                          "relative_tolerance": 1e-6,
                          "absolute_tolerance": 1e-8,
                          "maximum_iterations": 80,
                          "relaxation_parameter": 0.7,
                          "error_on_nonconvergence": False}})

# ------ Split solution and post-processing ------
(u_sol, p_sol) = U.split(deepcopy=True)

# Wall shear stress on downstream top wall (y = 2H)
Vdg = FunctionSpace(mesh, "DG", 0)   # facetwise scalar
tau_w = Function(Vdg)

grad_u = grad(u_sol)
tau_expr = mu * grad_u[0, 1]          # ∂u_x/∂y
tau_proj = project(tau_expr, Vdg, solver_type="cg")
tau_w.assign(tau_proj)

# Keep cells adjacent to the downstream top wall (boundary id 4)
tau_w_top = Function(Vdg)
tau_w_top.vector()[:] = 0.0
dg0_dofmap = Vdg.dofmap()
top_wall_samples = []
for facet in facets(mesh):
    if boundaries[facet.index()] == 4:
        x_mid = float(facet.midpoint().x())
        for cell in cells(facet):
            dof = dg0_dofmap.cell_dofs(cell.index())[0]
            tau_w_top.vector()[dof] = tau_w.vector()[dof]
            top_wall_samples.append((x_mid, float(tau_w.vector()[dof])))

# Reattachment point: first downstream zero crossing of wall shear on y = 2H.
tau_samples = sorted(top_wall_samples, key=lambda item: item[0])

reattachment_x = None
for (x0, tau0), (x1, tau1) in zip(tau_samples[:-1], tau_samples[1:]):
    if abs(tau0) <= 1.0e-14:
        reattachment_x = x0
        break
    if tau0*tau1 < 0.0:
        reattachment_x = x0 - tau0*(x1 - x0)/(tau1 - tau0)
        break
if reattachment_x is None and tau_samples and abs(tau_samples[-1][1]) <= 1.0e-14:
    reattachment_x = tau_samples[-1][0]

with open("reattachment_point.json", "w", encoding="utf-8") as f:
    json.dump({
        "quantity": "downstream_top_wall_shear_zero_crossing",
        "found": reattachment_x is not None,
        "x_m": None if reattachment_x is None else float(reattachment_x),
        "wall": "y = 2H downstream of the step",
        "method": "linear interpolation between adjacent downstream-top-wall facet samples",
        "sample_count": len(tau_samples),
    }, f, indent=2)
    f.write("\n")

# ------ Output ------
# Velocity magnitude PNG
u_mag = sqrt(dot(u_sol, u_sol))
u_mag = project(u_mag, FunctionSpace(mesh, "CG", 2))
plt.figure()
p = plot(u_mag, title="Velocity magnitude")
plt.colorbar(p)
plt.savefig("q6_u.png", dpi=300)

# XDMF solution
with XDMFFile(mesh.mpi_comm(), "q6_soln.xdmf") as xdmf:
    xdmf.write(u_sol, 0.0)
    xdmf.write(p_sol, 0.0)

# Shear stress on downstream top wall
with XDMFFile(mesh.mpi_comm(), "q6_tau_top.xdmf") as xdmf:
    xdmf.write(tau_w_top, 0.0)

# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
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
