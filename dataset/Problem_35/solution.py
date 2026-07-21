# Reference solution created from the benchmark statement for dataset Problem_35.
# Original benchmark: Multiphysics Problem 4 (Hard).

from __future__ import print_function

from dolfin import *
import json
import math
import numpy as np


# -----------------------------------------------------------------------------
# Geometry and mesh
# -----------------------------------------------------------------------------
Lx = math.pi
y_interface = 0.0

# Structured grid aligned with the Stokes-Darcy interface y = 0.
nx, ny = 80, 80
mesh = RectangleMesh(Point(0.0, -1.0), Point(Lx, 1.0), nx, ny, "crossed")


# -----------------------------------------------------------------------------
# Parameters from the benchmark statement
# -----------------------------------------------------------------------------
g_value = 1.0
rho_value = 1.0
nu_value = 1.0
k_value = 1.0
alpha_value = 1.0
mu_value = k_value * rho_value * g_value
K_value = k_value * rho_value * g_value / mu_value
darcy_lambda_value = k_value / mu_value

g = Constant(g_value)
rho = Constant(rho_value)
nu = Constant(nu_value)
darcy_lambda = Constant(darcy_lambda_value)
alpha_over_sqrt_k = Constant(alpha_value / math.sqrt(k_value))


# -----------------------------------------------------------------------------
# Cell and facet markers
# -----------------------------------------------------------------------------
STOKES, DARCY = 1, 2
cell_markers = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)


class StokesRegion(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] >= y_interface - DOLFIN_EPS


class DarcyRegion(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] <= y_interface + DOLFIN_EPS


StokesRegion().mark(cell_markers, STOKES)
DarcyRegion().mark(cell_markers, DARCY)
dxm = Measure("dx", domain=mesh, subdomain_data=cell_markers)

INTERFACE = 10
facet_markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)


class Interface(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], y_interface, 1.0e-12)


Interface().mark(facet_markers, INTERFACE)
dS = Measure("dS", domain=mesh, subdomain_data=facet_markers)


# -----------------------------------------------------------------------------
# Function space: global mixed space with pointwise pins outside each subdomain
#
# The physical weak form is integrated only on the correct subdomains. Extra
# pointwise zero constraints remove unused degrees of freedom from the other
# half of the mesh without constraining the interface trace.
# -----------------------------------------------------------------------------
P2v = VectorElement("Lagrange", mesh.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, MixedElement([P2v, P1, P1]))

w = Function(W)
(u_S, p_S, p_D) = TrialFunctions(W)
(v_S, q_S, psi_D) = TestFunctions(W)


# -----------------------------------------------------------------------------
# Exact expressions used for external Dirichlet data and diagnostics
# -----------------------------------------------------------------------------
w_expr = "-K - (g*x[1])/(2.0*nu) + (K/2.0 - alpha*g/(4.0*nu*nu))*x[1]*x[1]"
dw_expr = "-g/(2.0*nu) + (K - alpha*g/(2.0*nu*nu))*x[1]"

u_S_bc_expr = Expression(
    (
        "(%s)*cos(x[0])" % dw_expr,
        "(%s)*sin(x[0])" % w_expr,
    ),
    degree=4,
    K=K_value,
    g=g_value,
    nu=nu_value,
    alpha=alpha_value,
)

p_D_bc_expr = Expression(
    "rho*g*exp(x[1])*sin(x[0])",
    degree=4,
    rho=rho_value,
    g=g_value,
)

u_S_exact_restricted = Expression(
    (
        "x[1] >= -tol ? (%s)*cos(x[0]) : 0.0" % dw_expr,
        "x[1] >= -tol ? (%s)*sin(x[0]) : 0.0" % w_expr,
    ),
    degree=4,
    tol=1.0e-12,
    K=K_value,
    g=g_value,
    nu=nu_value,
    alpha=alpha_value,
)

p_D_exact_restricted = Expression(
    "x[1] <= tol ? rho*g*exp(x[1])*sin(x[0]) : 0.0",
    degree=4,
    tol=1.0e-12,
    rho=rho_value,
    g=g_value,
)

b_expr = Expression(
    (
        "((nu*K - (alpha*g)/(2.0*nu))*x[1] - g/2.0)*cos(x[0])",
        "(((nu*K)/2.0 - (alpha*g)/(4.0*nu))*x[1]*x[1] - (g/2.0)*x[1] + ((alpha*g)/(2.0*nu) - 2.0*nu*K))*sin(x[0])",
    ),
    degree=4,
    nu=nu_value,
    K=K_value,
    alpha=alpha_value,
    g=g_value,
)


# -----------------------------------------------------------------------------
# Boundary conditions
# -----------------------------------------------------------------------------
btol = 1.0e-12


class StokesExternalBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (
            near(x[1], 1.0, btol)
            or (near(x[0], 0.0, btol) and x[1] >= -btol)
            or (near(x[0], Lx, btol) and x[1] >= -btol)
        )


class DarcyExternalBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (
            near(x[1], -1.0, btol)
            or (near(x[0], 0.0, btol) and x[1] <= btol)
            or (near(x[0], Lx, btol) and x[1] <= btol)
        )


class LowerUnusedForStokes(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < -btol


class LowerUnusedForStokesPressure(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < -btol


class UpperUnusedForDarcyPressure(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] > btol


class StokesPressureGauge(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0, 1.0e-10) and near(x[1], 0.5, 1.0e-10)


zero_vector = Constant((0.0, 0.0))
zero_scalar = Constant(0.0)

bcs = [
    DirichletBC(W.sub(0), u_S_bc_expr, StokesExternalBoundary()),
    DirichletBC(W.sub(2), p_D_bc_expr, DarcyExternalBoundary()),
    DirichletBC(W.sub(0), zero_vector, LowerUnusedForStokes(), method="pointwise"),
    DirichletBC(W.sub(1), zero_scalar, LowerUnusedForStokesPressure(), method="pointwise"),
    DirichletBC(W.sub(2), zero_scalar, UpperUnusedForDarcyPressure(), method="pointwise"),
    DirichletBC(W.sub(1), zero_scalar, StokesPressureGauge(), method="pointwise"),
]


# -----------------------------------------------------------------------------
# Coupled weak form for the statement as written
#
# Stokes in Omega_S:
#   -div(sigma(u_S, p_S)) = b, div(u_S) = 0
#
# Darcy in Omega_D:
#   u_D = -(k/mu) grad(p_D), div(u_D) = 0
#
# Interface Gamma, n = (0, -1), t = (1, 0):
#   u_S.n = u_D.n
#   n.sigma.n = -p_D/rho
#   (alpha/sqrt(k)) u_S.t = -t.sigma.n
# -----------------------------------------------------------------------------
def eps(a):
    return sym(grad(a))


n_vec = Constant((0.0, -1.0))
t_vec = Constant((1.0, 0.0))

a = (
    2.0 * nu * inner(eps(u_S), eps(v_S)) * dxm(STOKES)
    - p_S * div(v_S) * dxm(STOKES)
    + q_S * div(u_S) * dxm(STOKES)
    + darcy_lambda * inner(grad(p_D), grad(psi_D)) * dxm(DARCY)
    + (1.0 / rho) * avg(p_D) * avg(dot(v_S, n_vec)) * dS(INTERFACE)
    + alpha_over_sqrt_k
    * avg(dot(u_S, t_vec))
    * avg(dot(v_S, t_vec)) * dS(INTERFACE)
    - avg(dot(u_S, n_vec)) * avg(psi_D) * dS(INTERFACE)
)

L_form = inner(b_expr, v_S) * dxm(STOKES)


try:
    linear_solver = "mumps" if has_lu_solver_method("mumps") else "lu"
except Exception:
    linear_solver = "lu"

print("Solving coupled Stokes-Darcy weak form...")
solve(a == L_form, w, bcs, solver_parameters={"linear_solver": linear_solver})

u_S_sol, p_S_sol, p_D_sol = w.split(deepcopy=True)
u_S_sol.rename("u_S", "Stokes velocity")
p_S_sol.rename("p_S", "Stokes pressure")
p_D_sol.rename("p_D", "Darcy pressure")

VD = VectorFunctionSpace(mesh, "DG", 0)
x_coord = SpatialCoordinate(mesh)
u_D_sol = project(
    conditional(
        le(x_coord[1], Constant(y_interface + 1.0e-12)),
        -darcy_lambda * grad(p_D_sol),
        as_vector((0.0, 0.0)),
    ),
    VD,
)
u_D_sol.rename("u_D", "Darcy velocity")


# -----------------------------------------------------------------------------
# Benchmark outputs: Stokes velocity and Darcy pressure in XDMF format
# -----------------------------------------------------------------------------
stokes_file = XDMFFile(mesh.mpi_comm(), "stokes_velocity.xdmf")
stokes_file.parameters["flush_output"] = True
stokes_file.write(u_S_sol, 0.0)
stokes_file.close()

darcy_file = XDMFFile(mesh.mpi_comm(), "darcy_pressure.xdmf")
darcy_file.parameters["flush_output"] = True
darcy_file.write(p_D_sol, 0.0)
darcy_file.close()


# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

solution_file = XDMFFile(mesh.mpi_comm(), "solution.xdmf")
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u_S_sol, 0.0)
solution_file.write(p_S_sol, 0.0)
solution_file.write(p_D_sol, 0.0)
solution_file.write(u_D_sol, 0.0)
solution_file.close()

checkpoint_file = XDMFFile(mesh.mpi_comm(), "read_checkpoint.xdmf")
checkpoint_file.write_checkpoint(u_S_sol, "u_S", 0.0, XDMFFile.Encoding.HDF5, False)
checkpoint_file.write_checkpoint(p_S_sol, "p_S", 0.0, XDMFFile.Encoding.HDF5, True)
checkpoint_file.write_checkpoint(p_D_sol, "p_D", 0.0, XDMFFile.Encoding.HDF5, True)
checkpoint_file.write_checkpoint(u_D_sol, "u_D", 0.0, XDMFFile.Encoding.HDF5, True)
checkpoint_file.close()

metadata = {
    "checkpoint_file": "read_checkpoint.xdmf",
    "solution_file": "solution.xdmf",
    "mesh_file": "mesh.xdmf",
    "time_index": -1,
    "fields": [
        {
            "checkpoint_name": "u_S",
            "source_symbol": "u_S_sol",
            "meaning": "Stokes velocity field on the upper subdomain, zero-pinned outside Omega_S",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'u_S', -1)",
        },
        {
            "checkpoint_name": "p_S",
            "source_symbol": "p_S_sol",
            "meaning": "Stokes pressure field on the upper subdomain, zero-pinned outside Omega_S",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'p_S', -1)",
        },
        {
            "checkpoint_name": "p_D",
            "source_symbol": "p_D_sol",
            "meaning": "Darcy pressure field on the lower subdomain, zero-pinned outside Omega_D",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'p_D', -1)",
        },
        {
            "checkpoint_name": "u_D",
            "source_symbol": "u_D_sol",
            "meaning": "Darcy velocity derived from u_D = -(k/mu) grad(p_D), zero outside Omega_D",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'u_D', -1)",
        },
    ],
}

with open("read_checkpoint.json", "w", encoding="utf-8") as handle:
    json.dump(metadata, handle, indent=2)
    handle.write("\n")


# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------
V_exact = VectorFunctionSpace(mesh, "Lagrange", 2)
Q_exact = FunctionSpace(mesh, "Lagrange", 1)
u_S_exact_fn = interpolate(u_S_exact_restricted, V_exact)
p_D_exact_fn = interpolate(p_D_exact_restricted, Q_exact)

u_D_exact = Expression(
    ("-rho*g*exp(x[1])*cos(x[0])", "-rho*g*exp(x[1])*sin(x[0])"),
    degree=4,
    rho=rho_value,
    g=g_value,
)

u_error = math.sqrt(assemble(inner(u_S_sol - u_S_bc_expr, u_S_sol - u_S_bc_expr) * dxm(STOKES)))
pS_error = math.sqrt(assemble(p_S_sol * p_S_sol * dxm(STOKES)))
pD_error = math.sqrt(assemble((p_D_sol - p_D_bc_expr) ** 2 * dxm(DARCY)))
uD_error = math.sqrt(assemble(inner(u_D_sol - u_D_exact, u_D_sol - u_D_exact) * dxm(DARCY)))
uS_linf = u_S_sol.vector().norm("linf")
pS_linf = p_S_sol.vector().norm("linf")
pD_linf = p_D_sol.vector().norm("linf")
uD_linf = u_D_sol.vector().norm("linf")

xs = np.linspace(0.0, Lx, 201)
w0 = -K_value
dw0 = -g_value / (2.0 * nu_value)
d2w0 = K_value - alpha_value * g_value / (2.0 * nu_value * nu_value)
exact_uS_n = -w0 * np.sin(xs)
exact_uD_n = (k_value / mu_value) * rho_value * g_value * np.sin(xs)
exact_pD = rho_value * g_value * np.sin(xs)
exact_normal_traction = 2.0 * nu_value * dw0 * np.sin(xs)
exact_sigma_xy = nu_value * (d2w0 + w0) * np.cos(xs)
exact_t_sigma_n = -exact_sigma_xy
exact_uS_t = dw0 * np.cos(xs)

manufactured_mass_residual = float(np.max(np.abs(exact_uS_n - exact_uD_n)))
manufactured_normal_residual = float(np.max(np.abs(exact_normal_traction + exact_pD / rho_value)))
manufactured_saffman_residual = float(
    np.max(np.abs((alpha_value / math.sqrt(k_value)) * exact_uS_t + exact_t_sigma_n))
)

print("Saved stokes_velocity.xdmf and darcy_pressure.xdmf")
print("Saved canonical dataset outputs")
print("Weak form solved with the intended Saffman interface condition.")
print("Physical-domain L2 error in u_S = %.12e" % u_error)
print("Physical-domain L2 error in p_S = %.12e" % pS_error)
print("Physical-domain L2 error in p_D = %.12e" % pD_error)
print("Physical-domain L2 error in derived u_D = %.12e" % uD_error)
print("Solution ||u_S||_linf = %.12e" % uS_linf)
print("Solution ||p_S||_linf = %.12e" % pS_linf)
print("Solution ||p_D||_linf = %.12e" % pD_linf)
print("Solution ||u_D||_linf = %.12e" % uD_linf)
print("Manufactured mass residual max = %.12e" % manufactured_mass_residual)
print("Manufactured normal-traction residual max = %.12e" % manufactured_normal_residual)
print("Manufactured Saffman residual max = %.12e" % manufactured_saffman_residual)
