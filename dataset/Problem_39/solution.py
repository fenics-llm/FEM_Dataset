# Reference solution created from the benchmark statement for dataset Problem_39.
# Original benchmark: Multiphysics Problem 8 (Hard).

from __future__ import print_function

from dolfin import *
import csv
import json
import numpy as np


# -----------------------------------------------------------------------------
# Geometry and mesh
# -----------------------------------------------------------------------------
L = 1.0
H = 0.20
x_interface = 0.60

# Structured grid aligned with x = 0.60.
nx, ny = 160, 32
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), nx, ny, "crossed")


# -----------------------------------------------------------------------------
# Physical parameters, SI units
# -----------------------------------------------------------------------------
mu_value = 0.02
K_value = 1.0e-6
Ubar_value = 0.1

mu = Constant(mu_value)
darcy_drag = Constant(mu_value / K_value)


# -----------------------------------------------------------------------------
# Cell markers: 0 = free fluid Omega_f, 1 = porous medium Omega_p
# -----------------------------------------------------------------------------
FLUID, POROUS = 0, 1
cell_markers = MeshFunction("size_t", mesh, mesh.topology().dim(), FLUID)


class PorousMedium(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] >= x_interface - DOLFIN_EPS


PorousMedium().mark(cell_markers, POROUS)
dxm = Measure("dx", domain=mesh, subdomain_data=cell_markers)


# -----------------------------------------------------------------------------
# Function space: Taylor-Hood P2/P1 mixed element for combined (u, p)
# -----------------------------------------------------------------------------
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, MixedElement([P2, P1]))

w = Function(W)
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)


# -----------------------------------------------------------------------------
# Boundary and interface conditions
# -----------------------------------------------------------------------------
tol = 1.0e-12


class Inlet(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0, tol)


class Outlet(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], L, tol)


class FluidBottom(SubDomain):
    def inside(self, x, on_boundary):
        return (
            on_boundary
            and near(x[1], 0.0, tol)
            and x[0] <= x_interface + tol
        )


class FluidTop(SubDomain):
    def inside(self, x, on_boundary):
        return (
            on_boundary
            and near(x[1], H, tol)
            and x[0] <= x_interface + tol
        )


class PorousBottom(SubDomain):
    def inside(self, x, on_boundary):
        return (
            on_boundary
            and near(x[1], 0.0, tol)
            and x[0] >= x_interface - tol
        )


class PorousTop(SubDomain):
    def inside(self, x, on_boundary):
        return (
            on_boundary
            and near(x[1], H, tol)
            and x[0] >= x_interface - tol
        )


class Interface(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], x_interface, tol)


inlet_velocity = Expression(
    ("6.0*Ubar*x[1]*(H - x[1])/(H*H)", "0.0"),
    degree=2,
    Ubar=Ubar_value,
    H=H,
)

zero_velocity = Constant((0.0, 0.0))
zero_scalar = Constant(0.0)

bcs = [
    DirichletBC(W.sub(0), inlet_velocity, Inlet()),
    DirichletBC(W.sub(0), zero_velocity, FluidBottom()),
    DirichletBC(W.sub(0), zero_velocity, FluidTop()),
    DirichletBC(W.sub(0).sub(1), zero_scalar, PorousBottom()),
    DirichletBC(W.sub(0).sub(1), zero_scalar, PorousTop()),
    # Tangential direction on the vertical interface is y, so u_y = 0.
    DirichletBC(W.sub(0).sub(1), zero_scalar, Interface(), method="pointwise"),
    DirichletBC(W.sub(1), zero_scalar, Outlet()),
]


# -----------------------------------------------------------------------------
# Weak form
#
# The benchmark asks for Stokes flow in Omega_f and Darcy flow in Omega_p.
# In this legacy single-mesh FEniCS implementation, u and p are conforming
# fields over both subdomains. This strongly enforces pressure continuity and
# normal-velocity continuity across x = 0.60. The porous subdomain uses Darcy's
# law only: (mu/K) u + grad(p) = 0, written in mixed weak form after integrating
# the pressure-gradient term by parts.
# -----------------------------------------------------------------------------
def eps(a):
    return sym(grad(a))


a = (
    2.0 * mu * inner(eps(u), eps(v)) * dxm(FLUID)
    + darcy_drag * inner(u, v) * dxm(POROUS)
    - p * div(v) * dxm
    + q * div(u) * dxm
)
rhs = Constant(0.0) * q * dxm


try:
    linear_solver = "mumps" if has_lu_solver_method("mumps") else "lu"
except Exception:
    linear_solver = "lu"

print("Solving coupled Stokes-Darcy system...")
solve(a == rhs, w, bcs, solver_parameters={"linear_solver": linear_solver})

u_sol, p_sol = w.split(deepcopy=True)
u_sol.rename("u", "combined velocity")
p_sol.rename("p", "pressure")


# -----------------------------------------------------------------------------
# Benchmark outputs
# -----------------------------------------------------------------------------
q15 = XDMFFile(mesh.mpi_comm(), "q15_solution.xdmf")
q15.parameters["flush_output"] = True
q15.parameters["functions_share_mesh"] = True
q15.write(u_sol, 0.0)
q15.write(p_sol, 0.0)
q15.close()
print("Saved q15_solution.xdmf")

sample_count = 201
ys = np.linspace(0.0, H, sample_count)
ux_values = []
uy_values = []

with open("q15_interface.csv", "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["y", "u_x", "u_y"])
    for y in ys:
        value = u_sol(Point(x_interface, float(y)))
        ux_values.append(float(value[0]))
        uy_values.append(float(value[1]))
        writer.writerow([float(y), float(value[0]), float(value[1])])

print("Saved q15_interface.csv")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation

    coordinates = mesh.coordinates()
    cells = mesh.cells()
    pressure_values = p_sol.compute_vertex_values(mesh)
    triangulation = Triangulation(coordinates[:, 0], coordinates[:, 1], cells)

    fig, ax = plt.subplots(figsize=(9.0, 2.4))
    image = ax.tripcolor(triangulation, pressure_values, shading="gouraud")
    ax.axvline(x_interface, color="white", linewidth=0.8)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(0.0, L)
    ax.set_ylim(0.0, H)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("Pressure field p")
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("p (Pa)")
    fig.tight_layout()
    fig.savefig("q15_p.png", dpi=220)
    plt.close(fig)
    print("Saved q15_p.png")
except Exception as exc:
    print("Failed to save q15_p.png: %s" % exc)


# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

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

metadata = {
    "checkpoint_file": "read_checkpoint.xdmf",
    "solution_file": "solution.xdmf",
    "mesh_file": "mesh.xdmf",
    "time_index": -1,
    "fields": [
        {
            "checkpoint_name": "u",
            "source_symbol": "u_sol",
            "meaning": "combined velocity field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'u', -1)",
        },
        {
            "checkpoint_name": "p",
            "source_symbol": "p_sol",
            "meaning": "pressure field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'p', -1)",
        },
    ],
}

with open("read_checkpoint.json", "w", encoding="utf-8") as handle:
    json.dump(metadata, handle, indent=2)
    handle.write("\n")

print("Saved canonical dataset outputs")


# -----------------------------------------------------------------------------
# Console summary
# -----------------------------------------------------------------------------
ux_values = np.array(ux_values)
uy_values = np.array(uy_values)
interface_flux = float(np.trapz(ux_values, ys))
mean_interface_ux = interface_flux / H
pressure_dofs = p_sol.vector().get_local()

print("Interface flux integral int u_x dy = %.12e m^2/s" % interface_flux)
print("Mean interface u_x = %.12e m/s" % mean_interface_ux)
print("Interface u_x min/max = %.12e / %.12e m/s" % (ux_values.min(), ux_values.max()))
print("Interface |u_y| max = %.12e m/s" % np.abs(uy_values).max())
print("Pressure dof min/max = %.12e / %.12e Pa" % (pressure_dofs.min(), pressure_dofs.max()))
