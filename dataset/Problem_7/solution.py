# Copied from Results_ALL-FEM-main/reference solutions/solid/7/solid7.py for dataset Problem_7.
# Original benchmark: Solid Mechanics Problem 7 (Easy).

from dolfin import *

# --- Geometry & mesh (structured 80 x 16 over (0,1.0) x (0,0.20) m) ---
Lx, Ly = 1.0, 0.20
mesh = RectangleMesh(Point(0.0, 0.0), Point(Lx, Ly), 80, 16, "left/right")

# --- Subdomains via y-split (top: Aluminium, bottom: Steel) ---
y_split = 0.10  # interface at y = 0.10 m

# --- Material properties (piecewise via Expressions) ---
E_Al, nu_Al = 70e9, 0.30
E_St, nu_St = 200e9, 0.30

# Spatially varying Young's modulus and Poisson's ratio using the y-split
E  = Expression("x[1] >= y0 ? E_al : E_st", degree=0, y0=y_split, E_al=E_Al, E_st=E_St)
nu = Expression("x[1] >= y0 ? nu_al : nu_st", degree=0, y0=y_split, nu_al=nu_Al, nu_st=nu_St)

# --- Function space ---
V = VectorFunctionSpace(mesh, "CG", 1)

# --- Kinematics & plane-stress constitutive law ---
def eps(u):
    return sym(grad(u))

# Plane stress: mu = E/(2(1+nu)), lambda_ps = 2*mu*nu/(1 - nu) = E*nu/(1 - nu*nu)
def sigma(u):
    mu  = E/(2.0*(1.0 + nu))
    lam = 2.0*mu*nu/(1.0 - nu)  # plane stress effective lambda
    return 2.0*mu*eps(u) + lam*tr(eps(u))*Identity(2)

# --- Boundary markers for Neumann data ---
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0)

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], Lx)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.0)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], Ly)

left = Left(); right = Right(); bottom = Bottom(); top = Top()
left.mark(boundaries, 1)
right.mark(boundaries, 2)
bottom.mark(boundaries, 3)
top.mark(boundaries, 4)

ds = Measure("ds", domain=mesh, subdomain_data=boundaries)

# --- Dirichlet BCs (left edge fixed) ---
u0 = Constant((0.0, 0.0))
bc_left = DirichletBC(V, u0, boundaries, 1)

# --- Neumann traction on right edge: (0, -5000) N/m ---
t_right = Constant((0.0, -5000.0))

# --- Variational problem ---
u = TrialFunction(V)
v = TestFunction(V)

a = inner(sigma(u), eps(v)) * dx
L = dot(t_right, v) * ds(2)  # only on the right edge

# --- Solve ---
u_h = Function(V, name="displacement")
solve(a == L, u_h, [bc_left], solver_parameters={"linear_solver": "mumps"})

# --- Output 1: XDMF of displacement field ---
with XDMFFile(mesh.mpi_comm(), "q7_disp.xdmf") as xdmf:
    xdmf.parameters["flush_output"] = True
    xdmf.parameters["functions_share_mesh"] = True
    xdmf.write(u_h)

# --- Output 2: colour map PNG of |u| ---
# Compute displacement magnitude
Vsc = FunctionSpace(mesh, "CG", 1)
umag = project(sqrt(dot(u_h, u_h)), Vsc, solver_type="cg", preconditioner_type="amg")
umag.rename("u_magnitude", "")

# Save a high-resolution PNG using matplotlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.figure()
c = plot(umag, title="Displacement magnitude |u|")
plt.colorbar(c)
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.tight_layout()
plt.savefig("q7_disp.png", dpi=300)

# Optional: also write the scalar field to XDMF (handy for post-processing)
with XDMFFile(mesh.mpi_comm(), "q7_umag.xdmf") as xdmf_umag:
    xdmf_umag.parameters["flush_output"] = True
    xdmf_umag.parameters["functions_share_mesh"] = True
    xdmf_umag.write(umag)

print("Solved. Outputs: q7_disp.xdmf (displacement), q7_disp.png (colour map), q7_umag.xdmf (|u|).")

# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
import json

mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

u_h.rename("u", "displacement field")

solution_file = XDMFFile(mesh.mpi_comm(), "solution.xdmf")
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u_h, 0.0)
solution_file.close()

checkpoint_file = XDMFFile(mesh.mpi_comm(), "read_checkpoint.xdmf")
checkpoint_file.write_checkpoint(u_h, "u", 0.0, XDMFFile.Encoding.HDF5, False)
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
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")
