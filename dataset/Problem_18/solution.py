# Copied from Results_ALL-FEM-main/reference solutions/fluid/2/fluid2.py for dataset Problem_18.
# Original benchmark: Fluid Mechanics Problem 2 (Easy).

# stokes_q2.py
from dolfin import *
import matplotlib.pyplot as plt

# ----- Geometry & mesh -----
L, H = 2.0, 0.20
nx, ny = 120, 12
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), nx, ny)

# ----- Materials / parameters -----
mu  = Constant(1.0)    # dynamic viscosity [Pa·s]
rho = Constant(1.0)    # density [kg·m^-3] (not used in steady Stokes)
f   = Constant((1.0, 0.0))  # body force [N·m^-3]

# ----- Function spaces (Taylor–Hood P2–P1) -----
Ve = VectorElement("Lagrange", mesh.ufl_cell(), 2)  # velocity
Qe = FiniteElement("Lagrange", mesh.ufl_cell(), 1)  # pressure
W  = FunctionSpace(mesh, MixedElement([Ve, Qe]))

# ----- Boundary conditions -----
# Walls: y = 0 and y = H  -> no slip
walls = CompiledSubDomain("near(x[1], 0.0) || near(x[1], H)", H=H)
bc_u_walls = DirichletBC(W.sub(0), Constant((0.0, 0.0)), walls)

# Pressure pin at a single point to remove constant nullspace
# (keeps traction-free inlet/outlet as natural conditions)
# Here we pin p = 0 at the bottom-left corner (0, 0).
bc_p_pin = DirichletBC(W.sub(1), Constant(0.0),
                       "near(x[0], 0.0) && near(x[1], 0.0)", method="pointwise")

bcs = [bc_u_walls, bc_p_pin]

# ----- Variational problem -----
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

# Symmetric gradient
def eps(w):
    return sym(grad(w))

a = (2.0*mu*inner(eps(u), eps(v)) - div(v)*p + q*div(u)) * dx
Lform = inner(f, v) * dx

# ----- Solve -----
w = Function(W, name="w")
solve(a == Lform, w, bcs,
      solver_parameters={"linear_solver": "mumps"})

(u_h, p_h) = w.split(deepcopy=True)
u_h.rename("u", "velocity")
p_h.rename("p", "pressure")

# ----- Save solution fields to XDMF -----
with XDMFFile(mesh.mpi_comm(), "q2_solution.xdmf") as xdmf:
    xdmf.parameters["flush_output"] = True
    xdmf.parameters["functions_share_mesh"] = True
    xdmf.write(u_h, 0.0)
    xdmf.write(p_h, 0.0)

# ----- Save a colour map of speed |u| -----
speed = project(sqrt(dot(u_h, u_h)), FunctionSpace(mesh, "CG", 1))
speed.set_allow_extrapolation(True)

plt.figure()
c = plot(speed)
plt.colorbar(c)
plt.title("|u| (speed)")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.tight_layout()
plt.savefig("q2_speed.png", dpi=200)
plt.close()

print("Done. Wrote q2_solution.xdmf (u, p) and q2_speed.png.")

# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
import json

mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

u_h.rename("u", "velocity field")
p_h.rename("p", "pressure field")

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
            "meaning": "velocity field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'u', -1)",
        },
        {
            "symbol": "p_h",
            "checkpoint_name": "p",
            "meaning": "pressure field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'p', -1)",
        },
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")
