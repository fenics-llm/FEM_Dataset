# Copied from Results_ALL-FEM-main/reference solutions/solid/15/solid15.py for dataset Problem_15.
# Original benchmark: Solid Mechanics Problem 15 (Hard).

# q15_svk_strip_3holes.py
#
# Legacy FEniCS (dolfin + mshr) solution:
# Finite-strain Saint-Venant–Kirchhoff (SVK), plane strain,
# perforated strip with 3 holes, displacement-driven loading,
# load stepping with Newton, stop when max principal GL strain <= 0.03.

from dolfin import *
from mshr import *
import numpy as np

# Non-interactive plotting
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------
# 0) Problem parameters
# ----------------------------
Lx, Ly = 1.20, 0.20
a = 0.03
hole_centers = [(0.30, 0.10), (0.60, 0.10), (0.90, 0.10)]

# Material (convert MPa -> Pa for SI consistency)
lmbda = Constant(5.769e6)
mu    = Constant(3.846e6)

# Prescribed displacement on right edge
ux_total = 0.012  # m

# Strain cap
Emax_cap = 0.03

# Mesh resolution (increase if you need more accuracy around holes)
resolution = 50

# Newton solver controls
newton_abs_tol = 1e-8
newton_rel_tol = 1e-7
newton_max_it  = 30

# Load stepping controls
t_end   = 1.0
dt      = 0.10     # initial step
dt_min  = 1e-4     # smallest allowed step
dt_max  = 0.20     # maximum step (optional growth)


# ----------------------------
# 1) Mesh: rectangle minus 3 circles
# ----------------------------
rect = Rectangle(Point(0.0, 0.0), Point(Lx, Ly))

domain = rect
for (cx, cy) in hole_centers:
    domain = domain - Circle(Point(cx, cy), a, 64)

mesh = generate_mesh(domain, resolution)

# ----------------------------
# 2) Function space
# ----------------------------
V = VectorFunctionSpace(mesh, "CG", 2)  # quadratic displacement is common for hyperelasticity
u = Function(V, name="u")               # unknown
du = TrialFunction(V)                   # increment
v  = TestFunction(V)                    # test


# ----------------------------
# 3) Boundary definitions
# ----------------------------
class LeftBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0, 1e-8)

class RightBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], Lx, 1e-8)

left = LeftBoundary()
right = RightBoundary()

# Time/load-dependent right displacement (only u_x changes)
u_R = Expression(("t*ux_total", "0.0"), t=0.0, ux_total=ux_total, degree=1)

bc_left      = DirichletBC(V, Constant((0.0, 0.0)), left)
bc_right     = DirichletBC(V, u_R, right)

bcs = [bc_left, bc_right]


# ----------------------------
# 4) Kinematics + SVK energy
# ----------------------------
I2 = Identity(2)
F2 = I2 + grad(u)            # deformation gradient (2D)

# Embed into 3D for True Plane Strain (matches Reference)
F3 = as_tensor([[F2[0, 0], F2[0, 1], 0.0],
                [F2[1, 0], F2[1, 1], 0.0],
                [0.0,      0.0,      1.0]])

I3 = Identity(3)
E3 = 0.5*(F3.T*F3 - I3)      # 3x3 Green–Lagrange strain
S3 = lmbda*tr(E3)*I3 + 2.0*mu*E3     # 3x3 2nd PK stress
P3 = F3*S3                           # 3x3 1st PK stress

# In-plane part used in the weak form
P2 = as_tensor([[P3[0, 0], P3[0, 1]],
                [P3[1, 0], P3[1, 1]]])

# Residual and consistent tangent (automatic differentiation)
R = inner(P2, grad(v))*dx
J = derivative(R, u, du)


# ----------------------------
# 5) Nonlinear solver setup
# ----------------------------
problem = NonlinearVariationalProblem(R, u, bcs, J)
solver  = NonlinearVariationalSolver(problem)

prm = solver.parameters["newton_solver"]
prm["absolute_tolerance"] = newton_abs_tol
prm["relative_tolerance"] = newton_rel_tol
prm["maximum_iterations"] = newton_max_it
prm["relaxation_parameter"] = 1.0
prm["error_on_nonconvergence"] = True
# Use a robust direct solver if available
prm["linear_solver"] = "mumps"

set_log_level(LogLevel.PROGRESS)


# ----------------------------
# 6) Post-processing expressions
# ----------------------------
# Max principal Green-Lagrange strain in 2D (matches Reference E1 logic):
E_xx = E3[0, 0]
E_yy = E3[1, 1]
E_xy = E3[0, 1]
rad = sqrt(((E_xx - E_yy) / 2.0)**2 + E_xy**2)
E1_expr = (E_xx + E_yy) / 2.0 + rad
Emax_expr = conditional(gt(E1_expr, 0.0), E1_expr, 0.0)

# For von Mises of PK2 stress
s3 = S3 - (1.0/3.0)*tr(S3)*I3
vmS_expr = sqrt(1.5*inner(s3, s3))

# Spaces for projection
V0 = FunctionSpace(mesh, "DG", 0)


def compute_max_Emax():
    """Project Emax to DG0 and return (Emax_function, max_value)."""
    Emax_fun = project(Emax_expr, V0)
    max_val = Emax_fun.vector().get_local().max()
    return Emax_fun, max_val


# ----------------------------
# 7) Load stepping with strain cap
# ----------------------------
u_prev = Function(V)
u_prev.assign(u)

t = 0.0
accepted = 0

print("\n--- Load stepping (cap E_max <= %.5f) ---" % Emax_cap)
while t < t_end - 1e-14 and dt >= dt_min:
    t_trial = min(t + dt, t_end)

    # Apply trial displacement multiplier (t in [0,1])
    u_R.t = t_trial

    # Good practice: start Newton from last accepted solution
    u.assign(u_prev)

    converged = True
    try:
        solver.solve()
    except RuntimeError as e:
        converged = False

    if converged:
        Emax_fun, maxE = compute_max_Emax()
        print("t_trial = %.6f | dt = %.3e | max(E_max) = %.6f" % (t_trial, dt, maxE))
    else:
        maxE = np.inf
        print("t_trial = %.6f | dt = %.3e | Newton did not converge" % (t_trial, dt))

    # Accept/reject rule
    if converged and maxE <= Emax_cap + 1e-12:
        # Accept
        t = t_trial
        u_prev.assign(u)
        accepted += 1

        # Optional: grow dt a bit after success
        dt = min(1.2*dt, dt_max)
    else:
        # Reject: do NOT advance load, reduce dt
        u.assign(u_prev)
        u_R.t = t
        dt *= 0.5

print("\nAccepted steps:", accepted)
print("Final accepted load factor t = %.6f" % t)

# Ensure u holds the accepted solution
u.assign(u_prev)
u_R.t = t

# Final fields for output
Emax_fun = project(Emax_expr, V0)
vmS_fun  = project(vmS_expr, V0)
Emax_fun.rename("E_max", "max_principal_Green_Lagrange_strain")
vmS_fun.rename("vmS", "vonMises_of_PK2")


# ----------------------------
# 8) Save figures (PNG)
# ----------------------------
# Deformed configuration plot
plt.figure(figsize=(10, 2.2))
plot(u, mode="displacement")
plt.gca().set_aspect("equal", "box")
plt.title("Deformed configuration (u applied, t=%.4f)" % t)
plt.tight_layout()
plt.savefig("q15_def.png", dpi=300)
plt.close()

# E_max colormap
plt.figure(figsize=(10, 2.2))
p = plot(Emax_fun)
plt.gca().set_aspect("equal", "box")
plt.title("Max principal Green–Lagrange strain E_max (t=%.4f)" % t)
plt.colorbar(p)
plt.tight_layout()
plt.savefig("q15_Emax.png", dpi=300)
plt.close()

# von Mises of PK2 stress colormap
plt.figure(figsize=(10, 2.2))
p = plot(vmS_fun)
plt.gca().set_aspect("equal", "box")
plt.title("von Mises of PK2 stress σ_vm(S) (t=%.4f)" % t)
plt.colorbar(p)
plt.tight_layout()
plt.savefig("q15_vmS.png", dpi=300)
plt.close()


# ----------------------------
# 9) Export XDMF (final u and E_max)
# ----------------------------
# Displacement
xdmf_u = XDMFFile(mesh.mpi_comm(), "q15_u.xdmf")
xdmf_u.parameters["flush_output"] = True
xdmf_u.write(u, 0.0)
xdmf_u.close()

# E_max
xdmf_e = XDMFFile(mesh.mpi_comm(), "q15_Emax.xdmf")
xdmf_e.parameters["flush_output"] = True
xdmf_e.write(Emax_fun, 0.0)
xdmf_e.close()

print("\nWrote:")
print("  q15_def.png")
print("  q15_Emax.png")
print("  q15_vmS.png")
print("  q15_u.xdmf (+ .h5)")
print("  q15_Emax.xdmf (+ .h5)")

# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
import json

mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

u.rename("u", "displacement field")

solution_file = XDMFFile(mesh.mpi_comm(), "solution.xdmf")
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u, 0.0)
solution_file.close()

checkpoint_file = XDMFFile(mesh.mpi_comm(), "read_checkpoint.xdmf")
checkpoint_file.write_checkpoint(u, "u", 0.0, XDMFFile.Encoding.HDF5, False)
checkpoint_file.close()

read_checkpoint_metadata = {
    "format": "FEniCS XDMFFile.write_checkpoint",
    "time": 0.0,
    "fields": [
        {
            "symbol": "u",
            "checkpoint_name": "u",
            "meaning": "displacement field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'u', -1)",
        },
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")
