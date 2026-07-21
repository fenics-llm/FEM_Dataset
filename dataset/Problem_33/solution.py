# Copied from Results_ALL-FEM-main/reference solutions/multiphysics - Copy/2/multi-2.py for dataset Problem_33.
# Original benchmark: Multiphysics Problem 2 (Hard).

# Allen-Cahn curvature-flow equation on the unit square
# Legacy FEniCS / DOLFIN version

from __future__ import print_function

from dolfin import *
import math


# -----------------------------------------------------------------------------
# 1. Mesh and function space
# -----------------------------------------------------------------------------
N = 200
mesh = UnitSquareMesh(N, N)
V = FunctionSpace(mesh, "Lagrange", 1)


# -----------------------------------------------------------------------------
# 2. Model parameters
# -----------------------------------------------------------------------------
eps_val = 0.01
M_val = 1.0
dt_val = 1.0e-3
T = 0.20
num_steps = int(round(T / dt_val))

eps = Constant(eps_val)
M = Constant(M_val)
dt = Constant(dt_val)


# -----------------------------------------------------------------------------
# 3. Initial condition
#
# Signed Euclidean distance to the surface of the centered square:
# center = (0.5, 0.5), side length = 0.5, half side = 0.25.
# The sign convention is negative inside and positive outside.
# -----------------------------------------------------------------------------
class InitialPhi(UserExpression):
    def __init__(self, eps_value, **kwargs):
        super(InitialPhi, self).__init__(**kwargs)
        self.eps_value = float(eps_value)

    def eval(self, values, x):
        half_side = 0.25

        qx = abs(x[0] - 0.5) - half_side
        qy = abs(x[1] - 0.5) - half_side

        # Euclidean signed distance to an axis-aligned square.
        outside = math.sqrt(max(qx, 0.0)**2 + max(qy, 0.0)**2)
        inside = min(max(qx, qy), 0.0)
        d_rect = outside + inside

        values[0] = math.tanh(d_rect / (math.sqrt(2.0) * self.eps_value))

    def value_shape(self):
        return ()


phi_n = Function(V)
phi = Function(V)

phi_n.interpolate(InitialPhi(eps_val, degree=4))
phi_n.rename("phi", "phase field")

phi.assign(phi_n)
phi.rename("phi", "phase field")


# -----------------------------------------------------------------------------
# 4. Weak form: fully implicit backward Euler
#
# PDE:
#   phi_t = -M * ((1/eps) * W'(phi) - eps * Laplacian(phi))
#
# with
#   W'(phi) = phi^3 - phi.
#
# Equivalent strong form:
#   phi_t + (M/eps) * W'(phi) - M*eps*Laplacian(phi) = 0.
#
# Weak form after integration by parts:
#   int ((phi - phi_n)/dt) v dx
# + int (M/eps) (phi^3 - phi) v dx
# + int M eps grad(phi).grad(v) dx = 0.
#
# The boundary term vanishes because grad(phi).n = 0.
# -----------------------------------------------------------------------------
v = TestFunction(V)
dphi = TrialFunction(V)

Wprime = phi**3 - phi

F = ((phi - phi_n) / dt) * v * dx \
    + (M / eps) * Wprime * v * dx \
    + M * eps * dot(grad(phi), grad(v)) * dx

J = derivative(F, phi, dphi)

# Homogeneous Neumann boundary conditions are natural here.
bcs = []

problem = NonlinearVariationalProblem(F, phi, bcs, J)
solver = NonlinearVariationalSolver(problem)

prm = solver.parameters
prm["nonlinear_solver"] = "newton"
prm["newton_solver"]["absolute_tolerance"] = 1.0e-10
prm["newton_solver"]["relative_tolerance"] = 1.0e-8
prm["newton_solver"]["maximum_iterations"] = 25
prm["newton_solver"]["relaxation_parameter"] = 1.0

try:
    prm["newton_solver"]["error_on_nonconvergence"] = True
except Exception:
    pass


# -----------------------------------------------------------------------------
# 5. Output setup
# -----------------------------------------------------------------------------
output_steps = {
    0: "phi_t0_00.xdmf",
    int(round(0.05 / dt_val)): "phi_t0_05.xdmf",
    int(round(0.10 / dt_val)): "phi_t0_10.xdmf",
    int(round(0.20 / dt_val)): "phi_t0_20.xdmf",
}
output_step_times = {
    0: 0.00,
    int(round(0.05 / dt_val)): 0.05,
    int(round(0.10 / dt_val)): 0.10,
    int(round(0.20 / dt_val)): 0.20,
}
snapshots = {}
snapshots[0] = Function(V)
snapshots[0].assign(phi_n)

mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

phi_n.rename("phi", "phase field")

print("Saved phi_t0.00.xdmf")


# -----------------------------------------------------------------------------
# 6. Time stepping
# -----------------------------------------------------------------------------
for step in range(1, num_steps + 1):
    t = step * dt_val

    # Use previous time step as Newton initial guess.
    phi.assign(phi_n)

    solver.solve()

    if step in output_steps:
        phi.rename("phi", "phase field")
        output_time = output_step_times[step]
        snapshots[step] = Function(V)
        snapshots[step].assign(phi)
        phi.rename("phi", "phase field")
        print("Saved %s at t = %.3f" % (output_steps[step], output_time))

    phi_n.assign(phi)

for step in sorted(output_steps):
    snap = snapshots[step]
    output_time = output_step_times[step]
    snap.rename("phi", "phase field at t = %.2f" % output_time)
    phi_file = XDMFFile(mesh.mpi_comm(), output_steps[step])
    phi_file.parameters["flush_output"] = True
    phi_file.write(snap, output_time)
    phi_file.close()

phi_n.rename("phi", "phase field at t = 0.20")

solution_file = XDMFFile(mesh.mpi_comm(), "solution.xdmf")
solution_file.parameters["flush_output"] = True
solution_file.write(phi_n, 0.20)
solution_file.close()

checkpoint_file = XDMFFile(mesh.mpi_comm(), "read_checkpoint.xdmf")
checkpoint_file.write_checkpoint(phi_n, "phi", 0.20, XDMFFile.Encoding.HDF5, False)
checkpoint_file.close()

import json
read_checkpoint_metadata = {
    "format": "FEniCS XDMFFile.write_checkpoint",
    "time": 0.20,
    "fields": [
        {
            "symbol": "phi_n",
            "checkpoint_name": "phi",
            "time": 0.20,
            "meaning": "phase field at t = 0.20",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'phi', -1)",
        },
    ],
    "normal_outputs": [
        {"time": 0.00, "file": "phi_t0_00.xdmf"},
        {"time": 0.05, "file": "phi_t0_05.xdmf"},
        {"time": 0.10, "file": "phi_t0_10.xdmf"},
        {"time": 0.20, "file": "phi_t0_20.xdmf"},
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")

print("Done.")
