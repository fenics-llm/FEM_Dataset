# Copied from Results_ALL-FEM-main/reference solutions/fluid/15/fluid15.py for dataset Problem_31.
# Original benchmark: Fluid Mechanics Problem 15 (Medium).

# filename: navier_stokes_periodic.py
from __future__ import print_function
from dolfin import *
import numpy as np

# ----------------------------------------------------------------------
# 1. Periodic boundary definition
# ----------------------------------------------------------------------
class PeriodicBoundary(SubDomain):
    """
    Periodic on unit square:
      x = 1 maps to x = 0
      y = 1 maps to y = 0
      corner (1,1) maps to (0,0)
    """

    def inside(self, x, on_boundary):
        return bool(
            on_boundary and
            (
                near(x[0], 0.0) or near(x[1], 0.0)
            ) and
            not (
                near(x[0], 1.0) or near(x[1], 1.0)
            )
        )

    def map(self, x, y):
        if near(x[0], 1.0) and near(x[1], 1.0):
            y[0] = x[0] - 1.0
            y[1] = x[1] - 1.0
        elif near(x[0], 1.0):
            y[0] = x[0] - 1.0
            y[1] = x[1]
        elif near(x[1], 1.0):
            y[0] = x[0]
            y[1] = x[1] - 1.0
        else:
            y[0] = x[0]
            y[1] = x[1]

# ----------------------------------------------------------------------
# 2. Mesh and periodic function space
# ----------------------------------------------------------------------
N = 32
mesh = UnitSquareMesh(N, N, "crossed")
pbc = PeriodicBoundary()

V_el = VectorElement("CG", mesh.ufl_cell(), 2)
Q_el = FiniteElement("CG", mesh.ufl_cell(), 1)
W_el = MixedElement([V_el, Q_el])

W = FunctionSpace(mesh, W_el, constrained_domain=pbc)
V = FunctionSpace(mesh, V_el, constrained_domain=pbc)
Q = FunctionSpace(mesh, Q_el, constrained_domain=pbc)

# ----------------------------------------------------------------------
# 3. Unknowns and tests
# ----------------------------------------------------------------------
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

w = Function(W)
w0 = Function(W)

u0, p0 = split(w0)

# ----------------------------------------------------------------------
# 4. Parameters
# ----------------------------------------------------------------------
rho = Constant(1.0)
nu = Constant(1.0e-3)

# ----------------------------------------------------------------------
# 5. Initial condition
# ----------------------------------------------------------------------
u0_expr = Expression(
    (
        "sin(2*pi*x[0])*cos(2*pi*x[1])",
        "-cos(2*pi*x[0])*sin(2*pi*x[1])"
    ),
    degree=5,
    pi=np.pi
)

u_init = interpolate(u0_expr, V)
p_init = interpolate(Constant(0.0), Q)

assign(w0.sub(0), u_init)
assign(w0.sub(1), p_init)

# ----------------------------------------------------------------------
# 6. Time stepping
# ----------------------------------------------------------------------
T = 1.0
dt = 0.0025
num_steps = int(round(T/dt))

output_times = [0.0, 0.25, 0.5, 1.0]
output_tol = 0.5*dt + 1.0e-12
checkpoint_names = {
    0.0: "u_t0_00",
    0.25: "u_t0_25",
    0.5: "u_t0_50",
    1.0: "u_t1_00",
}

# ----------------------------------------------------------------------
# 7. Pressure gauge
# ----------------------------------------------------------------------
class PressureGauge(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0) and near(x[1], 0.0)

bc_p = DirichletBC(W.sub(1), Constant(0.0), PressureGauge(), method="pointwise")
bcs = [bc_p]

# ----------------------------------------------------------------------
# 8. Weak form: semi-implicit backward Euler
#
# Unknown u = u^{n+1}
# Previous velocity u0 = u^n
# Convection linearized as (u0 · grad) u
# ----------------------------------------------------------------------
def epsilon(a):
    return sym(grad(a))

F = (
    rho*dot((u - u0)/dt, v)*dx
    + rho*dot(dot(u0, nabla_grad(u)), v)*dx
    + 2.0*nu*inner(epsilon(u), epsilon(v))*dx
    - p*div(v)*dx
    + q*div(u)*dx
)

a = lhs(F)
Lform = rhs(F)

# ----------------------------------------------------------------------
# 9. Output file
# ----------------------------------------------------------------------
u_init.rename(checkpoint_names[0.0], "velocity field at t = 0.0")

xdmf = XDMFFile(mesh.mpi_comm(), "velocity_periodic.xdmf")
xdmf.parameters["flush_output"] = True
xdmf.parameters["functions_share_mesh"] = True
xdmf.write(u_init, 0.0)

mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

solution_file = XDMFFile(mesh.mpi_comm(), "solution.xdmf")
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u_init, 0.0)

checkpoint_file = XDMFFile(mesh.mpi_comm(), "read_checkpoint.xdmf")
checkpoint_file.write_checkpoint(u_init, checkpoint_names[0.0], 0.0, XDMFFile.Encoding.HDF5, False)
checkpoint_append = True

# ----------------------------------------------------------------------
# 10. Time loop
# ----------------------------------------------------------------------
t = 0.0

for step in range(1, num_steps + 1):
    t = step*dt

    solve(
        a == Lform,
        w,
        bcs,
        solver_parameters={"linear_solver": "mumps"}
    )

    u_sol, p_sol = w.split(deepcopy=True)

    assign(w0.sub(0), u_sol)
    assign(w0.sub(1), p_sol)

    matched_output_time = None
    for tout in output_times:
        if abs(t - tout) <= output_tol:
            matched_output_time = tout
            break

    if matched_output_time is not None:
        print("Saving velocity at t = %.4f" % matched_output_time)
        u_sol.rename(
            checkpoint_names[matched_output_time],
            "velocity field at t = %.2f" % matched_output_time,
        )
        xdmf.write(u_sol, 0.0)
        solution_file.write(u_sol, 0.0)
        checkpoint_file.write_checkpoint(
            u_sol,
            checkpoint_names[matched_output_time],
            0.0,
            XDMFFile.Encoding.HDF5,
            checkpoint_append,
        )

xdmf.close()
solution_file.close()
checkpoint_file.close()

import json
read_checkpoint_metadata = {
    "format": "FEniCS XDMFFile.write_checkpoint",
    "times": output_times,
    "fields": [
        {
            "symbol": "u_init/u_sol",
            "checkpoint_name": checkpoint_names[0.0],
            "time": 0.0,
            "meaning": "velocity field at t = 0.0",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'u_t0_00', -1)",
        },
        {
            "symbol": "u_sol",
            "checkpoint_name": checkpoint_names[0.25],
            "time": 0.25,
            "meaning": "velocity field at t = 0.25",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'u_t0_25', -1)",
        },
        {
            "symbol": "u_sol",
            "checkpoint_name": checkpoint_names[0.5],
            "time": 0.5,
            "meaning": "velocity field at t = 0.5",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'u_t0_50', -1)",
        },
        {
            "symbol": "u_sol",
            "checkpoint_name": checkpoint_names[1.0],
            "time": 1.0,
            "meaning": "velocity field at t = 1.0",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'u_t1_00', -1)",
        },
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")

print("Simulation finished.")
