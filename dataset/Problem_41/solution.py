"""Reference solution for Problem 41: high-Womersley channel flow."""

from __future__ import print_function

import json
import math
import os

from dolfin import *
from mpi4py import MPI as MPI4Py


class PeriodicBoundary(SubDomain):
    """Map the right channel boundary to the left boundary."""

    def inside(self, x, on_boundary):
        return bool(on_boundary and near(x[0], 0.0))

    def map(self, x, y):
        y[0] = x[0] - 2.0
        y[1] = x[1]


class ChannelWalls(SubDomain):
    """Identify the walls at y = +/- 0.5."""

    def inside(self, x, on_boundary):
        return bool(
            on_boundary
            and (near(x[1], -0.5, DOLFIN_EPS) or near(x[1], 0.5, DOLFIN_EPS))
        )


def main():
    problem_dir = os.path.dirname(os.path.abspath(__file__))
    comm = MPI.comm_world
    start_time = MPI4Py.Wtime()

    solution_file = XDMFFile(comm, os.path.join(problem_dir, "solution.xdmf"))
    solution_file.parameters["flush_output"] = True
    solution_file.parameters["functions_share_mesh"] = True
    checkpoint_file = XDMFFile(
        comm, os.path.join(problem_dir, "read_checkpoint.xdmf")
    )

    L = 2.0
    a_channel = 0.5
    nx = 64
    ny = 160
    rho_value = 1.0
    mu_value = 1.0
    G0 = 1.0
    omega = 400.0
    period = 2.0 * math.pi / omega
    steps_per_period = 100
    num_periods = 80
    dt_value = period / steps_per_period
    dt = Constant(dt_value)
    num_steps = num_periods * steps_per_period
    output_steps = {
        7900: (79.0 * period, "79P"),
        7925: (79.25 * period, "79P_plus_1over4"),
        7950: (79.5 * period, "79P_plus_1over2"),
        7975: (79.75 * period, "79P_plus_3over4"),
        8000: (80.0 * period, "80P"),
    }

    mesh = RectangleMesh(
        Point(0.0, -a_channel), Point(L, a_channel), nx, ny, "crossed"
    )
    periodic_boundary = PeriodicBoundary()
    wall_id = 1
    facets = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)
    ChannelWalls().mark(facets, wall_id)
    dx = Measure("dx", domain=mesh)

    mesh_file = XDMFFile(comm, os.path.join(problem_dir, "mesh.xdmf"))
    mesh_file.write(mesh)
    mesh_file.close()

    velocity_element = VectorElement("Lagrange", mesh.ufl_cell(), 2)
    pressure_element = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    real_element = FiniteElement("Real", mesh.ufl_cell(), 0)
    W = FunctionSpace(
        mesh,
        MixedElement([velocity_element, pressure_element, real_element]),
        constrained_domain=periodic_boundary,
    )
    V_output = VectorFunctionSpace(mesh, "Lagrange", 2)
    Q_output = FunctionSpace(mesh, "Lagrange", 1)

    (u, p, pressure_constraint) = TrialFunctions(W)
    (v, q, pressure_constraint_test) = TestFunctions(W)
    w = Function(W)
    w_previous = Function(W)
    rho = Constant(rho_value)
    mu = Constant(mu_value)

    # Constant Crank--Nicolson matrix.
    a_form = (
        (rho / dt) * inner(u, v) * dx
        + 0.5 * mu * inner(grad(u), grad(v)) * dx
        - p * div(v) * dx
        + q * div(u) * dx
        + pressure_constraint * q * dx
        + pressure_constraint_test * p * dx
    )

    # Constant history operator and unit-amplitude force vector.
    history_form = (
        (rho / dt) * inner(u, v) * dx
        - 0.5 * mu * inner(grad(u), grad(v)) * dx
    )
    force_form = inner(Constant((G0, 0.0)), v) * dx
    no_slip = DirichletBC(W.sub(0), Constant((0.0, 0.0)), facets, wall_id)
    bcs = [no_slip]

    A = assemble(a_form)
    history_matrix = assemble(history_form)
    force_vector = assemble(force_form)
    for bc in bcs:
        bc.apply(A)
    linear_solver = LUSolver(A, "mumps")

    checkpoint_append = False
    checkpoint_fields = []
    output_summary = []
    domain_area = assemble(Constant(1.0) * dx)

    for step in range(1, num_steps + 1):
        t = step * dt_value
        t_midpoint = t - 0.5 * dt_value
        b = history_matrix * w_previous.vector()
        b.axpy(math.cos(omega * t_midpoint), force_vector)
        for bc in bcs:
            bc.apply(b)
        linear_solver.solve(w.vector(), b)
        w_previous.assign(w)

        if step % 1000 == 0 and MPI.rank(comm) == 0:
            print("Completed step %d / %d (t/P = %.1f)." % (step, num_steps, t / period))

        if step not in output_steps:
            continue

        requested_time, time_label = output_steps[step]
        u_periodic = w.sub(0, deepcopy=True)
        p_periodic = w.sub(1, deepcopy=True)
        u_output = Function(V_output)
        p_output = Function(Q_output)
        LagrangeInterpolator.interpolate(u_output, u_periodic)
        LagrangeInterpolator.interpolate(p_output, p_periodic)
        u_output.rename("u", "velocity")
        p_output.rename("p", "pressure")

        solution_file.write(u_output, requested_time)
        solution_file.write(p_output, requested_time)

        u_checkpoint_name = "u_t" + time_label
        p_checkpoint_name = "p_t" + time_label
        checkpoint_file.write_checkpoint(
            u_output,
            u_checkpoint_name,
            requested_time,
            XDMFFile.Encoding.HDF5,
            checkpoint_append,
        )
        checkpoint_append = True
        checkpoint_file.write_checkpoint(
            p_output,
            p_checkpoint_name,
            requested_time,
            XDMFFile.Encoding.HDF5,
            True,
        )
        checkpoint_fields.extend(
            [
                {
                    "symbol": "u_output",
                    "checkpoint_name": u_checkpoint_name,
                    "time": requested_time,
                    "meaning": "velocity field",
                    "read_example": (
                        "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf')"
                        ".read_checkpoint(function, '%s', -1)"
                        % u_checkpoint_name
                    ),
                },
                {
                    "symbol": "p_output",
                    "checkpoint_name": p_checkpoint_name,
                    "time": requested_time,
                    "meaning": "mean-zero pressure field",
                    "read_example": (
                        "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf')"
                        ".read_checkpoint(function, '%s', -1)"
                        % p_checkpoint_name
                    ),
                },
            ]
        )

        mean_pressure = assemble(p_output * dx) / domain_area
        mean_velocity = assemble(u_output[0] * dx) / domain_area
        output_summary.append(
            {
                "time": requested_time,
                "periods": requested_time / period,
                "mean_streamwise_velocity": mean_velocity,
                "mean_pressure": mean_pressure,
            }
        )
        print(
            "Saved t/P = %.2f: <u_x> = %.12e, <p> = %.12e"
            % (requested_time / period, mean_velocity, mean_pressure)
        )

    solution_file.close()
    checkpoint_file.close()

    metadata = {
        "format": "FEniCS XDMFFile.write_checkpoint",
        "period": period,
        "times": [output_steps[step][0] for step in sorted(output_steps)],
        "fields": checkpoint_fields,
    }
    with open(
        os.path.join(problem_dir, "read_checkpoint.json"), "w", encoding="utf-8"
    ) as metadata_file:
        json.dump(metadata, metadata_file, indent=2)
        metadata_file.write("\n")

    elapsed = MPI4Py.Wtime() - start_time
    if MPI.rank(comm) == 0:
        print("Completed %d steps in %.3f s." % (num_steps, elapsed))
        print(json.dumps(output_summary, indent=2))


if __name__ == "__main__":
    main()
