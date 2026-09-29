"""Reference solution for Problem 40: low-Womersley oscillatory channel flow."""

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
    """Identify the no-slip walls at y = +/- 0.5."""

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
    nx = 128
    ny = 64
    rho = Constant(1.0)
    mu = Constant(1.0)
    G0 = 1.0
    omega = 1.0
    dt_value = math.pi / 100.0
    dt = Constant(dt_value)
    num_steps = 600
    output_steps = {
        400: (4.0 * math.pi, "4pi"),
        450: (4.5 * math.pi, "9pi_over_2"),
        500: (5.0 * math.pi, "5pi"),
        550: (5.5 * math.pi, "11pi_over_2"),
        600: (6.0 * math.pi, "6pi"),
    }

    mesh = RectangleMesh(
        Point(0.0, -a_channel), Point(L, a_channel), nx, ny, "crossed"
    )
    periodic_boundary = PeriodicBoundary()
    wall_id = 1
    facets = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)
    ChannelWalls().mark(facets, wall_id)
    dx = Measure("dx", domain=mesh)
    ds = Measure("ds", domain=mesh, subdomain_data=facets)

    mesh_file = XDMFFile(comm, os.path.join(problem_dir, "mesh.xdmf"))
    mesh_file.write(mesh)
    mesh_file.close()

    velocity_element = VectorElement("Lagrange", mesh.ufl_cell(), 2)
    pressure_element = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    real_element = FiniteElement("Real", mesh.ufl_cell(), 0)
    mixed_element = MixedElement(
        [velocity_element, pressure_element, real_element]
    )
    W = FunctionSpace(mesh, mixed_element, constrained_domain=periodic_boundary)
    V = VectorFunctionSpace(
        mesh, "Lagrange", 2, constrained_domain=periodic_boundary
    )
    V_output = VectorFunctionSpace(mesh, "Lagrange", 2)
    Q_output = FunctionSpace(mesh, "Lagrange", 1)

    (u, p, pressure_constraint) = TrialFunctions(W)
    (v, q, pressure_constraint_test) = TestFunctions(W)
    w = Function(W)
    u_previous = Function(V)
    u_midpoint = 0.5 * (u + u_previous)
    body_force = Expression(
        ("G0*cos(omega*t)", "0.0"),
        degree=2,
        G0=G0,
        omega=omega,
        t=0.0,
    )

    # Crank--Nicolson weak form of the question's Laplacian Stokes operator.
    F = (
        rho * inner((u - u_previous) / dt, v) * dx
        + mu * inner(grad(u_midpoint), grad(v)) * dx
        - p * div(v) * dx
        + q * div(u) * dx
        + pressure_constraint * q * dx
        + pressure_constraint_test * p * dx
        - inner(body_force, v) * dx
    )
    a_form = lhs(F)
    L_form = rhs(F)
    no_slip = DirichletBC(W.sub(0), Constant((0.0, 0.0)), facets, wall_id)
    bcs = [no_slip]

    A = assemble(a_form)
    for bc in bcs:
        bc.apply(A)
    linear_solver = LUSolver(A, "mumps")

    checkpoint_append = False
    checkpoint_fields = []
    output_summary = []
    domain_area = assemble(Constant(1.0) * dx)

    for step in range(1, num_steps + 1):
        t = step * dt_value
        body_force.t = t - 0.5 * dt_value
        b = assemble(L_form)
        for bc in bcs:
            bc.apply(b)
        linear_solver.solve(w.vector(), b)

        u_periodic = w.sub(0, deepcopy=True)
        u_previous.assign(u_periodic)

        if step not in output_steps:
            continue

        requested_time, time_label = output_steps[step]
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
                "mean_streamwise_velocity": mean_velocity,
                "mean_pressure": mean_pressure,
            }
        )
        print(
            "Saved t = %.12g: <u_x> = %.12e, <p> = %.12e"
            % (requested_time, mean_velocity, mean_pressure)
        )

    solution_file.close()
    checkpoint_file.close()

    metadata = {
        "format": "FEniCS XDMFFile.write_checkpoint",
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
        print("Completed %d Crank--Nicolson steps in %.3f s." % (num_steps, elapsed))
        print(json.dumps(output_summary, indent=2))


if __name__ == "__main__":
    main()
