"""Reference implementation for dataset Problem_34.

Mixed, fully implicit P1 Cahn--Hilliard solve with periodic concentration and
chemical-potential fields.  The formulation is documented in
weak_form_and_implementation.md.
"""

from __future__ import print_function

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from dolfin import *


class PeriodicBoundary(SubDomain):
    """Map right/top unit-square boundaries to left/bottom boundaries."""

    def inside(self, x, on_boundary):
        return bool(
            on_boundary
            and (near(x[0], 0.0) or near(x[1], 0.0))
            and not (
                (near(x[0], 0.0) and near(x[1], 1.0))
                or (near(x[0], 1.0) and near(x[1], 0.0))
            )
        )

    def map(self, x, y):
        if near(x[0], 1.0) and near(x[1], 1.0):
            y[0] = 0.0
            y[1] = 0.0
        elif near(x[0], 1.0):
            y[0] = 0.0
            y[1] = x[1]
        else:
            y[0] = x[0]
            y[1] = 0.0


def periodic_trace_error(field, sample_count):
    """Return the largest opposite-side trace difference for a scalar field."""

    maximum = 0.0
    for index in range(sample_count + 1):
        coordinate = float(index) / sample_count
        maximum = max(
            maximum,
            abs(float(field(Point(0.0, coordinate))) - float(field(Point(1.0, coordinate)))),
            abs(float(field(Point(coordinate, 0.0))) - float(field(Point(coordinate, 1.0)))),
        )
    return maximum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--final-time",
        type=float,
        default=1.0e-4,
        help="End time for this run (default: paper-comparison target, 1e-4).",
    )
    args = parser.parse_args()

    case_directory = Path(__file__).resolve().parent
    validation_directory = case_directory / "validation"
    validation_directory.mkdir(exist_ok=True)

    comm = MPI.comm_world
    # INFO is required during development so DOLFIN prints the absolute and
    # relative nonlinear residual at every Newton iteration.
    set_log_level(LogLevel.INFO)

    # Model and time-step controls.  The mesh is an initial baseline because
    # the benchmark does not prescribe a resolution.
    mesh_resolution = 128
    random_seed = 340063
    theta_value = 1.5
    alpha_value = 3000.0
    final_time = args.final_time
    initial_dt = 1.0e-7
    minimum_dt = 1.0e-9
    maximum_dt = 1.0e-5
    requested_output_times = [0.0, 3.0e-6, 1.0e-4]
    output_times = [time for time in requested_output_times if time <= final_time + 1.0e-14]
    if not output_times or abs(output_times[-1] - final_time) > 1.0e-14:
        output_times.append(final_time)
    record_every_accepted_step = True
    trace_sample_count = 64

    if final_time <= 0.0 or final_time > requested_output_times[-1]:
        raise ValueError("--final-time must be in (0, 1e-4].")
    if MPI.size(comm) != 1:
        raise RuntimeError("Problem 34 currently supports serial execution only.")

    mesh = UnitSquareMesh(mesh_resolution, mesh_resolution)
    dx = Measure("dx", domain=mesh)

    periodic_boundary = PeriodicBoundary()
    P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    W = FunctionSpace(
        mesh,
        MixedElement([P1, P1]),
        constrained_domain=periodic_boundary,
    )
    V = FunctionSpace(mesh, "Lagrange", 1, constrained_domain=periodic_boundary)

    expected_periodic_dimension = mesh_resolution**2
    if V.dim() != expected_periodic_dimension:
        raise RuntimeError(
            "Periodic P1 space has dimension %d; expected %d."
            % (V.dim(), expected_periodic_dimension)
        )
    print(
        "Verified periodic P1 dimension: %d" % V.dim(),
        flush=True,
    )

    w = Function(W)
    w_previous = Function(W)
    dw = TrialFunction(W)
    q, v = TestFunctions(W)
    c, mu = split(w)
    c_previous, _ = split(w_previous)

    # Seeded uniform perturbation.  Correct its FE integral to zero mean, then
    # rescale if necessary so the recorded perturbation remains in [-0.05, 0.05].
    random_generator = np.random.RandomState(random_seed)
    c_initial = Function(V)
    c_initial.vector()[:] = 0.63 + random_generator.uniform(
        -0.05, 0.05, V.dim()
    )
    c_initial.vector().apply("insert")
    initial_mean = assemble(c_initial * dx) / assemble(Constant(1.0) * dx)
    c_initial.vector()[:] -= initial_mean - 0.63
    c_initial.vector().apply("insert")
    maximum_perturbation = max(
        abs(float(c_initial.vector().min()) - 0.63),
        abs(float(c_initial.vector().max()) - 0.63),
    )
    if maximum_perturbation > 0.05:
        c_initial.vector()[:] = 0.63 + (0.05 / maximum_perturbation) * (
            c_initial.vector().get_local() - 0.63
        )
        c_initial.vector().apply("insert")
    assign(w.sub(0), c_initial)
    assign(w_previous.sub(0), c_initial)

    # A consistent-enough potential initial guess helps the first Newton step;
    # the coupled residual determines the converged chemical potential.
    c_for_mu = variable(c_initial)
    mu_c_initial = (1.0 / (2.0 * theta_value)) * ln(
        c_for_mu / (1.0 - c_for_mu)
    ) + 1.0 - 2.0 * c_for_mu
    mu_initial = project(3.0 * alpha_value * mu_c_initial, V, solver_type="mumps")
    assign(w.sub(1), mu_initial)
    assign(w_previous.sub(1), mu_initial)

    dt = Constant(initial_dt)
    mobility = c * (1.0 - c)
    mu_c = (1.0 / (2.0 * theta_value)) * ln(c / (1.0 - c)) + 1.0 - 2.0 * c

    # Weak form for weak_form_and_implementation.md Eq. (4):
    # (c-c_n, q)/dt + (M(c) grad(mu), grad(q)) = 0
    # (mu, v) - (3 alpha mu_c(c), v) - (grad(c), grad(v)) = 0.
    residual = (
        ((c - c_previous) / dt) * q * dx
        + dot(mobility * grad(mu), grad(q)) * dx
        + mu * v * dx
        - 3.0 * alpha_value * mu_c * v * dx
        - dot(grad(c), grad(v)) * dx
    )
    jacobian = derivative(residual, w, dw)

    nonlinear_problem = NonlinearVariationalProblem(residual, w, [], jacobian)
    nonlinear_solver = NonlinearVariationalSolver(nonlinear_problem)
    solver_parameters = nonlinear_solver.parameters["newton_solver"]
    solver_parameters["linear_solver"] = "mumps"
    solver_parameters["absolute_tolerance"] = 1.0e-10
    solver_parameters["relative_tolerance"] = 1.0e-8
    solver_parameters["maximum_iterations"] = 20
    solver_parameters["report"] = True
    solver_parameters["error_on_nonconvergence"] = True

    mesh_file = XDMFFile(comm, str(case_directory / "mesh.xdmf"))
    mesh_file.write(mesh)
    mesh_file.close()

    benchmark_file = XDMFFile(comm, str(case_directory / "cahn_hilliard.xdmf"))
    benchmark_file.parameters["flush_output"] = True
    benchmark_file.parameters["functions_share_mesh"] = True

    initial_mass = float(assemble(c_initial * dx))
    initial_energy = float(
        assemble(
            (
                c_initial * ln(c_initial)
                + (1.0 - c_initial) * ln(1.0 - c_initial)
                + 2.0 * theta_value * c_initial * (1.0 - c_initial)
                + (theta_value / (3.0 * alpha_value))
                * dot(grad(c_initial), grad(c_initial))
            )
            * dx
        )
    )
    mass_history = []
    snapshot_diagnostics = []
    mass_history_fieldnames = [
        "step",
        "time",
        "dt",
        "newton_iterations",
        "mass",
        "mass_error",
        "free_energy",
        "energy_change_from_previous_step",
        "concentration_minimum",
        "concentration_maximum",
        "concentration_periodic_trace_error",
        "chemical_potential_periodic_trace_error",
    ]
    mass_history_file = (validation_directory / "mass_history.csv").open(
        "w", newline="", encoding="utf-8"
    )
    mass_history_writer = csv.DictWriter(
        mass_history_file, fieldnames=mass_history_fieldnames
    )
    mass_history_writer.writeheader()
    mass_history_file.flush()

    def record_snapshot(time_value):
        c_output, mu_output = w.split(deepcopy=True)
        c_output.rename("c", "concentration")
        mu_output.rename("mu", "chemical potential")
        benchmark_file.write(c_output, time_value)
        benchmark_file.write(mu_output, time_value)
        snapshot_diagnostics.append(
            {
                "time": time_value,
                "mass": float(assemble(c_output * dx)),
                "concentration_minimum": float(c_output.vector().min()),
                "concentration_maximum": float(c_output.vector().max()),
                "concentration_periodic_trace_error": periodic_trace_error(
                    c_output, trace_sample_count
                ),
                "chemical_potential_periodic_trace_error": periodic_trace_error(
                    mu_output, trace_sample_count
                ),
            }
        )
        print("Saved snapshot at t = %.8g" % time_value, flush=True)

    record_snapshot(0.0)
    time_value = 0.0
    dt_value = initial_dt
    previous_energy = initial_energy
    output_index = 1
    accepted_steps = 0
    rejected_steps = 0

    while time_value < final_time - 1.0e-14:
        next_output_time = output_times[output_index]
        proposed_dt = dt_value
        dt_trial = min(proposed_dt, next_output_time - time_value)
        if dt_trial < minimum_dt:
            raise RuntimeError("Adaptive time step fell below the configured minimum.")

        w.assign(w_previous)
        dt.assign(dt_trial)
        print(
            "Starting step %d: t = %.8g -> %.8g (dt = %.3e)"
            % (accepted_steps + 1, time_value, time_value + dt_trial, dt_trial),
            flush=True,
        )
        try:
            newton_iterations, converged = nonlinear_solver.solve()
            c_trial, _ = w.split(deepcopy=True)
            c_minimum = float(c_trial.vector().min())
            c_maximum = float(c_trial.vector().max())
            accepted = converged and c_minimum > 0.0 and c_maximum < 1.0
        except RuntimeError:
            accepted = False
            newton_iterations = solver_parameters["maximum_iterations"]

        if not accepted:
            rejected_steps += 1
            dt_value = 0.5 * dt_trial
            w.assign(w_previous)
            print(
                "Rejected step at t = %.8g; retrying with dt = %.3e"
                % (time_value, dt_value),
                flush=True,
            )
            continue

        time_value += dt_trial
        accepted_steps += 1
        w_previous.assign(w)
        mass = float(assemble(c_trial * dx))
        free_energy = float(
            assemble(
                (
                    c_trial * ln(c_trial)
                    + (1.0 - c_trial) * ln(1.0 - c_trial)
                    + 2.0 * theta_value * c_trial * (1.0 - c_trial)
                    + (theta_value / (3.0 * alpha_value))
                    * dot(grad(c_trial), grad(c_trial))
                )
                * dx
            )
        )
        c_trace_error = periodic_trace_error(c_trial, trace_sample_count)
        _, mu_trial = w.split(deepcopy=True)
        mu_trace_error = periodic_trace_error(mu_trial, trace_sample_count)
        mass_history_row = {
            "step": accepted_steps,
            "time": time_value,
            "dt": dt_trial,
            "newton_iterations": int(newton_iterations),
            "mass": mass,
            "mass_error": mass - initial_mass,
            "free_energy": free_energy,
            "energy_change_from_previous_step": free_energy - previous_energy,
            "concentration_minimum": c_minimum,
            "concentration_maximum": c_maximum,
            "concentration_periodic_trace_error": c_trace_error,
            "chemical_potential_periodic_trace_error": mu_trace_error,
        }
        previous_energy = free_energy
        mass_history.append(mass_history_row)
        mass_history_writer.writerow(mass_history_row)
        mass_history_file.flush()
        if accepted_steps % 25 == 0:
            print(
                "Accepted step %d at t = %.8g (dt = %.3e, Newton = %d)"
                % (accepted_steps, time_value, dt_trial, newton_iterations),
                flush=True,
            )

        if record_every_accepted_step:
            record_snapshot(time_value)

        if abs(time_value - next_output_time) <= 1.0e-13:
            time_value = next_output_time
            if not record_every_accepted_step:
                record_snapshot(time_value)
            output_index += 1

        # Development controller: Newton convergence controls cost/robustness;
        # exact requested output times always cap the proposed step above.
        if newton_iterations <= 3:
            dt_value = min(1.25 * proposed_dt, maximum_dt)
        elif newton_iterations == 4:
            dt_value = proposed_dt
        else:
            dt_value = max(0.5 * proposed_dt, minimum_dt)

    benchmark_file.close()
    mass_history_file.close()

    c_final, mu_final = w.split(deepcopy=True)
    c_final.rename("c", "concentration")
    mu_final.rename("mu", "chemical potential")
    solution_file = XDMFFile(comm, str(case_directory / "solution.xdmf"))
    solution_file.parameters["flush_output"] = True
    solution_file.parameters["functions_share_mesh"] = True
    solution_file.write(c_final, final_time)
    solution_file.write(mu_final, final_time)
    solution_file.close()

    checkpoint_file = XDMFFile(comm, str(case_directory / "read_checkpoint.xdmf"))
    checkpoint_file.write_checkpoint(c_final, "c", final_time, XDMFFile.Encoding.HDF5, False)
    checkpoint_file.write_checkpoint(mu_final, "mu", final_time, XDMFFile.Encoding.HDF5, True)
    checkpoint_file.close()

    with (case_directory / "read_checkpoint.json").open("w", encoding="utf-8") as output_file:
        json.dump(
            {
                "format": "FEniCS XDMFFile.write_checkpoint",
                "time": final_time,
                "fields": [
                    {
                        "symbol": "c_final",
                        "checkpoint_name": "c",
                        "time": final_time,
                        "meaning": "concentration at the configured final time",
                    },
                    {
                        "symbol": "mu_final",
                        "checkpoint_name": "mu",
                        "time": final_time,
                        "meaning": "chemical potential at the configured final time",
                    },
                ],
                "normal_outputs": [
                    {"file": "cahn_hilliard.xdmf", "time": item}
                    for item in output_times
                ],
            },
            output_file,
            indent=2,
        )
        output_file.write("\n")

    maximum_mass_error = max(abs(item["mass_error"]) for item in mass_history)
    maximum_c_trace_error = max(
        item["concentration_periodic_trace_error"] for item in mass_history
    )
    maximum_mu_trace_error = max(
        item["chemical_potential_periodic_trace_error"] for item in mass_history
    )
    validation = {
        "problem": "Problem_34",
        "mesh": {
            "type": "UnitSquareMesh",
            "resolution": mesh_resolution,
            "periodic_scalar_dimension": V.dim(),
        },
        "initial_condition": {
            "distribution": "seeded uniform perturbation",
            "random_seed": random_seed,
            "mean": initial_mass,
            "perturbation_bounds": [-0.05, 0.05],
        },
        "time": {
            "scheme": "backward Euler with adaptive time step",
            "initial_dt": initial_dt,
            "minimum_dt": minimum_dt,
            "maximum_dt": maximum_dt,
            "final_time": final_time,
            "record_every_accepted_step": record_every_accepted_step,
            "accepted_steps": accepted_steps,
            "rejected_steps": rejected_steps,
        },
        "mass_conservation": {
            "initial_mass": initial_mass,
            "maximum_absolute_mass_error": maximum_mass_error,
            "relative_mass_error": maximum_mass_error / abs(initial_mass),
            "criterion_relative_error_at_most": 1.0e-8,
        },
        "free_energy": {
            "definition": "paper Eq. (11)",
            "initial": initial_energy,
            "final": mass_history[-1]["free_energy"],
            "maximum_stepwise_increase": max(
                0.0,
                max(item["energy_change_from_previous_step"] for item in mass_history),
            ),
        },
        "periodic_fields": {
            "sample_points_per_boundary": trace_sample_count + 1,
            "maximum_concentration_trace_error": maximum_c_trace_error,
            "maximum_chemical_potential_trace_error": maximum_mu_trace_error,
            "criterion_absolute_error_at_most": 1.0e-10,
        },
        "snapshots": snapshot_diagnostics,
    }
    validation["passed"] = bool(
        validation["mass_conservation"]["relative_mass_error"] <= 1.0e-8
        and maximum_c_trace_error <= 1.0e-10
        and maximum_mu_trace_error <= 1.0e-10
    )
    with (validation_directory / "mass_and_periodicity_validation.json").open(
        "w", encoding="utf-8"
    ) as output_file:
        json.dump(validation, output_file, indent=2)
        output_file.write("\n")

    print("Accepted steps: %d; rejected steps: %d" % (accepted_steps, rejected_steps))
    print("Maximum absolute mass error: %.6e" % maximum_mass_error)
    print("Maximum concentration periodic trace error: %.6e" % maximum_c_trace_error)
    print("Maximum chemical-potential periodic trace error: %.6e" % maximum_mu_trace_error)
    print("Mass/periodicity validation passed: %s" % validation["passed"])


if __name__ == "__main__":
    main()
