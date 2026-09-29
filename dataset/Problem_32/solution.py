"""Problem 32: ALE diffusion-reaction on an expanding disk."""

from __future__ import print_function

import csv
import json
import math

from dolfin import *
from mshr import Circle, generate_mesh


class RadialVector(UserExpression):
    """Evaluate scale*x/|x|, with the vector set to zero at the disk centre."""

    def __init__(self, scale, **kwargs):
        super().__init__(**kwargs)
        self.scale = scale

    def eval(self, values, x):
        radius = math.hypot(x[0], x[1])
        if radius > 1.0e-14:
            values[0] = self.scale * x[0] / radius
            values[1] = self.scale * x[1] / radius
        else:
            values[0] = 0.0
            values[1] = 0.0

    def value_shape(self):
        return (2,)


def main():
    # -------------------------------------------------------------------------
    # Parameters.
    # -------------------------------------------------------------------------
    radius_initial = 0.05
    expansion_speed = 1.0e-3
    diffusivity_value = 1.0e-5
    kappa_value = 1.0e-4
    dt_value = 0.01
    final_time = 10.0
    num_steps = int(round(final_time / dt_value))
    mesh_resolution = 50  # Approximately h = 1.0e-3 m across the initial disk.

    comm = MPI.comm_world

    # -------------------------------------------------------------------------
    # Initial mesh, measures, and scalar P1 concentration space.
    # -------------------------------------------------------------------------
    mesh = generate_mesh(Circle(Point(0.0, 0.0), radius_initial), mesh_resolution)
    dx = Measure("dx", domain=mesh)
    ds = Measure("ds", domain=mesh)
    n = FacetNormal(mesh)
    V = FunctionSpace(mesh, "Lagrange", 1)

    c_old = Function(V)
    c_old.interpolate(Constant(1.0))
    c_new = Function(V)
    c_trial = TrialFunction(V)
    v = TestFunction(V)

    D = Constant(diffusivity_value)
    kappa = Constant(kappa_value)
    dt = Constant(dt_value)

    # The prescribed radial mesh velocity.  The point x = 0 is assigned zero
    # velocity, exactly as in the benchmark statement.
    mesh_velocity = RadialVector(scale=expansion_speed, degree=1)
    displacement = RadialVector(scale=expansion_speed * dt_value, degree=1)

    # WF-32a from weak_form.md.  Integration by parts of diffusion plus
    # (-D grad(c) - w c).n = 0 produces the final, essential w*c boundary term.
    a = (
        c_trial * v / dt * dx
        - v * dot(mesh_velocity, grad(c_trial)) * dx
        + D * dot(grad(c_trial), grad(v)) * dx
        + kappa * c_trial * v * dx
        + v * c_trial * dot(mesh_velocity, n) * ds
    )
    L = c_old * v / dt * dx

    # -------------------------------------------------------------------------
    # Mass is recorded every step.  The final concentration is written once to
    # the canonical solution.xdmf output below.
    # -------------------------------------------------------------------------

    mass_rows = []
    c_old.rename("c", "chemical concentration")
    mass_initial = assemble(c_old * dx)
    mass_rows.append((0, 0.0, mass_initial, 1.0, 0.0))
    print("step %4d, t = %.2f s, mass = %.16e" % (0, 0.0, mass_initial))

    # The canonical agent-visible mesh is the initial ALE mesh.  Write it after
    # the first form assembly (required by this legacy DOLFIN/MSHR environment)
    # but before any ALE motion, then keep it immutable.
    mesh_file = XDMFFile(comm, "mesh.xdmf")
    mesh_file.write(mesh)
    mesh_file.close()

    # -------------------------------------------------------------------------
    # Time stepping.  ALE.move applies the prescribed per-step radial
    # displacement.  The degrees of freedom of c_old remain attached to the
    # moved nodes, providing the required ALE transport of the prior field.
    # -------------------------------------------------------------------------
    for step in range(1, num_steps + 1):
        ALE.move(mesh, displacement)

        A = assemble(a)
        b = assemble(L)
        solve(A, c_new.vector(), b, "mumps")

        time_value = step * dt_value
        c_new.rename("c", "chemical concentration")
        mass = assemble(c_new * dx)
        expected_mass = mass_initial * math.exp(-kappa_value * time_value)
        relative_error = (mass - expected_mass) / mass_initial
        mass_rows.append((step, time_value, mass, expected_mass, relative_error))
        if step % 100 == 0:
            print(
                "step %4d, t = %.2f s, mass = %.16e, relative error = %.3e"
                % (step, time_value, mass, relative_error)
            )

        c_old.assign(c_new)

    # -------------------------------------------------------------------------
    # Standard dataset outputs.  The final moved geometry is embedded in both
    # solution files; mesh.xdmf remains the initial radius-0.05 input mesh.
    # -------------------------------------------------------------------------
    c_old.rename("c", "chemical concentration")
    solution_file = XDMFFile(comm, "solution.xdmf")
    solution_file.parameters["flush_output"] = True
    solution_file.parameters["functions_share_mesh"] = True
    solution_file.write(c_old, final_time)
    solution_file.close()

    checkpoint_file = XDMFFile(comm, "read_checkpoint.xdmf")
    checkpoint_file.write_checkpoint(
        c_old, "c", final_time, XDMFFile.Encoding.HDF5, False
    )
    checkpoint_file.close()

    if MPI.rank(comm) == 0:
        with open("total_mass.csv", "w", newline="", encoding="utf-8") as output_file:
            writer = csv.writer(output_file, lineterminator="\n")
            writer.writerow(
                [
                    "step",
                    "time_s",
                    "total_mass",
                    "expected_mass",
                    "relative_error_vs_expected",
                ]
            )
            writer.writerows(mass_rows)

        summary = {
            "validation_case": "ALE diffusion-reaction moving-boundary run",
            "expected_total_mass": "M(0) * exp(-kappa * t)",
            "kappa_per_s": kappa_value,
            "initial_mass": mass_initial,
            "final_mass": mass_rows[-1][2],
            "maximum_absolute_relative_mass_error": max(
                abs(row[4]) for row in mass_rows
            ),
            "num_time_steps": num_steps,
            "dt_s": dt_value,
            "final_time_s": final_time,
        }
        with open("mass_validation.json", "w", encoding="utf-8") as output_file:
            json.dump(summary, output_file, indent=2)
            output_file.write("\n")

        checkpoint_metadata = {
            "format": "FEniCS XDMFFile.write_checkpoint",
            "checkpoint_policy": "Only the final concentration at t = 10.0 s is stored.",
            "mesh_policy": (
                "mesh.xdmf stores the initial radius-0.05 ALE mesh; the final "
                "radius-0.06 geometry is embedded in the solution and checkpoint files."
            ),
            "time": final_time,
            "fields": [
                {
                    "symbol": "c_old",
                    "checkpoint_name": "c",
                    "time": final_time,
                    "meaning": "chemical concentration at the final time",
                    "read_example": (
                        "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf')"
                        ".read_checkpoint(function, 'c', -1)"
                    ),
                }
            ],
            "normal_outputs": [
                {"file": "total_mass.csv", "meaning": "total mass at every time step"},
            ],
        }
        with open("read_checkpoint.json", "w", encoding="utf-8") as output_file:
            json.dump(checkpoint_metadata, output_file, indent=2)
            output_file.write("\n")

    print("Finished Problem 32 kappa = %.1e validation run." % kappa_value)


if __name__ == "__main__":
    main()
