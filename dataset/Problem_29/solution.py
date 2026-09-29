"""SUPG solution of the revised Problem 29 steady advection--diffusion case.

Provenance:
    Originally copied from
    Results_ALL-FEM-main/reference solutions/fluid/13/fluid13.py.
    Rewritten for the owner-approved outlet-layer scale, symmetric mesh,
    diffusion-aware SUPG parameter, canonical outputs, and diagnostics.
"""

from __future__ import print_function

import json
import math
import os
import time

from dolfin import *


def main():
    start_time = time.time()
    comm = MPI.comm_world
    rank = MPI.rank(comm)

    # Physical parameters and owner-selected outlet-layer target.
    L = 1.0
    H = 0.10
    Umax = 1.0e-2
    diffusivity = 1.0e-5
    nx = 440
    ny = 44

    hx = L / float(nx)
    hy = H / float(ny)
    delta = diffusivity / Umax
    global_peclet = Umax * L / diffusivity
    streamwise_cell_peclet = Umax * hx / (2.0 * diffusivity)
    visible_transition_width = math.log(100.0) * delta
    transition_intervals = visible_transition_width / hx

    # Crossed triangles preserve reflection symmetry about y=H/2.
    mesh = RectangleMesh(
        Point(0.0, 0.0), Point(L, H), nx, ny, diagonal="crossed"
    )
    dx_measure = Measure("dx", domain=mesh)

    V = FunctionSpace(mesh, "CG", 1)
    c = TrialFunction(V)
    v = TestFunction(V)

    velocity = Expression(
        ("4.0*Umax*x[1]*(H-x[1])/(H*H)", "0.0"),
        Umax=Umax,
        H=H,
        degree=2,
    )
    velocity_vector = as_vector((velocity[0], velocity[1]))

    class Inlet(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], 0.0)

    class Outlet(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], L)

    bcs = [
        DirichletBC(V, Constant(0.0), Inlet()),
        DirichletBC(V, Constant(1.0), Outlet()),
    ]

    # Weak form for the revised Problem 29 model.
    a_galerkin = (
        dot(velocity_vector, grad(c)) * v
        + diffusivity * dot(grad(c), grad(v))
    ) * dx_measure
    linear_form = Constant(0.0) * v * dx_measure

    # Diffusion-aware SUPG parameter. It approaches h/(2|u|) for high Pe and
    # h^2/(12D) for low Pe, without a singularity where the wall velocity is 0.
    cell_size = CellDiameter(mesh)
    speed = sqrt(dot(velocity_vector, velocity_vector))
    tau = 1.0 / sqrt(
        (2.0 * speed / cell_size) ** 2
        + (12.0 * diffusivity / cell_size ** 2) ** 2
    )
    strong_residual = dot(velocity_vector, grad(c)) - diffusivity * div(grad(c))
    a_supg = tau * dot(velocity_vector, grad(v)) * strong_residual * dx_measure
    a = a_galerkin + a_supg

    A = assemble(a)
    b = assemble(linear_form)
    for bc in bcs:
        bc.apply(A, b)

    c_h = Function(V, name="c")
    solve(A, c_h.vector(), b, "mumps")
    c_h.rename("c", "concentration field")

    mesh_file = XDMFFile(comm, "mesh.xdmf")
    mesh_file.write(mesh)
    mesh_file.close()

    solution_file = XDMFFile(comm, "solution.xdmf")
    solution_file.parameters["flush_output"] = True
    solution_file.parameters["functions_share_mesh"] = True
    solution_file.write(c_h, 0.0)
    solution_file.close()

    checkpoint_file = XDMFFile(comm, "read_checkpoint.xdmf")
    checkpoint_file.write_checkpoint(
        c_h, "c", 0.0, XDMFFile.Encoding.HDF5, False
    )
    checkpoint_file.close()

    checkpoint_metadata = {
        "format": "FEniCS XDMFFile.write_checkpoint",
        "time": 0.0,
        "fields": [
            {
                "symbol": "c_h",
                "checkpoint_name": "c",
                "meaning": "concentration field",
                "read_example": (
                    "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf')"
                    ".read_checkpoint(function, 'c', -1)"
                ),
            }
        ],
    }
    if rank == 0:
        with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
            json.dump(checkpoint_metadata, metadata_file, indent=2)
            metadata_file.write("\n")

        artifact_dir = "reference solution artifacts"
        os.makedirs(artifact_dir, exist_ok=True)
        diagnostics = {
            "maximum_velocity_m_per_s": Umax,
            "diffusivity_m2_per_s": diffusivity,
            "mesh_subdivisions": [nx, ny],
            "streamwise_spacing_m": hx,
            "transverse_spacing_m": hy,
            "global_Peclet_number": global_peclet,
            "streamwise_cell_Peclet_number": streamwise_cell_peclet,
            "one_dimensional_e_folding_layer_m": delta,
            "layer_fraction_of_channel_length": delta / L,
            "one_percent_to_outlet_transition_width_m": visible_transition_width,
            "streamwise_intervals_across_transition": transition_intervals,
            "SUPG_enabled": True,
            "SUPG_tau": "[(2|u|/h_K)^2 + (12D/h_K^2)^2]^(-1/2)",
            "elapsed_seconds": time.time() - start_time,
        }
        with open(
            os.path.join(artifact_dir, "solver_diagnostics.json"),
            "w",
            encoding="utf-8",
        ) as output_file:
            json.dump(diagnostics, output_file, indent=2)
            output_file.write("\n")

        values = c_h.vector().get_local()
        print("Problem 29 revised SUPG solve complete")
        print("  mesh subdivisions: {} x {}".format(nx, ny))
        print("  Pe_L: {:.6f}".format(global_peclet))
        print("  streamwise Pe_h: {:.6f}".format(streamwise_cell_peclet))
        print("  delta/L: {:.6e}".format(delta / L))
        print("  transition intervals: {:.6f}".format(transition_intervals))
        print(
            "  nodal concentration range: [{:.12e}, {:.12e}]".format(
                values.min(), values.max()
            )
        )


if __name__ == "__main__":
    main()
