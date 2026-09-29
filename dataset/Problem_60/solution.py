"""Reference solution for Problem 60, adapted from source dataset entry 242."""

from __future__ import print_function

import json
import math
import os
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dolfin import *
from mshr import Circle, Rectangle, generate_mesh


problem_dir = os.path.dirname(os.path.abspath(__file__))
artifact_dir = os.path.join(problem_dir, "reference solution artifacts")
os.makedirs(artifact_dir, exist_ok=True)
start_time = time.time()
comm = MPI.comm_world

domain = Rectangle(Point(-1.0, -1.0), Point(1.0, 1.0)) - Circle(
    Point(0.5, 0.5), 0.25
)
mesh_resolution = 40
mesh = generate_mesh(domain, mesh_resolution)
V = FunctionSpace(mesh, "Lagrange", 1)
dx_mesh = Measure("dx", domain=mesh)
domain_area = assemble(Constant(1.0) * dx_mesh)


def outer_boundary(x, on_boundary):
    return on_boundary and (
        near(abs(x[0]), 1.0) or near(abs(x[1]), 1.0)
    )


def hole_boundary(x, on_boundary):
    distance = math.sqrt((x[0] - 0.5) ** 2 + (x[1] - 0.5) ** 2)
    return on_boundary and distance < 0.275


bcs = [
    DirichletBC(V, Constant(10.0), outer_boundary),
    DirichletBC(V, Constant(100.0), hole_boundary),
]

dt = 0.25
step_count = 20
k = Constant(1.0)
u = TrialFunction(V)
v = TestFunction(V)
u_old = interpolate(Constant(40.0), V)
a = u * v * dx + Constant(dt) * k * inner(grad(u), grad(v)) * dx
L = u_old * v * dx
u_new = Function(V, name="temperature")

time_history = [{"step": 0, "time": 0.0, "mean_temperature": 40.0}]
for step in range(1, step_count + 1):
    solve(
        a == L,
        u_new,
        bcs,
        solver_parameters={"linear_solver": "mumps"},
    )
    u_old.assign(u_new)
    current_time = step * dt
    time_history.append(
        {
            "step": step,
            "time": current_time,
            "mean_temperature": assemble(u_old * dx_mesh) / domain_area,
        }
    )
    if step % 2 == 0:
        plt.figure(figsize=(6, 5))
        plot(u_old)
        plt.title("Temperature at t = %.2f" % current_time)
        plt.tight_layout()
        plt.savefig(
            os.path.join(artifact_dir, "heat_implicit_%02d.png" % step),
            dpi=150,
        )
        plt.close()

u_h = Function(V, name="temperature")
u_h.assign(u_old)
u_h.rename("temperature", "temperature at t=5")

mesh_file = XDMFFile(comm, os.path.join(problem_dir, "mesh.xdmf"))
mesh_file.write(mesh)
mesh_file.close()
solution_file = XDMFFile(comm, os.path.join(problem_dir, "solution.xdmf"))
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u_h, step_count * dt)
solution_file.close()
checkpoint_file = XDMFFile(
    comm, os.path.join(problem_dir, "read_checkpoint.xdmf")
)
checkpoint_file.write_checkpoint(
    u_h,
    "temperature",
    step_count * dt,
    XDMFFile.Encoding.HDF5,
    False,
)
checkpoint_file.close()

if MPI.rank(comm) == 0:
    metadata = {
        "format": "FEniCS XDMFFile.write_checkpoint",
        "time": step_count * dt,
        "fields": [
            {
                "symbol": "u_h",
                "checkpoint_name": "temperature",
                "meaning": "temperature at t=5",
                "element_family": "Lagrange",
                "element_degree": 1,
                "value_shape": [],
            }
        ],
    }
    with open(
        os.path.join(problem_dir, "read_checkpoint.json"),
        "w",
        encoding="utf-8",
    ) as metadata_file:
        json.dump(metadata, metadata_file, indent=2)
        metadata_file.write("\n")

    diagnostics = {
        "source_dataset": "orange67/dataset_fenics_experiment_v7",
        "source_entry_one_based": 242,
        "mesh_resolution": mesh_resolution,
        "time_step": dt,
        "step_count": step_count,
        "final_time": step_count * dt,
        "cells": mesh.num_cells(),
        "vertices": mesh.num_vertices(),
        "final_minimum": u_h.vector().min(),
        "final_maximum": u_h.vector().max(),
        "time_history": time_history,
        "elapsed_seconds": time.time() - start_time,
    }
    with open(
        os.path.join(artifact_dir, "solver_diagnostics.json"),
        "w",
        encoding="utf-8",
    ) as diagnostics_file:
        json.dump(diagnostics, diagnostics_file, indent=2)
        diagnostics_file.write("\n")
    print("Problem 60 transient heat solve complete")
