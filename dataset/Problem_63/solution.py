"""Reference solution for Problem 63, adapted from source dataset entry 257."""

from __future__ import print_function

import json
import os
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dolfin import *


problem_dir = os.path.dirname(os.path.abspath(__file__))
artifact_dir = os.path.join(problem_dir, "reference solution artifacts")
os.makedirs(artifact_dir, exist_ok=True)
start_time = time.time()
comm = MPI.comm_world

mesh = RectangleMesh(Point(0.0, 0.0), Point(5.0, 5.0), 50, 50)
V = FunctionSpace(mesh, "CG", 1)


def left_boundary(x, on_boundary):
    return on_boundary and near(x[0], 0.0)


def top_bottom_boundary(x, on_boundary):
    return on_boundary and (near(x[1], 0.0) or near(x[1], 5.0))


bcs = [
    DirichletBC(
        V,
        Expression("100*x[1]*(5.0-x[1])", degree=3),
        left_boundary,
    ),
    DirichletBC(V, Constant(0.0), top_bottom_boundary),
]

u = TrialFunction(V)
v = TestFunction(V)
k = Expression("1.0 + 100*x[1]*(5.0-x[1])", degree=3)
a = k * inner(grad(u), grad(v)) * dx
L = Constant(0.0) * v * dx

u_h = Function(V, name="temperature")
solve(a == L, u_h, bcs, solver_parameters={"linear_solver": "mumps"})
u_h.rename("temperature", "variable-diffusivity temperature")

mesh_file = XDMFFile(comm, os.path.join(problem_dir, "mesh.xdmf"))
mesh_file.write(mesh)
mesh_file.close()
solution_file = XDMFFile(comm, os.path.join(problem_dir, "solution.xdmf"))
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u_h, 0.0)
solution_file.close()
checkpoint_file = XDMFFile(
    comm, os.path.join(problem_dir, "read_checkpoint.xdmf")
)
checkpoint_file.write_checkpoint(
    u_h, "temperature", 0.0, XDMFFile.Encoding.HDF5, False
)
checkpoint_file.close()

if MPI.rank(comm) == 0:
    metadata = {
        "format": "FEniCS XDMFFile.write_checkpoint",
        "time": 0.0,
        "fields": [
            {
                "symbol": "u_h",
                "checkpoint_name": "temperature",
                "meaning": "steady variable-diffusivity temperature",
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

    plt.figure(figsize=(6, 5))
    image = plot(u_h)
    plt.colorbar(image, label="temperature")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.savefig(os.path.join(artifact_dir, "heat_steady_vark.png"), dpi=180)
    plt.close()

    dx_mesh = Measure("dx", domain=mesh)
    diagnostics = {
        "source_dataset": "orange67/dataset_fenics_experiment_v7",
        "source_entry_one_based": 257,
        "mesh_subdivisions": [50, 50],
        "cells": mesh.num_cells(),
        "vertices": mesh.num_vertices(),
        "diffusivity_range": [1.0, 626.0],
        "minimum_temperature": u_h.vector().min(),
        "maximum_temperature": u_h.vector().max(),
        "domain_integral_temperature": assemble(u_h * dx_mesh),
        "elapsed_seconds": time.time() - start_time,
    }
    with open(
        os.path.join(artifact_dir, "solver_diagnostics.json"),
        "w",
        encoding="utf-8",
    ) as diagnostics_file:
        json.dump(diagnostics, diagnostics_file, indent=2)
        diagnostics_file.write("\n")
    print("Problem 63 variable-diffusivity solve complete")
