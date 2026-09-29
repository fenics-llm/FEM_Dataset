"""Reference solution for Problem 59, adapted from source dataset entry 208."""

from __future__ import print_function

import json
import math
import os
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from dolfin import *


problem_dir = os.path.dirname(os.path.abspath(__file__))
artifact_dir = os.path.join(problem_dir, "reference solution artifacts")
os.makedirs(artifact_dir, exist_ok=True)
start_time = time.time()
comm = MPI.comm_world

element_count = 10
mesh = UnitIntervalMesh(element_count)
V = FunctionSpace(mesh, "CG", 1)
exact_solution = Expression(
    "x[0] - sinh(x[0])/sinh(1.0)",
    degree=10,
)
g_neumann = Constant(1.0 - math.cosh(1.0) / math.sinh(1.0))
bc = DirichletBC(V, exact_solution, "near(x[0], 0.0) && on_boundary")

u = TrialFunction(V)
v = TestFunction(V)
a = (inner(grad(u), grad(v)) + u * v) * dx
L = Expression("x[0]", degree=1) * v * dx + g_neumann * v * ds

u_h = Function(V, name="u")
solve(a == L, u_h, bc, solver_parameters={"linear_solver": "mumps"})
u_h.rename("u", "mixed-boundary diffusion-reaction solution")

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
    u_h, "u", 0.0, XDMFFile.Encoding.HDF5, False
)
checkpoint_file.close()

if MPI.rank(comm) == 0:
    metadata = {
        "format": "FEniCS XDMFFile.write_checkpoint",
        "time": 0.0,
        "fields": [
            {
                "symbol": "u_h",
                "checkpoint_name": "u",
                "meaning": "mixed-boundary diffusion-reaction solution",
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

    coordinates = V.tabulate_dof_coordinates().reshape(
        (-1, mesh.geometry().dim())
    )[:, 0]
    values = u_h.vector().get_local()
    order = np.argsort(coordinates)
    dense_x = np.linspace(0.0, 1.0, 401)
    exact_values = dense_x - np.sinh(dense_x) / np.sinh(1.0)
    plt.figure(figsize=(7, 5))
    plt.plot(dense_x, exact_values, label="analytical")
    plt.plot(coordinates[order], values[order], "o-", label="P1 solution")
    plt.xlabel("x")
    plt.ylabel("u")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(artifact_dir, "mixed_bvp_solution.png"), dpi=180)
    plt.close()

    diagnostics = {
        "source_dataset": "orange67/dataset_fenics_experiment_v7",
        "source_entry_one_based": 208,
        "element_count": element_count,
        "neumann_value": float(g_neumann),
        "cells": mesh.num_cells(),
        "vertices": mesh.num_vertices(),
        "elapsed_seconds": time.time() - start_time,
    }
    with open(
        os.path.join(artifact_dir, "solver_diagnostics.json"),
        "w",
        encoding="utf-8",
    ) as diagnostics_file:
        json.dump(diagnostics, diagnostics_file, indent=2)
        diagnostics_file.write("\n")
    print("Problem 59 mixed-boundary solve complete")
