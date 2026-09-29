"""Reference solution for Problem 62, adapted from source dataset entry 252."""

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

mesh = RectangleMesh(Point(-1.0, -1.0), Point(1.0, 1.0), 100, 100)
V = FunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)
a = inner(grad(u), grad(v)) * dx
L = Constant(2.0) * v * dx
bc = DirichletBC(V, Constant(0.0), DomainBoundary())

A, rhs = assemble_system(a, L, bc)
point_load = PointSource(V, Point(0.0, 0.0), 100.0)
point_load.apply(rhs)

u_h = Function(V, name="temperature")
solve(A, u_h.vector(), rhs, "mumps")
u_h.rename("temperature", "point-source temperature")

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
                "meaning": "temperature from uniform and point sources",
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
    plt.savefig(os.path.join(artifact_dir, "heat_steady_point.png"), dpi=180)
    plt.close()

    dx_mesh = Measure("dx", domain=mesh)
    diagnostics = {
        "source_dataset": "orange67/dataset_fenics_experiment_v7",
        "source_entry_one_based": 252,
        "mesh_subdivisions": [100, 100],
        "cells": mesh.num_cells(),
        "vertices": mesh.num_vertices(),
        "source_strength": 100.0,
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
    print("Problem 62 point-source heat solve complete")
