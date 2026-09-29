"""Reference solution for Problem 65, adapted from source dataset entry 281."""

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
p_value = 10
scale = float(2 ** (4 * p_value))

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, "CG", 1)
u0 = Expression(
    "scale*pow(x[0],p)*pow(1.0-x[0],p)"
    "*pow(x[1],p)*pow(1.0-x[1],p)",
    scale=scale,
    p=p_value,
    degree=10,
)
bc = DirichletBC(V, u0, DomainBoundary())

g_x = "pow(x[0],p)*pow(1.0-x[0],p)"
g_y = "pow(x[1],p)*pow(1.0-x[1],p)"
g_x_second = (
    "(p*(p-1)*pow(x[0],p-2)*pow(1.0-x[0],p)"
    "-2*p*p*pow(x[0],p-1)*pow(1.0-x[0],p-1)"
    "+p*(p-1)*pow(x[0],p)*pow(1.0-x[0],p-2))"
)
g_y_second = (
    "(p*(p-1)*pow(x[1],p-2)*pow(1.0-x[1],p)"
    "-2*p*p*pow(x[1],p-1)*pow(1.0-x[1],p-1)"
    "+p*(p-1)*pow(x[1],p)*pow(1.0-x[1],p-2))"
)
f = Expression(
    "-scale*((" + g_x_second + ")*(" + g_y + ")"
    "+(" + g_x + ")*(" + g_y_second + "))",
    scale=scale,
    p=p_value,
    degree=10,
)

u = TrialFunction(V)
v = TestFunction(V)
dx_high = Measure(
    "dx", domain=mesh, metadata={"quadrature_degree": 12}
)
a = inner(grad(u), grad(v)) * dx_high
L = f * v * dx_high

u_h = Function(V, name="u")
solve(a == L, u_h, bc, solver_parameters={"linear_solver": "mumps"})
u_h.rename("u", "high-degree manufactured Poisson solution")

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
                "meaning": "scalar manufactured Poisson solution",
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
    plot(mesh)
    plt.tight_layout()
    plt.savefig(os.path.join(artifact_dir, "mitchell1_mesh.png"), dpi=180)
    plt.close()
    plt.figure(figsize=(6, 5))
    image = plot(u_h)
    plt.colorbar(image)
    plt.tight_layout()
    plt.savefig(os.path.join(artifact_dir, "mitchell1_solution.png"), dpi=180)
    plt.close()

    diagnostics = {
        "source_dataset": "orange67/dataset_fenics_experiment_v7",
        "source_entry_one_based": 281,
        "polynomial_parameter": p_value,
        "quadrature_degree": 12,
        "mesh_subdivisions": [10, 10],
        "cells": mesh.num_cells(),
        "vertices": mesh.num_vertices(),
        "minimum": u_h.vector().min(),
        "maximum": u_h.vector().max(),
        "elapsed_seconds": time.time() - start_time,
    }
    with open(
        os.path.join(artifact_dir, "solver_diagnostics.json"),
        "w",
        encoding="utf-8",
    ) as diagnostics_file:
        json.dump(diagnostics, diagnostics_file, indent=2)
        diagnostics_file.write("\n")
    print("Problem 65 manufactured Poisson solve complete")
