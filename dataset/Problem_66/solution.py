"""Reference solution for Problem 66, adapted from source dataset entry 296."""

from __future__ import print_function

import json
import math
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
alpha = 2.5

mesh = RectangleMesh(Point(-1.0, -1.0), Point(1.0, 1.0), 4, 4)
V = FunctionSpace(mesh, "CG", 1)
u0 = Expression(
    "(x[0] <= 0.0) ? cos(pi*x[1]/2.0)"
    " : cos(pi*x[1]/2.0) + pow(x[0], alpha)",
    alpha=alpha,
    pi=math.pi,
    degree=5,
)
f = Expression(
    "(x[0] <= 0.0)"
    " ? (pi*pi/4.0)*cos(pi*x[1]/2.0)"
    " : (pi*pi/4.0)*cos(pi*x[1]/2.0)"
    " - alpha*(alpha-1.0)*pow(x[0], alpha-2.0)",
    alpha=alpha,
    pi=math.pi,
    degree=5,
)
bc = DirichletBC(V, u0, DomainBoundary())
u = TrialFunction(V)
v = TestFunction(V)
dx_piecewise = Measure(
    "dx", domain=mesh, metadata={"quadrature_degree": 8}
)
a = inner(grad(u), grad(v)) * dx_piecewise
L = f * v * dx_piecewise

u_h = Function(V, name="u")
solve(a == L, u_h, bc, solver_parameters={"linear_solver": "mumps"})
u_h.rename("u", "limited-regularity Poisson solution")

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
                "meaning": "piecewise manufactured Poisson solution",
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
    plt.savefig(os.path.join(artifact_dir, "mitchell10_mesh.png"), dpi=180)
    plt.close()
    plt.figure(figsize=(6, 5))
    image = plot(u_h)
    plt.colorbar(image)
    plt.tight_layout()
    plt.savefig(os.path.join(artifact_dir, "mitchell10_solution.png"), dpi=180)
    plt.close()

    diagnostics = {
        "source_dataset": "orange67/dataset_fenics_experiment_v7",
        "source_entry_one_based": 296,
        "alpha": alpha,
        "interface_x": 0.0,
        "quadrature_degree": 8,
        "mesh_subdivisions": [4, 4],
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
    print("Problem 66 limited-regularity Poisson solve complete")
