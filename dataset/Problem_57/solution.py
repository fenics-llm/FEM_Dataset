"""Reference solution for Problem 57, adapted from source dataset entry 196."""

from __future__ import print_function

import json
import os
import time

from dolfin import *


problem_dir = os.path.dirname(os.path.abspath(__file__))
start_time = time.time()
comm = MPI.comm_world

mesh = UnitSquareMesh(32, 32)
V = FunctionSpace(mesh, "Lagrange", 1)
dx = Measure("dx", domain=mesh)
ds = Measure("ds", domain=mesh)


def dirichlet_boundary(x, on_boundary):
    return on_boundary and (
        near(x[0], 0.0, DOLFIN_EPS) or near(x[0], 1.0, DOLFIN_EPS)
    )


bc = DirichletBC(V, Constant(0.0), dirichlet_boundary)
u = TrialFunction(V)
v = TestFunction(V)
f = Expression(
    "10*exp(-(pow(x[0]-0.5, 2) + pow(x[1]-0.5, 2))/0.02)",
    degree=6,
)
g = Expression("sin(5*x[0])", degree=6)
a = inner(grad(u), grad(v)) * dx
L = f * v * dx + g * v * ds

A = assemble(a)
b = assemble(L)
bc.apply(A, b)
u_h = Function(V, name="u")
solve(A, u_h.vector(), b, "mumps")
u_h.rename("u", "mixed-boundary Poisson solution")

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
                "meaning": "scalar mixed-boundary Poisson solution",
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
        "source_entry_one_based": 196,
        "mesh_subdivisions": [32, 32],
        "cells": mesh.num_cells(),
        "vertices": mesh.num_vertices(),
        "elapsed_seconds": time.time() - start_time,
    }
    artifact_dir = os.path.join(problem_dir, "reference solution artifacts")
    with open(
        os.path.join(artifact_dir, "solver_diagnostics.json"),
        "w",
        encoding="utf-8",
    ) as diagnostics_file:
        json.dump(diagnostics, diagnostics_file, indent=2)
        diagnostics_file.write("\n")
    print("Problem 57 mixed-boundary Poisson solve complete")
