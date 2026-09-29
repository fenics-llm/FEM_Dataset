"""Reference solution for Problem 55, adapted from source dataset entry 186."""

from __future__ import print_function

import json
import os
import time

from dolfin import *


problem_dir = os.path.dirname(os.path.abspath(__file__))
start_time = time.time()
comm = MPI.comm_world

# Source mesh and quadratic continuous space.
mesh = UnitSquareMesh(32, 32)
V = FunctionSpace(mesh, "CG", 2)
dx = Measure("dx", domain=mesh)
dS = Measure("dS", domain=mesh)

bc = DirichletBC(V, Constant(0.0), "on_boundary")
u = TrialFunction(V)
v = TestFunction(V)

h = CellDiameter(mesh)
h_avg = (h("+") + h("-")) / 2.0
n = FacetNormal(mesh)
f = Expression("4*pow(pi,4)*sin(pi*x[0])*sin(pi*x[1])", degree=6)
alpha = Constant(8.0)

# Retain the source solver's symmetric C0 interior-penalty form.
a = (
    inner(div(grad(u)), div(grad(v))) * dx
    - inner(avg(div(grad(u))), jump(grad(v), n)) * dS
    - inner(jump(grad(u), n), avg(div(grad(v)))) * dS
    + alpha("+") / h_avg * inner(jump(grad(u), n), jump(grad(v), n)) * dS
)
L = f * v * dx

A = assemble(a)
b = assemble(L)
bc.apply(A, b)

u_h = Function(V, name="u")
solve(A, u_h.vector(), b, "mumps")
u_h.rename("u", "biharmonic solution")

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
                "meaning": "scalar biharmonic solution",
                "element_family": "Lagrange",
                "element_degree": 2,
                "value_shape": [],
                "read_example": (
                    "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf')"
                    ".read_checkpoint(function, 'u', -1)"
                ),
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
        "source_entry_one_based": 186,
        "mesh_subdivisions": [32, 32],
        "cells": mesh.num_cells(),
        "vertices": mesh.num_vertices(),
        "element": "CG2",
        "penalty_parameter": 8.0,
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

    print("Problem 55 C0 interior-penalty solve complete")
