"""Reference solution for Problem 56, adapted from source dataset entry 193."""

from __future__ import print_function

import json
import os
import time

from dolfin import *


problem_dir = os.path.dirname(os.path.abspath(__file__))
start_time = time.time()
comm = MPI.comm_world

mesh = UnitSquareMesh(64, 64)
P1 = FiniteElement("CG", mesh.ufl_cell(), 1)
R = FiniteElement("Real", mesh.ufl_cell(), 0)
W = FunctionSpace(mesh, MixedElement([P1, R]))
V = FunctionSpace(mesh, "CG", 1)
C = FunctionSpace(mesh, "DG", 0)
dx = Measure("dx", domain=mesh)
ds = Measure("ds", domain=mesh)

# Correct the source typo TrialFunction(W) to TrialFunctions(W).
(u, c) = TrialFunctions(W)
(v, d) = TestFunctions(W)
f = Expression(
    "10*exp(-(pow(x[0]-0.5, 2) + pow(x[1]-0.5, 2))/0.02)",
    degree=6,
)
g = Expression("-sin(5*x[0])", degree=6)

# Source mixed form: c enforces compatibility and d enforces mean(u)=0.
a = (inner(grad(u), grad(v)) + c * v + u * d) * dx
L = f * v * dx + g * v * ds

A = assemble(a)
b = assemble(L)
w_h = Function(W)
solve(A, w_h.vector(), b, "mumps")
u_mixed, c_mixed = w_h.split(deepcopy=True)

u_h = Function(V, name="u")
LagrangeInterpolator.interpolate(u_h, u_mixed)
u_h.rename("u", "mean-zero Neumann solution")
c_value = float(c_mixed.vector().get_local()[0])
c_h = interpolate(Constant(c_value), C)
c_h.rename("c", "constant compatibility multiplier")

mesh_file = XDMFFile(comm, os.path.join(problem_dir, "mesh.xdmf"))
mesh_file.write(mesh)
mesh_file.close()

solution_file = XDMFFile(comm, os.path.join(problem_dir, "solution.xdmf"))
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u_h, 0.0)
solution_file.write(c_h, 0.0)
solution_file.close()

checkpoint_file = XDMFFile(
    comm, os.path.join(problem_dir, "read_checkpoint.xdmf")
)
checkpoint_file.write_checkpoint(
    u_h, "u", 0.0, XDMFFile.Encoding.HDF5, False
)
checkpoint_file.write_checkpoint(
    c_h, "c", 0.0, XDMFFile.Encoding.HDF5, True
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
                "meaning": "mean-zero scalar field",
                "element_family": "Lagrange",
                "element_degree": 1,
                "value_shape": [],
            },
            {
                "symbol": "c_h",
                "checkpoint_name": "c",
                "meaning": "spatially constant compatibility multiplier",
                "element_family": "Discontinuous Lagrange",
                "element_degree": 0,
                "value_shape": [],
            },
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
        "source_entry_one_based": 193,
        "mesh_subdivisions": [64, 64],
        "cells": mesh.num_cells(),
        "vertices": mesh.num_vertices(),
        "compatibility_multiplier": c_value,
        "mean_u": assemble(u_h * dx),
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

    print("Problem 56 solve complete: c = %.12e" % c_value)
