"""Reference solution for Problem 52, adapted from source dataset entry 2."""

from __future__ import print_function

import json
import os
import time

from dolfin import *


problem_dir = os.path.dirname(os.path.abspath(__file__))
start_time = time.time()
comm = MPI.comm_world

# Create mesh and define function space.
mesh = UnitSquareMesh(6, 6)
V = FunctionSpace(mesh, "Lagrange", 1)
dx = Measure("dx", domain=mesh)

# Define the boundary condition from the source solver.
u0 = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)
bc = DirichletBC(V, u0, "on_boundary")

# Define the source solver's Poisson weak form.
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v)) * dx
L = f * v * dx

A = assemble(a)
b = assemble(L)
bc.apply(A, b)

# Retain the requested conjugate-gradient method and ILU preconditioner.
u_h = Function(V, name="u")
linear_solver = KrylovSolver("cg", "ilu")
linear_solver.parameters["relative_tolerance"] = 1.0e-12
linear_solver.parameters["absolute_tolerance"] = 1.0e-14
linear_solver.parameters["maximum_iterations"] = 1000
linear_solver.parameters["error_on_nonconvergence"] = True
iterations = linear_solver.solve(A, u_h.vector(), b)
u_h.rename("u", "Poisson solution")

# Canonical mesh, visualization, and checkpoint outputs.
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
    checkpoint_metadata = {
        "format": "FEniCS XDMFFile.write_checkpoint",
        "time": 0.0,
        "fields": [
            {
                "symbol": "u_h",
                "checkpoint_name": "u",
                "meaning": "scalar Poisson solution",
                "element_family": "Lagrange",
                "element_degree": 1,
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
        json.dump(checkpoint_metadata, metadata_file, indent=2)
        metadata_file.write("\n")

    diagnostics = {
        "source_dataset": "orange67/dataset_fenics_experiment_v7",
        "source_entry_one_based": 2,
        "mesh_subdivisions": [6, 6],
        "cells": mesh.num_cells(),
        "vertices": mesh.num_vertices(),
        "linear_solver": "cg",
        "preconditioner": "ilu",
        "iterations": int(iterations),
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

    print("Problem 52 solve complete: CG iterations = %d" % iterations)
