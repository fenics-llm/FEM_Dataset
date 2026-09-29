"""Reference solution for Problem 58, adapted from source dataset entry 197."""

from __future__ import print_function

import json
import os
import time

from dolfin import *


problem_dir = os.path.dirname(os.path.abspath(__file__))
start_time = time.time()
comm = MPI.comm_world

mesh = UnitCubeMesh(16, 16, 16)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
W = FunctionSpace(mesh, MixedElement([V.ufl_element(), Q.ufl_element()]))


def right(x, on_boundary):
    return on_boundary and near(x[0], 1.0)


def left(x, on_boundary):
    return on_boundary and near(x[0], 0.0)


def top_bottom(x, on_boundary):
    return on_boundary and (near(x[1], 0.0) or near(x[1], 1.0))


noslip = Constant((0.0, 0.0, 0.0))
inflow = Expression(("-sin(pi*x[1])", "0.0", "0.0"), degree=5)
zero = Constant(0.0)
bcs = [
    DirichletBC(W.sub(0), noslip, top_bottom),
    DirichletBC(W.sub(0), inflow, right),
    DirichletBC(W.sub(1), zero, left),
]

(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
f = Constant((0.0, 0.0, 0.0))

# Symmetric Stokes form for -Delta(u) + grad(p) = f and div(u) = 0.
a = (
    inner(grad(u), grad(v))
    - div(v) * p
    - q * div(u)
) * dx
L = inner(f, v) * dx

# Source-style block preconditioner.
preconditioner_form = (
    inner(grad(u), grad(v)) + inner(u, v) + p * q
) * dx
A, rhs = assemble_system(a, L, bcs)
P, _ = assemble_system(preconditioner_form, L, bcs)

mixed_solution = Function(W)
solver = KrylovSolver("minres", "amg")
solver.parameters["relative_tolerance"] = 1.0e-9
solver.parameters["absolute_tolerance"] = 1.0e-11
solver.parameters["maximum_iterations"] = 2000
solver.parameters["monitor_convergence"] = True
solver.set_operators(A, P)
iterations = solver.solve(mixed_solution.vector(), rhs)

u_mixed, p_mixed = mixed_solution.split(deepcopy=True)
u_h = Function(V, name="velocity")
p_h = Function(Q, name="pressure")
LagrangeInterpolator.interpolate(u_h, u_mixed)
LagrangeInterpolator.interpolate(p_h, p_mixed)
u_h.rename("velocity", "Stokes velocity")
p_h.rename("pressure", "Stokes pressure")

mesh_file = XDMFFile(comm, os.path.join(problem_dir, "mesh.xdmf"))
mesh_file.write(mesh)
mesh_file.close()

solution_file = XDMFFile(comm, os.path.join(problem_dir, "solution.xdmf"))
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u_h, 0.0)
solution_file.write(p_h, 0.0)
solution_file.close()

checkpoint_file = XDMFFile(
    comm, os.path.join(problem_dir, "read_checkpoint.xdmf")
)
checkpoint_file.write_checkpoint(
    u_h, "velocity", 0.0, XDMFFile.Encoding.HDF5, False
)
checkpoint_file.write_checkpoint(
    p_h, "pressure", 0.0, XDMFFile.Encoding.HDF5, True
)
checkpoint_file.close()

dx_mesh = Measure("dx", domain=mesh)
if MPI.rank(comm) == 0:
    metadata = {
        "format": "FEniCS XDMFFile.write_checkpoint",
        "time": 0.0,
        "fields": [
            {
                "symbol": "u_h",
                "checkpoint_name": "velocity",
                "meaning": "three-dimensional Stokes velocity",
                "element_family": "Lagrange",
                "element_degree": 2,
                "value_shape": [3],
            },
            {
                "symbol": "p_h",
                "checkpoint_name": "pressure",
                "meaning": "Stokes pressure",
                "element_family": "Lagrange",
                "element_degree": 1,
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
        "source_entry_one_based": 197,
        "mesh_subdivisions": [16, 16, 16],
        "cells": mesh.num_cells(),
        "vertices": mesh.num_vertices(),
        "velocity_dofs": V.dim(),
        "pressure_dofs": Q.dim(),
        "linear_iterations": iterations,
        "velocity_l2_norm": sqrt(assemble(inner(u_h, u_h) * dx_mesh)),
        "pressure_l2_norm": sqrt(assemble(p_h * p_h * dx_mesh)),
        "elapsed_seconds": time.time() - start_time,
    }
    artifact_dir = os.path.join(problem_dir, "reference solution artifacts")
    os.makedirs(artifact_dir, exist_ok=True)
    with open(
        os.path.join(artifact_dir, "solver_diagnostics.json"),
        "w",
        encoding="utf-8",
    ) as diagnostics_file:
        json.dump(diagnostics, diagnostics_file, indent=2)
        diagnostics_file.write("\n")
    print("Problem 58 three-dimensional Stokes solve complete")
