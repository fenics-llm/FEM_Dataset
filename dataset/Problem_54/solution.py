"""Reference solution for Problem 54, adapted from source dataset entry 184."""

from __future__ import print_function

import json
import os
import time

from dolfin import *


problem_dir = os.path.dirname(os.path.abspath(__file__))
start_time = time.time()
comm = MPI.comm_world

# Source solver: initial mesh, P1 space, and homogeneous side conditions.
mesh = UnitSquareMesh(8, 8)
initial_cells = mesh.num_cells()
V = FunctionSpace(mesh, "Lagrange", 1)
u0 = Constant(0.0)
bc = DirichletBC(
    V,
    u0,
    "on_boundary && (near(x[0], 0.0) || near(x[0], 1.0))",
)

u = TrialFunction(V)
v = TestFunction(V)
f = Expression(
    "10*exp(-(pow(x[0]-0.5, 2) + pow(x[1]-0.5, 2))/0.02)",
    degree=4,
)
g = Expression("sin(5*x[0])", degree=4)
dx = Measure("dx", domain=mesh)
ds = Measure("ds", domain=mesh)
a = inner(grad(u), grad(v)) * dx
L = f * v * dx + g * v * ds

u_h = Function(V, name="u")
goal = u_h * dx
tolerance = 1.0e-5
problem = LinearVariationalProblem(a, L, u_h, bc)
adaptive_solver = AdaptiveLinearVariationalSolver(problem, goal)
adaptive_solver.parameters["error_control"]["dual_variational_solver"][
    "linear_solver"
] = "mumps"
adaptive_solver.solve(tolerance)

# The adaptive Function stores a hierarchy; the leaf is the requested answer.
u_final = u_h.leaf_node()
u_final.rename("u", "adaptive Poisson solution")
final_mesh = u_final.function_space().mesh()
final_dx = Measure("dx", domain=final_mesh)

mesh_file = XDMFFile(comm, os.path.join(problem_dir, "mesh.xdmf"))
mesh_file.write(final_mesh)
mesh_file.close()

solution_file = XDMFFile(comm, os.path.join(problem_dir, "solution.xdmf"))
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u_final, 0.0)
solution_file.close()

checkpoint_file = XDMFFile(
    comm, os.path.join(problem_dir, "read_checkpoint.xdmf")
)
checkpoint_file.write_checkpoint(
    u_final, "u", 0.0, XDMFFile.Encoding.HDF5, False
)
checkpoint_file.close()

if MPI.rank(comm) == 0:
    metadata = {
        "format": "FEniCS XDMFFile.write_checkpoint",
        "time": 0.0,
        "fields": [
            {
                "symbol": "u_final",
                "checkpoint_name": "u",
                "meaning": "adaptive scalar Poisson solution",
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
        json.dump(metadata, metadata_file, indent=2)
        metadata_file.write("\n")

    diagnostics = {
        "source_dataset": "orange67/dataset_fenics_experiment_v7",
        "source_entry_one_based": 184,
        "initial_cells": initial_cells,
        "final_cells": final_mesh.num_cells(),
        "final_vertices": final_mesh.num_vertices(),
        "goal_tolerance": tolerance,
        "computed_goal": assemble(u_final * final_dx),
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

    print(
        "Problem 54 adaptive solve complete: %d -> %d cells"
        % (initial_cells, final_mesh.num_cells())
    )
