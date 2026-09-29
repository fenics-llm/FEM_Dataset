"""Reference solution for Problem 64, adapted from source dataset entry 270."""

from __future__ import print_function

import json
import math
import os
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dolfin import *
from mshr import Circle, generate_mesh


problem_dir = os.path.dirname(os.path.abspath(__file__))
artifact_dir = os.path.join(problem_dir, "reference solution artifacts")
os.makedirs(artifact_dir, exist_ok=True)
start_time = time.time()
comm = MPI.comm_world

T = 10.0
R = 0.3
theta = 0.2
x0 = 0.6 * R * math.cos(theta)
y0 = 0.6 * R * math.sin(theta)
sigma = 0.025

# Use the stated physical radius rather than a unit computational disk.
mesh = generate_mesh(Circle(Point(0.0, 0.0), R), 40)
V = FunctionSpace(mesh, "CG", 1)
bc = DirichletBC(V, Constant(0.0), DomainBoundary())
w = TrialFunction(V)
v = TestFunction(V)
a = Constant(T) * inner(grad(w), grad(v)) * dx
p = Expression(
    "4*exp(-0.5*(pow((x[0]-x0)/sigma,2)"
    "+pow((x[1]-y0)/sigma,2)))",
    x0=x0,
    y0=y0,
    sigma=sigma,
    degree=10,
)
L = p * v * dx

w_h = Function(V, name="deflection")
problem = LinearVariationalProblem(a, L, w_h, bc)
solver = LinearVariationalSolver(problem)
solver.parameters["linear_solver"] = "cg"
solver.parameters["preconditioner"] = "ilu"
solver.parameters["krylov_solver"]["relative_tolerance"] = 1.0e-10
solver.parameters["krylov_solver"]["absolute_tolerance"] = 1.0e-12
solver.parameters["krylov_solver"]["maximum_iterations"] = 2000
solver.parameters["krylov_solver"]["monitor_convergence"] = True
solver.solve()
w_h.rename("deflection", "circular membrane deflection")

mesh_file = XDMFFile(comm, os.path.join(problem_dir, "mesh.xdmf"))
mesh_file.write(mesh)
mesh_file.close()
solution_file = XDMFFile(comm, os.path.join(problem_dir, "solution.xdmf"))
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(w_h, 0.0)
solution_file.close()
checkpoint_file = XDMFFile(
    comm, os.path.join(problem_dir, "read_checkpoint.xdmf")
)
checkpoint_file.write_checkpoint(
    w_h, "deflection", 0.0, XDMFFile.Encoding.HDF5, False
)
checkpoint_file.close()

if MPI.rank(comm) == 0:
    metadata = {
        "format": "FEniCS XDMFFile.write_checkpoint",
        "time": 0.0,
        "fields": [
            {
                "symbol": "w_h",
                "checkpoint_name": "deflection",
                "meaning": "transverse membrane deflection",
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

    pressure_h = interpolate(p, V)
    plt.figure(figsize=(6, 5))
    plot(mesh)
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(os.path.join(artifact_dir, "membrane_mesh.png"), dpi=180)
    plt.close()
    plt.figure(figsize=(6, 5))
    image = plot(w_h)
    plt.colorbar(image)
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(os.path.join(artifact_dir, "membrane_solution.png"), dpi=180)
    plt.close()
    plt.figure(figsize=(6, 5))
    image = plot(pressure_h)
    plt.colorbar(image)
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(os.path.join(artifact_dir, "membrane_pressure.png"), dpi=180)
    plt.close()

    dx_mesh = Measure("dx", domain=mesh)
    diagnostics = {
        "source_dataset": "orange67/dataset_fenics_experiment_v7",
        "source_entry_one_based": 270,
        "mesh_resolution": 40,
        "cells": mesh.num_cells(),
        "vertices": mesh.num_vertices(),
        "radius": R,
        "tension": T,
        "source_center": [x0, y0],
        "source_sigma": sigma,
        "integrated_pressure": assemble(p * dx_mesh),
        "maximum_deflection": w_h.vector().max(),
        "elapsed_seconds": time.time() - start_time,
    }
    with open(
        os.path.join(artifact_dir, "solver_diagnostics.json"),
        "w",
        encoding="utf-8",
    ) as diagnostics_file:
        json.dump(diagnostics, diagnostics_file, indent=2)
        diagnostics_file.write("\n")
    print("Problem 64 circular membrane solve complete")
