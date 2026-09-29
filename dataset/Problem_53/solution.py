"""Reference solution for Problem 53, adapted from source dataset entry 20."""

from __future__ import print_function

import json
import os
import time

from dolfin import *
import numpy as np


problem_dir = os.path.dirname(os.path.abspath(__file__))
start_time = time.time()
comm = MPI.comm_world

# Create the source solver's mesh.
mesh = UnitSquareMesh(20, 20)


class Omega0(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] <= 0.5 + DOLFIN_EPS


class Omega1(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] > 0.5 + DOLFIN_EPS


subdomains = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)
for cell in cells(mesh):
    subdomains[cell] = 0 if cell.midpoint().y() <= 0.5 else 1
dx = Measure("dx", domain=mesh, subdomain_data=subdomains)

V = FunctionSpace(mesh, "CG", 1)
V0 = FunctionSpace(mesh, "DG", 0)


def bottom_boundary(x, on_boundary):
    return on_boundary and near(x[1], 0.0)


def top_boundary(x, on_boundary):
    return on_boundary and near(x[1], 1.0)


bcs = [
    DirichletBC(V, Constant(0.0), bottom_boundary),
    DirichletBC(V, Constant(1.0), top_boundary),
]

# Retain the source solver's cellwise material assignment.
k = Function(V0, name="k")
k_values = [1.5, 50.0]
k_local = k.vector().get_local()
for cell in cells(mesh):
    dof = V0.dofmap().cell_dofs(cell.index())[0]
    k_local[dof] = k_values[subdomains[cell]]
k.vector().set_local(k_local)
k.vector().apply("insert")

u = TrialFunction(V)
v = TestFunction(V)
a = inner(k * grad(u), grad(v)) * dx
L = Constant(0.0) * v * dx

A = assemble(a)
b = assemble(L)
for bc in bcs:
    bc.apply(A, b)

u_h = Function(V, name="u")
solve(A, u_h.vector(), b, "mumps")
u_h.rename("u", "layered diffusion solution")

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
                "meaning": "scalar diffusion solution",
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
        "source_entry_one_based": 20,
        "mesh_subdivisions": [20, 20],
        "cells": mesh.num_cells(),
        "vertices": mesh.num_vertices(),
        "material_coefficients": k_values,
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

    print("Problem 53 layered diffusion solve complete")
