# FEniCS PDE Dataset

This repository contains PDE problem statements and associated FEniCS solution
data. A subset of the problems was sourced or adapted from existing open FEM
resources, including [COMET-FEniCS](https://comet-fenics.readthedocs.io/en/latest/)
and [Example codes for coupled theories in solid
mechanics](https://github.com/SolidMechanicsCoupledTheories/example_codes).

## Dataset structure

```text
manifest.csv
dataset/
  Problem_N/
    question.md
    mesh.xdmf, mesh.h5
    solution.xdmf, solution.h5
    read_checkpoint.xdmf, read_checkpoint.h5, read_checkpoint.json
    solution.py
    ...
```

Each `Problem_N` directory contains the files associated with one PDE problem.
`question.md` states the problem. The primitive-variable solution for each
problem is stored in `solution.xdmf` with its accompanying HDF5 data. The
`read_checkpoint` files allow the solution to be read directly into a FEniCS
program; `solution.py` contains the corresponding FEniCS script.

## Manifest

`manifest.csv` indexes the dataset. Its fields are:

- **Problem number** — the problem identifier and directory name.
- **Classification** — the physical problem category.
- **Difficulty** — the assigned difficulty level.
- **Review status** — `Solution exists` or `Solution missing`.
- **One line statement of the problem** — a short description of the PDE
  problem.
