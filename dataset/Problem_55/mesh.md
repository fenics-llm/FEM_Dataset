# Mesh

The canonical mesh is:

    mesh = UnitSquareMesh(32, 32)

The entire exterior boundary is marked for the essential condition $u=0$.
The condition $\nabla^2u=0$ is natural in the chosen formulation.

Source: Hugging Face dataset orange67/dataset_fenics_experiment_v7, entry 186
(one-based indexing). The boundary conditions were made explicit from the
source solver and its manufactured load.
