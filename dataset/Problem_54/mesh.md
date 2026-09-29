# Mesh

The adaptive solve starts from:

    mesh = UnitSquareMesh(8, 8)

The canonical mesh.xdmf stores the final leaf mesh produced by goal-oriented
adaptive refinement. Dirichlet boundaries are $x=0$ and $x=1$; the remaining
horizontal boundaries carry the prescribed Neumann flux.

Source: Hugging Face dataset orange67/dataset_fenics_experiment_v7, entry 184
(one-based indexing). The source mathematical problem and adaptive goal were
retained.
