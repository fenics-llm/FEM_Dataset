# Mesh

The canonical mesh is:

    mesh = UnitSquareMesh(32, 32)

The vertical sides are Dirichlet boundaries. The horizontal sides carry the
prescribed Neumann flux.

Source: Hugging Face dataset orange67/dataset_fenics_experiment_v7, entry 196
(one-based indexing). The source's unspecified Neumann location was resolved
as the non-Dirichlet horizontal boundary, matching its weak form.
