# Mesh

The canonical mesh is:

    mesh = UnitSquareMesh(20, 20)

The material interface $y=0.5$ is aligned with mesh edges. Cells with midpoint
$y\leq0.5$ have marker 0 and cells above the interface have marker 1.

Source: Hugging Face dataset orange67/dataset_fenics_experiment_v7, entry 20
(one-based indexing). The PDE sign was written in the standard
$-\nabla\cdot(k\nabla u)=0$ form, which is equivalent to the source equation.
