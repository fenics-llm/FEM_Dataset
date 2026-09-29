# Mesh

The canonical mesh is the default diagonal subdivision produced by

```python
mesh = UnitSquareMesh(6, 6)
```

It contains 72 triangular cells and 49 vertices. The entire exterior boundary
is Dirichlet.

Source: Hugging Face dataset `orange67/dataset_fenics_experiment_v7`, entry 2
(one-based indexing). The source problem and solver were retained without a
mathematical change; canonical output and validation were added.
