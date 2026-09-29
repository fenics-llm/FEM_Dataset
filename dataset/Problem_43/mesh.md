# Mesh

The solver reproduces the source notebook's meridional mesh with legacy `mshr`. It subtracts circles of radii 9 and 11 and clips the result to the first quadrant, using `generate_mesh(domain, 40)`. The resulting unstructured triangular mesh is written to `mesh.xdmf` by `solution.py`.
