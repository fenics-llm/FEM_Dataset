# Mesh

No explicit mesh sentence is given in the benchmark statement.

## Benchmark Mesh Statement

No explicit mesh section is present in the benchmark statement.

## FEniCS Mesh Code

```python
# -------------------------
# Geometry and mesh (mshr CSG)
# -------------------------
domain = Rectangle(df.Point(0.0, 0.0), df.Point(Lx, Ly)) - Circle(df.Point(cx, cy), a, 64)
mesh = generate_mesh(domain, mesh_res)
```
