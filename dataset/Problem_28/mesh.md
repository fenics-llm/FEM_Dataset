# Mesh

No explicit mesh sentence is given in the benchmark statement.

## Benchmark Mesh Statement

No explicit mesh section is present in the benchmark statement.

## FEniCS Mesh Code

```python
# -----------------------
# Geometry and mesh
# -----------------------
L, H = 2.0, 0.20
nx, ny = 240, 24  # structured triangulation
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), nx, ny, "crossed")
```
