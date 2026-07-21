# Mesh

Structured mesh with $240 \times 24$ subdivisions.

## Benchmark Mesh Statement

Structured mesh with $240 \times 24$ subdivisions.

## FEniCS Mesh Code

```python
# -----------------------
# Geometry and mesh
# -----------------------
L, H = 2.0, 0.20
nx, ny = 240, 24
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), nx, ny, "crossed")  # triangles are robust
```
