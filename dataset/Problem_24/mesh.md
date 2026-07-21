# Mesh

Use a uniform structured mesh of $128 \times 32$ elements.

## Benchmark Mesh Statement

Use a uniform structured mesh of $128 \times 32$ elements.

## FEniCS Mesh Code

```python
# -----------------------
# Geometry and mesh
# -----------------------
Lx, Ly = 1.0, 0.20
nx, ny = 128, 32
mesh = RectangleMesh(Point(0.0, 0.0), Point(Lx, Ly), nx, ny, "crossed")  # triangles are robust
```
