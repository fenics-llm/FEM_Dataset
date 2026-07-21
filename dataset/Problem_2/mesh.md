# Mesh

Use uniform structured mesh with $40 \times 8$ subdivisions across $(x, y)$.

## Benchmark Mesh Statement

Use uniform structured mesh with $40 \times 8$ subdivisions across $(x, y)$.

## FEniCS Mesh Code

```python
# ----------------------------------------------------------------------
# Mesh
# ----------------------------------------------------------------------
Lx, Ly = 1.0, 0.20
nx, ny = 40, 8
mesh = RectangleMesh(Point(0.0, 0.0), Point(Lx, Ly), nx, ny)
```
