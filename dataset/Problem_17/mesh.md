# Mesh

Use a uniform structured mesh of $100 \times 10$ elements.

## Benchmark Mesh Statement

Use a uniform structured mesh of $100 \times 10$ elements.

## FEniCS Mesh Code

```python
# -------------------------------------------------
# 1. Geometry and mesh
# -------------------------------------------------
L, H = 2.0, 0.20
nx, ny = 100, 10
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), nx, ny, "crossed")
```
