# Mesh

Use a uniform structured mesh of $120 \times 12$ elements.

## Benchmark Mesh Statement

Use a uniform structured mesh of $120 \times 12$ elements.

## FEniCS Mesh Code

```python
# ----- Geometry & mesh -----
L, H = 2.0, 0.20
nx, ny = 120, 12
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), nx, ny)
```
