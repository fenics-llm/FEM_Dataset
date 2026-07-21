# Mesh

Use structured mesh $100 \times 20$ over $\Omega$.

## Benchmark Mesh Statement

Use structured mesh $100 \times 20$ over $\Omega$.

## FEniCS Mesh Code

```python
# ----------------------------
# Geometry and mesh (100 x 20)
# ----------------------------
Lx, Ly = 1.0, 0.20
nx, ny = 100, 20
mesh = RectangleMesh(Point(0.0, 0.0), Point(Lx, Ly), nx, ny, "right/left")
```
