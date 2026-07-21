# Mesh

No explicit mesh sentence is given in the benchmark statement.

## Benchmark Mesh Statement

No explicit mesh section is present in the benchmark statement.

## FEniCS Mesh Code

```python
# ------------------------------------------------------------
# Geometry & mesh: Omega = (0,1.0) x (0,0.20), structured 50 x 25
# ------------------------------------------------------------
Lx, Ly = 1.0, 0.20
nx, ny = 50, 25
mesh = RectangleMesh(Point(0.0, 0.0), Point(Lx, Ly), nx, ny)
```
