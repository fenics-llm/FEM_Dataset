# Mesh

No explicit mesh sentence is given in the benchmark statement.

## Benchmark Mesh Statement

No explicit mesh section is present in the benchmark statement.

## FEniCS Mesh Code

```python
# -----------------------------------------------------------------------------
# Mesh (structured 100 x 20 over (0,1) x (0,0.20))
# -----------------------------------------------------------------------------
Lx, Ly = 1.0, 0.20
nx, ny = 100, 20
mesh = RectangleMesh(Point(0.0, 0.0), Point(Lx, Ly), nx, ny, "right")
```
