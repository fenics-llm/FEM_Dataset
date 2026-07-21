# Mesh

Use structured mesh $80 \times 16$ over $\Omega$.

## Benchmark Mesh Statement

Use structured mesh $80 \times 16$ over $\Omega$.

## FEniCS Mesh Code

```python
# --- Geometry & mesh (structured 80 x 16 over (0,1.0) x (0,0.20) m) ---
Lx, Ly = 1.0, 0.20
mesh = RectangleMesh(Point(0.0, 0.0), Point(Lx, Ly), 80, 16, "left/right")
```
