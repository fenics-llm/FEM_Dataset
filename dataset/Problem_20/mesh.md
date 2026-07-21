# Mesh

Use a uniform structured mesh of $160 \times 16$ elements.

## Benchmark Mesh Statement

Use a uniform structured mesh of $160 \times 16$ elements.

## FEniCS Mesh Code

```python
# --- Geometry & mesh ---
L = 2.0
H = 0.20
nx, ny = 160, 16
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), nx, ny)
```
