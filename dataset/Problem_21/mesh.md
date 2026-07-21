# Mesh

Use a uniform structured mesh of $128 \times 128$ elements.

## Benchmark Mesh Statement

Use a uniform structured mesh of $128 \times 128$ elements.

## FEniCS Mesh Code

```python
# --- Mesh ---
nx = ny = 128
mesh = UnitSquareMesh(nx, ny)
```
