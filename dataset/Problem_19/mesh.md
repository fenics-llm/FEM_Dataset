# Mesh

Use a uniform structured mesh of $96 \times 96$ elements.

## Benchmark Mesh Statement

Use a uniform structured mesh of $96 \times 96$ elements.

## FEniCS Mesh Code

```python
# ------------------------------------------------------------
# 1. Mesh and geometry: unit square with 96 x 96 structured cells
# ------------------------------------------------------------
nx, ny = 96, 96
mesh = UnitSquareMesh(nx, ny)
```
