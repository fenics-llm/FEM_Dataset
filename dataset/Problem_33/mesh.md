# Mesh

Use a structured mesh with $200 \times 200$ elements.

## Benchmark Mesh Statement

Use a structured mesh with $200 \times 200$ elements.

## FEniCS Mesh Code

```python
# -----------------------------------------------------------------------------
# 1. Mesh and function space
# -----------------------------------------------------------------------------
N = 200
mesh = UnitSquareMesh(N, N)
```
