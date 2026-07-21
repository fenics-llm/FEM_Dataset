# Mesh

Uniform structured mesh with $20 \times 4$ subdivisions across $(x, y)$.

## Benchmark Mesh Statement

Uniform structured mesh with $20 \times 4$ subdivisions across $(x, y)$.

## FEniCS Mesh Code

```python
# --------------------------------------------------------------
# 1. Mesh
mesh = RectangleMesh(Point(0.0, 0.0), Point(1.0, 0.20), 40, 8, "crossed")
```
