# Mesh

Use a uniform mesh composed of $200 \times 20$ elements.

## Benchmark Mesh Statement

Use a uniform mesh composed of $200 \times 20$ elements.

## FEniCS Mesh Code

```python
# ------------------------------------------------------------
# 1. Mesh and geometry
# ------------------------------------------------------------
L, H = 2.0, 0.20
nx, ny = 200, 20
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), nx, ny, "crossed")
```
