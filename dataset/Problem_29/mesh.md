# Mesh

Use a uniform mesh composed of $100 \times 10$ elements.

## Benchmark Mesh Statement

Use a uniform mesh composed of $100 \times 10$ elements.

## FEniCS Mesh Code

```python
# ----------------------------
# Mesh and function space
# ----------------------------
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), nx, ny)
```
