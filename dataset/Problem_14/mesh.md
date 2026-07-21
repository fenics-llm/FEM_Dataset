# Mesh

No explicit mesh sentence is given in the benchmark statement.

## Benchmark Mesh Statement

No explicit mesh section is present in the benchmark statement.

## FEniCS Mesh Code

```python
# --------------------------------------------------------------
# Geometry and mesh
# --------------------------------------------------------------
L, H = 1.0, 0.20
a = 0.04
hole1 = Circle(Point(0.40, 0.10), a, 64)
hole2 = Circle(Point(0.60, 0.10), a, 64)
domain = Rectangle(Point(0.0, 0.0), Point(L, H)) - hole1 - hole2
mesh = generate_mesh(domain, 80)   # increase for finer resolution
```
