# Mesh

No explicit mesh sentence is given in the benchmark statement.

## Benchmark Mesh Statement

No explicit mesh section is present in the benchmark statement.

## FEniCS Mesh Code

```python
# ----------------------------------------------------------------------
# 1. Geometry (rectangle with a semicircular notch)
# ----------------------------------------------------------------------
a = 0.05
domain = Rectangle(Point(0.0, 0.0), Point(1.0, 0.20)) \
         - Circle(Point(0.5, 0.20), a, 64)
mesh = generate_mesh(domain, 128)
```
