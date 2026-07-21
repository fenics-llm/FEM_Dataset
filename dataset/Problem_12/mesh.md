# Mesh

No explicit mesh sentence is given in the benchmark statement.

## Benchmark Mesh Statement

No explicit mesh section is present in the benchmark statement.

## FEniCS Mesh Code

```python
domain = Rectangle(Point(0.0, 0.0), Point(Lx, Ly)) - Circle(hole_center, a)
mesh = generate_mesh(domain, 80)   # 80 ≈ mesh resolution (adjust if needed)
```
