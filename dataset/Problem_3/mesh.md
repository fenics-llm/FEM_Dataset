# Mesh

No explicit mesh sentence is given in the benchmark statement.

## Benchmark Mesh Statement

No explicit mesh section is present in the benchmark statement.

## FEniCS Mesh Code

```python
# Create domain: rectangle minus a circle
domain = Rectangle(Point(0.0, 0.0), Point(L, H)) - Circle(center, a, 64)
mesh = generate_mesh(domain, 64)   # mesh resolution (increase for finer results)
```
