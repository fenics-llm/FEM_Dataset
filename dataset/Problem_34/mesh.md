# Mesh

No explicit mesh sentence is given in the benchmark statement.

## Benchmark Mesh Statement

No explicit mesh section is present in the benchmark statement.

## FEniCS Mesh Code

```python
# Linear-element reference candidate selected to resolve the alpha=3000
# diffuse transition more clearly than the paper's quadratic 64-by-64 case.
mesh = UnitSquareMesh(128, 128)
```
