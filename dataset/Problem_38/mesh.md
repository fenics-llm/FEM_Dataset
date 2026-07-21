# Mesh

No explicit mesh sentence is given in the benchmark statement.

## Benchmark Mesh Statement

No explicit mesh section is present in the benchmark statement.

## FEniCS Mesh Code

```python
nx, ny = 200, 40                  # aligned with x = 0.4 and x = 0.6
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), nx, ny, "crossed")
```
