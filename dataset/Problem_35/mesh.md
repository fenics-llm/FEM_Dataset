# Mesh

No explicit mesh sentence is given in the benchmark statement. The solution
script uses a structured rectangular mesh on `(0, pi) x (-1, 1)` aligned with
the Stokes-Darcy interface at `y = 0`.

## Benchmark Mesh Statement

No explicit mesh section is present in the benchmark statement.

## FEniCS Mesh Code

```python
Lx = math.pi
y_interface = 0.0

# Structured grid aligned with the Stokes-Darcy interface y = 0.
nx, ny = 80, 80
mesh = RectangleMesh(Point(0.0, -1.0), Point(Lx, 1.0), nx, ny, "crossed")
```
