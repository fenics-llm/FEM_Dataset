# Mesh

No explicit mesh sentence is given in the benchmark statement.

## Benchmark Mesh Statement

No explicit mesh section is present in the benchmark statement.

## FEniCS Mesh Code

```python
upstream_channel = Rectangle(Point(-L_up, 0.0), Point(0.0, H))
downstream_channel = Rectangle(Point(0.0, 0.0), Point(L_down, 2.0*H))
mesh_resolution = 128
mesh = generate_mesh(upstream_channel + downstream_channel, mesh_resolution)
```

This is an L-shaped mesh for the actual backward-facing-step fluid domain.
Using a single full rectangle from `(-L_up, 0)` to `(L_down, 2H)` is invalid
for this problem because it includes the upstream upper solid block and makes
the upstream top wall and step wall interior lines instead of exterior no-slip
boundaries.
