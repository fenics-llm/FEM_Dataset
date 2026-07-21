# Mesh

No explicit mesh sentence is given in the benchmark statement.

## Benchmark Mesh Statement

No explicit mesh section is present in the benchmark statement.

## FEniCS Mesh Code

```python
# mesh resolution: finer near the obstacle
mesh = generate_mesh(domain, 320)   # 80 ≈ global cell size, adjust if needed
```
