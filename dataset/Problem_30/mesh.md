# Mesh

## Benchmark Mesh Statement

Problem 30 now follows the DFG 2D-3 benchmark geometry:

```text
[0, 2.2] x [0, 0.41] minus the closed disk of radius 0.05 centered at (0.2, 0.2).
```

The official reference uses a body-fitted quadrilateral coarse mesh followed by regular refinement. The problem does not prescribe one mandatory refinement level for the VMS discretization; the mesh must resolve the curved cylinder boundary, be refined near the cylinder and downstream wake, and support a mesh-convergence check of drag, lift, and pressure difference.

## FEniCS Mesh Code

```python
# Reference solution not available yet; add conforming DFG 2D-3 mesh
# construction with solution.py.
```
