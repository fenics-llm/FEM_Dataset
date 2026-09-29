# Mesh

Unstructured triangular mesh of $\Omega(0)$ with characteristic size $h_0 \approx 1.0 \times 10^{-3}$ m.

## Benchmark Mesh Statement

Unstructured triangular mesh of $\Omega(0)$ with characteristic size $h_0 \approx 1.0 \times 10^{-3}$ m.

The mesh motion is given by $w(x,t) = s x / \|x\|$ for $x \neq (0,0)$, $w(x,t) = 0$ for $x = 0$.

## FEniCS Mesh Code

```python
from dolfin import Point
from mshr import Circle, generate_mesh

mesh = generate_mesh(Circle(Point(0.0, 0.0), 0.05), 50)
```

The solver advances the prescribed mesh motion at every time step with
`ALE.move(mesh, displacement)`, where `displacement = s*dt*x/|x|` away from
the centre and is zero at the centre.  With the resolution parameter 50, the
initial mean element size is approximately `1e-3 m` across the disk diameter.

`mesh.xdmf` and `mesh.h5` store this initial radius-$0.05$ mesh. The solver
moves only its in-memory mesh; the final radius-$0.06$ geometry is embedded in
`solution.xdmf` and `read_checkpoint.xdmf` and does not overwrite the supplied
initial mesh.
