# Problem_8: Solid Mechanics Problem 8 (Medium)

Source: `main-v4.tex`, Benchmark Problems, label `sm_q8`.

Original benchmark category: Solid Q8.

Reference solution source: `Results_ALL-FEM-main/reference solutions/solid/8/solid8.py`.

**Geometry:**

Let $\Omega = (0, 1.0) \times (0, 0.20)$ meter be a rectangular plate modelled in 2D plane-stress.

Use structured mesh with $50 \times 25$ subdivisions over $\Omega$.

**Model:**

Small-strain, linear elasticity for displacement $u = (u_x, u_y)$ in $\Omega$.

$\sigma$ is the Cauchy stress tensor and $n$ is the unit outward normal.

**Material (orthotropic, plane-stress):**

The material’s principal axes (1--2) are rotated by $\theta = 30$ degrees anticlockwise with respect to the geometrical $x$--$y$ axes.

Orthotropic lamina properties in the local 1--2 axes:

$E_1 = 40$ GPa (Young's modulus along the principal 'grain' direction)

$E_2 = 10$ GPa (Young's modulus across the 'grain')

$G_{12} = 5$ GPa (Shear modulus in the 1--2 plane)

$\nu_{12} = 0.25$ (Poisson's ratio, for strain in direction 2 from a load in direction 1)

Use the standard plane-stress reduced stiffness $[Q]$ and rotate it to the global $x$--$y$ frame via the angle-transformed stiffness matrix at $\theta = 30$ degrees.

**Boundary conditions:**

Bottom edge ($y = 0$): Fixed, $u_x = 0, u_y = 0$.

Top edge ($y = 0.20$): Uniform downward traction $\sigma n = (0, -10 \text{ MPa} \cdot \text{m})$.

Vertical sides ($x = 0$ and $x = 1$): Traction-free, $\sigma n = (0, 0)$.

**Output:**

Save a color map of the horizontal displacement $u_x$ as `q8\_ux.png`. Also, save a color map of the von Mises stress as `q8\_vm.png`.

Also write fields to XDMF (`q8\_solution.xdmf`).

Also, save the displacement ($u$) and stress ($\sigma$) fields to `q8\_solution.xdmf`.
