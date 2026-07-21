# Problem_9: Solid Mechanics Problem 9 (Medium)

Source: `main-v4.tex`, Benchmark Problems, label `sm_q9`.

Original benchmark category: Solid Q9.

Reference solution source: `Results_ALL-FEM-main/reference solutions/solid/9/solid9.py`.

**Geometry:**

Let $\Omega = (0, 1.0) \times (0, 0.20)$ meter be a rectangular plate modelled in 2D plane-stress.

Use structured mesh with $100 \times 20$ subdivisions over $\Omega$.

**Model:**

Solve the linear elasticity equation for displacement $u = (u_x, u_y)$ in $\Omega$.

$\sigma$ is the Cauchy stress tensor and $n$ the unit outward normal.

**Material:**

Poisson's ratio $\nu = 0.30$.

Young's modulus varies with height: $E(y) = 100 \text{ GPa} + 100 \text{ GPa} \times (y / 0.20)$, for $y \in [0, 0.20]$ m.

**Boundary conditions and loads:**

Left edge ($x = 0$): fixed, $u_x = 0, u_y = 0$.

Right edge ($x = 1$): uniform traction $\sigma n = (2\times10^6 \text{ Pa} \cdot \text{m}, 0)$ per unit thickness.

Top ($y = 0.20$) and bottom ($y = 0$): traction-free, $\sigma n = 0$.

**Output:**

Save a color map of displacement magnitude $|u|$ as `q9\_disp.png`.

Save the resulting displacement field in XDMF format.
