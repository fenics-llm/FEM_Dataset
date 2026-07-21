# Problem_2: Solid Mechanics Problem 2 (Easy)

Source: `main-v4.tex`, Benchmark Problems, label `sm_q2`.

Original benchmark category: Solid Q2.

Reference solution source: `Results_ALL-FEM-main/reference solutions/solid/2/solid2.py`.

**Geometry:**

Let $\Omega = (0, 1.0) \times (0, 0.20)$ m be a rectangular plate.

**Mesh:**

Use uniform structured mesh with $40 \times 8$ subdivisions across $(x, y)$.

**Model:**

Plane-stress linear elasticity for displacement $u = (u_x, u_y)$ in $\Omega$.

$\sigma$ is the Cauchy stress tensor and $n$ is the unit outward normal.

**Material:**

Young's modulus $E = 200$ GPa, Poisson's ratio $\nu = 0.30$.

**Boundary conditions:**

Left edge ($x = 0$): fixed, $u_x = 0, u_y = 0$.

Top edge ($y = 0.20$): traction boundary condition is $t = \sigma n = (0, -2000)$ N/m.

Right edge ($x = 1.0$) and bottom edge ($y = 0$): traction-free, $\sigma n = 0$.

**Output:**

Save a color map of the vertical displacement $u_y$ as `q2\_uy.png`.

Save the resulting displacement field in XDMF format.
