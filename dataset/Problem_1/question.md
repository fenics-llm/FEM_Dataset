# Problem_1: Solid Mechanics Problem 1 (Easy)

Source: `main-v4.tex`, Benchmark Problems, label `sm_q1`.

Original benchmark category: Solid Q1.

Reference solution source: `Results_ALL-FEM-main/reference solutions/solid/1/ref_solid1.py`.

**Geometry:**

Let $\Omega = (0, 1.0) \times (0, 0.20)$ m be a rectangular plate.

**Mesh:**

Uniform structured mesh with $20 \times 4$ subdivisions across $(x, y)$.

**Model:**

Plane-stress linear elasticity for displacement $u = (u_x, u_y)$ in $\Omega$.

$\sigma$ is the Cauchy stress tensor and $n$ is the outward unit normal.

**Material:**

Young's modulus $E = 200$ GPa, Poisson's ratio $\nu = 0.30$.

**Boundary conditions:**

Left edge ($x = 0$): fixed, $u_x = 0, u_y = 0$.

Right edge ($x = 1$): prescribed displacement, $u_x = 0.001$ m, $u_y = 0$.

Top ($y = 0.20$) and bottom ($y = 0$): traction-free ($\sigma n = 0$).

**Output:**

Save a color map of the horizontal displacement $u_x$ as `q1\_ux.png`.

Save the resulting displacement field in XDMF format.
