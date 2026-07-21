# Problem_3: Solid Mechanics Problem 3 (Easy)

Source: `main-v4.tex`, Benchmark Problems, label `sm_q3`.

Original benchmark category: Solid Q3.

Reference solution source: `Results_ALL-FEM-main/reference solutions/solid/3/plate_with_hole.py`.

**Geometry:**

Let $\Omega = (0, 1.0) \times (0, 0.20)$ meter be a rectangular plate with a centered circular hole of radius $a = 0.05$ m at $(0.50, 0.10)$.

**Model:**

Plane-stress linear elasticity for displacement $u = (u_x, u_y)$ in $\Omega$.

$\sigma$ is the Cauchy stress tensor and $n$ is the outward unit normal.

**Material:**

Young's modulus $E = 200$ GPa, Poisson's ratio $\nu = 0.30$.

**Boundary conditions:**

Left edge ($x = 0$): fixed, $u_x = 0, u_y = 0$.

Right edge ($x = 1.0$): $\sigma n = (2 \text{ MPa} \cdot \text{m}, 0)$.

All other external boundaries (including the hole boundary): traction-free, $\sigma n = 0$.

**Output:**

Compute the von Mises equivalent stress field (plane-stress). Save a color map as `q3\_vm.png`.

Save the resulting displacement field in XDMF format.

Report the maximum von Mises stress at the hole boundary and the stress concentration factor ($K_t = \sigma_{\max}/2 \text{ MPa}$).
