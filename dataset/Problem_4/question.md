# Problem_4: Solid Mechanics Problem 4 (Easy)

Source: `main-v4.tex`, Benchmark Problems, label `sm_q4`.

Original benchmark category: Solid Q4.

Reference solution source: `Results_ALL-FEM-main/reference solutions/solid/4/ref_solid4.py`.

**Geometry:**

Let $\Omega = (0, 1.0) \times (0, 0.20)$ meter be a rectangular plate with two circular holes of radius $a = 0.04$ m centered at $(0.33, 0.10)$ m and $(0.67, 0.10)$ m.

**Model:**

Plane-stress linear elasticity for displacement $u = (u_x, u_y)$ in $\Omega$.

$\sigma$ the Cauchy stress tensor and $n$ is the unit outward normal.

**Material:**

Young's modulus $E = 200$ GPa, Poisson's ratio $\nu = 0.30$.

**Boundary conditions:**

Left edge ($x = 0$): clamped, $u_x = 0, u_y = 0$.

Right edge ($x = 1.0$): $\sigma n = (2 \text{ MPa} \cdot \text{m}, 0)$.

Remaining boundaries (top, bottom, and both hole boundaries): traction-free, $\sigma n = 0$.

**Output:**

Compute the von Mises equivalent stress and save a color map as `q4\_vm.png`.

Save the resulting displacement field in XDMF format.

Report the maximum von Mises stress at the hole boundary and the stress concentration factor ($K_t = \sigma_{\max}/2 \text{ MPa}$).
