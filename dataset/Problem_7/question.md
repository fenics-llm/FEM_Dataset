# Problem_7: Solid Mechanics Problem 7 (Easy)

Source: `main-v4.tex`, Benchmark Problems, label `sm_q7`.

Original benchmark category: Solid Q7.

Reference solution source: `Results_ALL-FEM-main/reference solutions/solid/7/solid7.py`.

**Geometry:**

Let $\Omega = (0, 1.0) \times (0, 0.20)$ meter be a rectangular plate in 2D plane-stress.

Define two material subdomains:

Top half $\Omega_{\text{Al}} = (0, 1.0) \times (0.10, 0.20)$ m: aluminum.

Bottom half $\Omega_{\text{Steel}} = (0, 1.0) \times (0.00, 0.10)$ m: steel.

**Mesh:**

Use structured mesh $80 \times 16$ over $\Omega$.

**Model:**

Small-strain, linear elasticity for displacement $u = (u_x, u_y)$ in $\Omega$.

$\sigma$ the Cauchy stress tensor and $n$ the unit outward normal.

There is a perfect bond at the material interface ($y = 0.10$ m). Enforce the continuity of displacement and traction across the interface.

**Material:**

$\Omega_{\text{Al}}$: Young's modulus $E = 70$ GPa, Poisson's ratio $\nu = 0.30$.

$\Omega_{\text{Steel}}$: Young's modulus $E = 200$ GPa, Poisson's ratio $\nu = 0.30$.

**Boundary conditions and loads:**

Left edge ($x = 0$): fixed, $u_x = 0, u_y = 0$.

Right edge ($x = 1$): uniform downward line traction $\sigma n = (0, -5000)$ N m$^{-1}$.

Top ($y = 0.20$) and bottom ($y = 0$): traction-free, $\sigma n = 0$.

**Output:**

Save a color map of displacement magnitude $|u|$ as `q7\_disp.png`.

Save the resulting displacement field in XDMF format.
