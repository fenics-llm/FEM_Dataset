# Problem_10: Solid Mechanics Problem 10 (Medium)

Source: `main-v4.tex`, Benchmark Problems, label `sm_q10`.

Original benchmark category: Solid Q10.

Reference solution source: `Results_ALL-FEM-main/reference solutions/solid/10/solid10.py`.

**Geometry:**

Let $\Omega = (0, 1.0) \times (0, 0.20)$ m be a rectangular strip analyzed in 2D.

**Mesh:**

Use structured mesh $100 \times 20$ over $\Omega$.

**Model:**

Small-strain linear elasticity for displacement $u = (u_x, u_y)$ on $\Omega$. Assume plane strain.

$\sigma$: Cauchy stress; $n$: outward unit normal.

**Material:**

Young's modulus $E = 5$ MPa, Poisson's ratio $\nu = 0.49$.

**Boundary conditions and loading:**

Left edge ($x = 0$): fixed, $u_x = 0, u_y = 0$.

Right edge ($x = 1.0$): prescribed displacement $u_x = 0.03$ m, $u_y = 0$.

Top ($y = 0.20$) and bottom ($y = 0$): traction-free, $\sigma n = 0$.

**Output:**

Save a color map of displacement magnitude $|u|$ as `q10\_disp.png`.

Save the resulting displacement field in XDMF format.

**Note:**

Because $\nu \approx 0.5$, pure displacement elements over-stiffen (volumetric locking). A mixed displacement-pressure formulation is standard practice for nearly incompressible solids.
