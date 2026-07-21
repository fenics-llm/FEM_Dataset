# Problem_11: Solid Mechanics Problem 11 (Medium)

Source: `main-v4.tex`, Benchmark Problems, label `sm_q11`.

Original benchmark category: Solid Q11.

Reference solution source: `Results_ALL-FEM-main/reference solutions/solid/11/solid11.py`.

**Geometry:**

Let $\Omega = (0, 1.0) \times (0, 0.20)$ m be a rectangular strip with a circular hole of radius $a = 0.04$ m, centered at $(0.50, 0.10)$.

**Model:**

Small-strain linear elasticity for displacement $u = (u_x, u_y)$ on $\Omega$.

$\sigma$ is the Cauchy stress tensor and $n$ is the outward unit normal.

**Material (plane strain, nearly incompressible):**

Young's modulus $E = 5$ MPa, Poisson's ratio $\nu = 0.49$.

**Boundary conditions and loading:**

Left edge ($x = 0$): fixed, $u_x = 0, u_y = 0$.

Right edge ($x = 1.0$): prescribed displacement, $u_x = 0.001$ m, $u_y = 0$.

Top ($y = 0.20$) and bottom ($y = 0$), and the circular hole boundary: traction-free $\sigma n = 0$.

**Output:**

Compute the von Mises equivalent stress and save a color map as `q11\_vm.png`.

Save a color map of the horizontal displacement $u_x$ as `q11\_ux.png`.

Save the resulting displacement field in XDMF format.

**Note:**

Because $\nu \approx 0.5$, pure displacement elements over-stiffen (volumetric locking). A mixed displacement-pressure formulation is standard practice for nearly incompressible solids.
