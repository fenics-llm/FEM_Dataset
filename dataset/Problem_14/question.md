# Problem_14: Solid Mechanics Problem 14 (Hard)

Source: `main-v4.tex`, Benchmark Problems, label `sm_q14`.

Original benchmark category: Solid Q14.

Reference solution source: `Results_ALL-FEM-main/reference solutions/solid/14/q14_neo_hookean.py`.

**Geometry:**

Let $\Omega = (0, 1.0) \times (0, 0.20)$ m be a rectangular strip with two circular holes of radius $a = 0.04$ m, centered at $(0.40, 0.10)$ m and $(0.60, 0.10)$ m.

**Model:**

The rectangular strip is modeled as large-deformation, isotropic quasi-incompressible Neo-Hookean solid with unknown displacement $u = (u_x, u_y)$ and hydrostatic pressure $p$.

$\sigma$ is the Cauchy stress tensor and $n$ is the outward unit normal to the boundary in the deformed configuration. The specimen is under plane strain conditions.

**Material:**

Take $E = 5$ MPa and $\nu = 0.49$.

**Boundary conditions and loading:**

Left edge ($x = 0$): fixed, $u_x = 0, u_y = 0$.

Right edge ($x = 1.0$): prescribed displacement, $u_x = +0.060$ m, $u_y = 0$.

Each hole boundary: uniform follower pressure $P_{\text{hole}} = 0.10$ MPa applied as a traction $\sigma n = -P_{\text{hole}} n$.

Top ($y = 0.20$) and bottom ($y = 0$): traction-free.

**Output:**

Save a color map of the pressure field $p$ as `q14\_p.png`.

Compute and save a color map of the von Mises stress from the Cauchy stress as `q14\_vm.png`.

Save the resulting displacement field in XDMF format.
