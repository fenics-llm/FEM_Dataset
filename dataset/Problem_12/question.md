# Problem_12: Solid Mechanics Problem 12 (Medium)

Source: `main-v4.tex`, Benchmark Problems, label `sm_q12`.

Original benchmark category: Solid Q12.

Reference solution source: `Results_ALL-FEM-main/reference solutions/solid/12/solid12.py`.

**Geometry:**

Let $\Omega = (0, 1.0) \times (0, 0.20)$ m be a rectangular strip with a circular hole of radius $a = 0.04$ m, centered at $(0.50, 0.10)$.

**Model:**

The rectangular strip is modeled as large-deformation, isotropic incompressible neo-Hookean solid, where the unknown fields are displacement $u = (u_x, u_y)$ and hydrostatic pressure $p$.

$\sigma$ is the Cauchy stress tensor and $n$ is the outward unit normal.

**Material (plane strain):**

Young's modulus $E = 5$ MPa and Poisson ratio $\nu = 0.5$. Use the standard incompressible neo-Hookean strain-energy with a pressure field $p$ to impose $J = 1$ (volume constraint) in a mixed ($u, p$) formulation.

**Boundary conditions and loading:**

Left edge ($x = 0$): fixed, $u_x = 0, u_y = 0$.

Right edge ($x = 1.0$): prescribed displacement, $u_x = 0.060$ m, $u_y = 0$.

Top ($y = 0.20$), bottom ($y = 0$), and the circular hole boundary: traction-free, $\sigma n = 0$.

**Output:**

Save a color map of hydrostatic pressure $p$ as `q12\_p.png`.

Compute and save a color map of von Mises stress (from the Cauchy stress) as `q12\_vm.png`.

Save the resulting displacement field in XDMF format.
