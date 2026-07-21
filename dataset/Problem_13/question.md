# Problem_13: Solid Mechanics Problem 13 (Medium)

Source: `main-v4.tex`, Benchmark Problems, label `sm_q13`.

Original benchmark category: Solid Q13.

Reference solution source: `Results_ALL-FEM-main/reference solutions/solid/13/solid13.py`.

**Geometry:**

Let $\Omega = (0, 1.0) \times (0, 0.20)$ m be a rectangular strip with a circular hole of radius $a = 0.04$ m, centered at $(0.50, 0.10)$.

**Model:**

The rectangular strip is modeled as large-deformation, isotropic incompressible Neo-Hookean solid, where the unknown fields are displacement $u = (u_x, u_y)$ and hydrostatic pressure $p$.

$\sigma$ the Cauchy stress tensor and $n$ is the outward unit normal in the deformed configuration.

**Material (plane strain):**

Take $E = 5$ MPa, $\nu = 0.5$. Use the standard incompressible neo-Hookean strain-energy with a pressure field $p$ to impose $J = 1$ (volume constraint) in a mixed ($u, p$) formulation.

**Boundary conditions and loading:**

Left edge ($x = 0$): fixed, $u_x = 0, u_y = 0$.

Hole boundary: uniform follower pressure $P_{\text{hole}} = 0.10$ MPa applied as a traction $\sigma n = -P_{\text{hole}} n$.

Right ($x = 1.0$), top ($y = 0.20$), and bottom ($y = 0$) edges: traction-free, $\sigma n = 0$.

**Output:**

Save a magnified ($\times 5$) plot of the deformed configuration as `q13\_def.png`.

Compute the von Mises equivalent stress from the Cauchy stress and save a color map as `q13\_vm.png`.

Save the resulting displacement field in XDMF format.
