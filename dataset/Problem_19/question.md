# Problem_19: Fluid Mechanics Problem 3 (Easy)

Source: `main-v4.tex`, Benchmark Problems, label `fm_q3`.

Original benchmark category: Fluid Q3.

Reference solution source: `Results_ALL-FEM-main/reference solutions/fluid/3/ref_fluid3.py`.

**Geometry:**

Let $\Omega = (0, 1) \times (0, 1)$ be a unit square cavity.
Boundary partition: $\Gamma_{\text{left}} = \{0\}\times[0, 1]$, $\Gamma_{\text{right}} = \{1\}\times[0, 1]$, $\Gamma_{\text{bottom}} = (0, 1)\times\{0\}$, $\Gamma_{\text{top}} = (0, 1)\times\{1\}$.

**Mesh:**

Use a uniform structured mesh of $96 \times 96$ elements.

**Model:**

Lid driven cavity flow for steady incompressible Stokes flow for velocity $u = (u_x, u_y)$ and pressure $p$ in $\Omega$. Use an inf-sup stable finite element pair (Taylor-Hood P2-P1).

**Parameters:**

Density $\rho = 1.0$ kg m$^{-3}$; dynamic viscosity $\mu = 1.0$ Pa$\cdot$s.

**Boundary conditions:**

Lid ($\Gamma_{\text{top}}$): prescribed tangential lid motion, $u = (1, 0)$.

Other walls ($\Gamma_{\text{left}} \cup \Gamma_{\text{right}} \cup \Gamma_{\text{bottom}}$): no-slip and no penetration, $u = (0, 0)$.

**Output:**

Save a color map of speed $|u|$ over $\Omega$ as `q3\_speed.png`.

Also, save the velocity field ($u$) and pressure field ($p$) to `q3\_soln.xdmf`.
