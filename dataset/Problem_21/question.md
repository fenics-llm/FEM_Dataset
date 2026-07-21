# Problem_21: Fluid Mechanics Problem 5 (Easy)

Source: `main-v4.tex`, Benchmark Problems, label `fm_q5`.

Original benchmark category: Fluid Q5.

Reference solution source: `Results_ALL-FEM-main/reference solutions/fluid/5/fluid5.py`.

**Geometry:**

Let $\Omega = (0, 1) \times (0, 1)$ be a unit square cavity.
Boundary partition: $\Gamma_{\text{left}} = \{0\}\times[0, 1]$, $\Gamma_{\text{right}} = \{1\}\times[0, 1]$, $\Gamma_{\text{bottom}} = (0, 1)\times\{0\}$, $\Gamma_{\text{top}} = (0, 1)\times\{1\}$.

**Mesh:**

Use a uniform structured mesh of $128 \times 128$ elements.

**Model:**

Lid driven cavity flow for steady incompressible Navier-Stokes for velocity $u = (u_x, u_y)$ and pressure $p$ in $\Omega$.

**Parameters:**

Given density $\rho = 1$ kg m$^{-3}$, Dynamic viscosity $\mu = 0.01$ Pa$\cdot$s.

**Boundary conditions:**

Lid ($\Gamma_{\text{top}}$): $u = (1, 0)$.

Walls ($\Gamma_{\text{left}} \cup \Gamma_{\text{right}} \cup \Gamma_{\text{bottom}}$): no-slip and no penetration.

Pressure gauge for uniqueness: fix $p = 0$ at $[0,0]$.

**Output:**

Save the color map of speed $|u|$ over $\Omega$ as `q5\_speed.png`.

Also, save the velocity field ($u$) and pressure field ($p$) to `q5\_soln.xdmf`.
