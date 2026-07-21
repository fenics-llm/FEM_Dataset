# Problem_20: Fluid Mechanics Problem 4 (Easy)

Source: `main-v4.tex`, Benchmark Problems, label `fm_q4`.

Original benchmark category: Fluid Q4.

Reference solution source: `Results_ALL-FEM-main/reference solutions/fluid/4/fluid4.py`.

**Geometry:**

Let $\Omega = (0, L) \times (0, H)$ be a rectangular channel with length $L = 2.0$ m and height $H = 0.20$ m.
Boundary partition: $\Gamma_{\text{in}} = \{0\}\times[0, H]$, $\Gamma_{\text{out}} = \{L\}\times[0, H]$, $\Gamma_w = (0, L)\times\{0\} \cup (0, L)\times\{H\}$.

**Mesh:**

Use a uniform structured mesh of $160 \times 16$ elements.

**Model:**

Steady incompressible Navier-Stokes for velocity $u = (u_x, u_y)$ and pressure $p$ in $\Omega$.

**Parameters:**

Dynamic viscosity $\mu = 0.01$ Pa s and $\rho = 1$ kg m$^{-3}$.

**Boundary conditions:**

Inlet ($\Gamma_{\text{in}}$): $u_x(y) = 6 \bar{U} (y/H) (1 - y/H)$, $u_y = 0$, with mean velocity $\bar{U} = 2.5$ m s$^{-1}$.

Walls ($\Gamma_w$): no-slip and no penetration.

Outlet ($\Gamma_{\text{out}}$): traction-free: $(-p I + \mu (\nabla u + \nabla u^T)) n = 0$.

where $n$ is outward normal unit vector.

**Output:**

Save a color map of $u_x$ over $\Omega$ as `q4\_ux.png`.

Also, save the velocity field ($u$) and pressure field ($p$) to `q4\_soln.xdmf`.
