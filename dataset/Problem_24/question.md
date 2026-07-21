# Problem_24: Fluid Mechanics Problem 8 (Easy)

Source: `main-v4.tex`, Benchmark Problems, label `fm_q8`.

Original benchmark category: Fluid Q8.

Reference solution source: `Results_ALL-FEM-main/reference solutions/fluid/8/fluid8.py`.

**Geometry:**

Let $\Omega = [0, 1] \times [0, 0.20]$ m be a rectangular channel.
Boundary partition: $\Gamma_{y0} = (0, 1)\times\{0\}$, $\Gamma_{yH} = (0, 1)\times\{0.20\}$.

**Mesh:**

Use a uniform structured mesh of $128 \times 32$ elements.

**Model:**

Steady incompressible Navier-Stokes with body force for velocity $u = (u_x, u_y)$ and pressure $p$ in $\Omega$.

Drive with a uniform body force $f = (G, 0)$.

**Boundary conditions:**

Periodic in $x$: $u$ and $p$ are periodic in $x$ direction.

Walls ($\Gamma_{y0}$ and $\Gamma_{yH}$): no-slip and no penetration.

Pressure normalization: the pressure is determined only up to an additive constant. After solving, subtract the domain-average pressure so that $\int_\Omega p \, dx = 0$.

**Parameters:**

Density $\rho = 1$ kg m$^{-3}$. Dynamic viscosity $\mu = 0.01$ Pa s. $G = 1$ N m$^{-3}$.

**Output:**

Save the velocity field ($u$) and pressure field ($p$) to `q9\_soln.xdmf`.
