# Problem_18: Fluid Mechanics Problem 2 (Easy)

Source: `main-v4.tex`, Benchmark Problems, label `fm_q2`.

Original benchmark category: Fluid Q2.

Reference solution source: `Results_ALL-FEM-main/reference solutions/fluid/2/fluid2.py`.

**Geometry:**

Let $\Omega = (0, L) \times (0, H)$ be a rectangular channel with length $L = 2.0$ m and height $H = 0.20$ m.
Boundary partition: $\Gamma_{\text{in}} = \{0\}\times[0, H]$, $\Gamma_{\text{out}} = \{L\}\times[0, H]$, $\Gamma_w = (0, L)\times\{0\} \cup (0, L)\times\{H\}$.

**Mesh:**

Use a uniform structured mesh of $120 \times 12$ elements.

**Model:**

Steady incompressible Stokes flow with body force for velocity $u = (u_x, u_y)$ and pressure $p$ in $\Omega$.

Uniform body force: $f = (1.0, 0.0)$ N m$^{-3}$.

**Parameters:**

Dynamic viscosity $\mu = 1.0$ Pa$\cdot$s; density $\rho = 1.0$ kg m$^{-3}$.

**Boundary conditions:**

Walls ($\Gamma_w$): no-slip and no penetration.

Inlet and outlet ($\Gamma_{\text{in}} \cup \Gamma_{\text{out}}$): traction-free natural condition: $(-p I + \mu(\nabla u + \nabla u^T)) n = 0$.

where $n$ denotes the unit outward normal vector.

**Output:**

Save a color map of speed $|u|$ over $\Omega$ as `q2\_speed.png`.

Also, save the velocity field ($u$) and pressure field ($p$) to `q2\_solution.xdmf`.

**Validation note:** Compare only the velocity field (including the maximum speed), not the pressure field.
