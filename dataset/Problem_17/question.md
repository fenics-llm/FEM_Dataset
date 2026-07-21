# Problem_17: Fluid Mechanics Problem 1 (Easy)

Source: `main-v4.tex`, Benchmark Problems, label `fm_q1`.

Original benchmark category: Fluid Q1.

Reference solution source: `Results_ALL-FEM-main/reference solutions/fluid/1/q1_stokes.py`.

**Geometry:**

Let $\Omega = (0, L) \times (0, H)$ be a rectangular channel with length $L = 2.0$ m and height $H = 0.20$ m.
Boundary partition: $\Gamma_{\text{in}} = \{0\}\times[0, H]$, $\Gamma_{\text{out}} = \{L\}\times[0, H]$, $\Gamma_w = (0, L)\times\{0\} \cup (0, L)\times\{H\}$.

**Mesh:**

Use a uniform structured mesh of $100 \times 10$ elements.

**Model:**

Steady incompressible Stokes flow in $\Omega$. Let velocity be $u = (u_x, u_y)$, pressure $p$ and Cauchy stress tensor $\sigma$.

**Parameters:**

Dynamic viscosity $\mu = 1.0$ Pa$\cdot$s; density $\rho = 1.0$ kg$\cdot$m$^{-3}$.

**Boundary conditions:**

Walls ($\Gamma_w$): no-slip and no penetration.

Inlet ($\Gamma_{\text{in}}$): traction (natural) boundary, $\sigma(u,p) n = - p_{\text{in}} n$ on $\Gamma_{\text{in}}$, with $p_{\text{in}} = 1.0$ Pa.

Outlet ($\Gamma_{\text{out}}$): traction (natural) boundary, $\sigma(u,p) n = - p_{\text{out}} n$ on $\Gamma_{\text{out}}$, with $p_{\text{out}} = 0$ Pa.

$n$ denotes the unit outward normal vector.

**Output:**

Save a color map of the speed $|u|$ over $\Omega$ to `q1\_speed.png`.

Save the velocity field ($u$) and pressure field ($p$) to `q1\_soln.xdmf`.
