# Problem_22: Fluid Mechanics Problem 6 (Medium)

Source: `main-v4.tex`, Benchmark Problems, label `fm_q6`.

Original benchmark category: Fluid Q6.

Reference solution source: `Results_ALL-FEM-main/reference solutions/fluid/6/fluid6.py`.

**Geometry:**

Let $H = 1.0$ m. The domain $\Omega$ is a 2D channel with a backward-facing step at $x = 0$. It is defined as the union of two rectangular regions:

Upstream channel: $\Omega_1 = \{ (x, y) : -3H \le x \le 0, \ 0 \le y \le H \}$

Downstream channel: $\Omega_2 = \{ (x, y) : 0 < x \le 20H, \ 0 \le y \le 2H \}$

The total domain $\Omega = \Omega_1 \cup \Omega_2$ has a 2:1 expansion ratio.

**Model:**

Steady, incompressible Navier-Stokes equations for velocity $u = (u_x, u_y)$ and pressure $p$ in $\Omega$.

**Boundary conditions:**

Inlet (at $x = -3H$, for $y \in [0, H]$): Prescribed velocity: $u_x(y) = 6 \bar{U} (y/H) (1 - y/H)$, $u_y = 0$.

Solid Walls (No-slip and no penetration): $u_x = 0, u_y = 0$ on the following boundaries:

Bottom Wall: $\{ (x, y) : y = 0, \ -3H \le x \le 20H \}$

Top Wall: $\{ (x, y) : y = H, \ -3H \le x \le 0 \} \cup \{ (x, y) : y = 2H, \ 0 \le x \le 20H \}$

Step Wall: $\{ (x, y) : x = 0, \ H \le y \le 2H \}$

Outlet (at $x = 20H$, for $y \in [0, 2H]$): Traction-free outflow: $(-p I + \mu(\nabla u + \nabla u^T)) n = 0$.

**Parameters:**

Density $\rho = 1$ kg m$^{-3}$; Dynamic viscosity $\mu = 0.01$ Pa$\cdot$s; Mean inlet speed $\bar{U} = 1.0$ m s$^{-1}$.

**Output:**

Compute the wall shear stress ($\tau_w$) on the downstream top wall for $x \in [0, 20H]$: $\tau_w(x,2H) = \mu (\partial u_x / \partial y)$.

The wall shear stress is used to find the re-attachment point, where $\tau_w(x,2H) = 0$.

Save velocity ($u$) field as `q6\_u.png`.

Also, save the velocity ($u$) and pressure ($p$) solution fields in XDMF format as `q6\_soln.xdmf`.
