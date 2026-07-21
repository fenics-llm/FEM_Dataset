# Problem_39: Multiphysics Problem 8 (Hard)

Source: `main-v4.tex`, Benchmark Problems, label `mf_q8`.

Original benchmark category: Multiphysics Q8.

Reference solution source: `None`.

**Geometry:**

The domain $[0, 1] \times [0, H]$ m is split into two adjacent rectangular subdomains:

Free Fluid ($\Omega_f$): $[0, 0.6] \times [0, H]$ m

Porous Medium ($\Omega_p$): $[0.6, 1.0] \times [0, H]$ m

Interface ($\Gamma$): The boundary between them at $x = 0.6$.

**Model:**

This is a coupled problem solving for velocity and pressure in each subdomain.

In $\Omega_f$ (Fluid): Steady Stokes equations for velocity $u_f$ and pressure $p_f$.

In $\Omega_p$ (Porous): Steady Darcy flow for velocity $u_p$ and pressure $p_p$.

$u_p = -(K/\mu) \nabla p_p$

$\nabla \cdot u_p = 0$

**Interface Conditions (on $\Gamma$ at $x = 0.6$)**

The two models are coupled by three conditions:

Continuity of Normal Velocity: $u_f \cdot n = u_p \cdot n$, where $n$ is the normal vector.

Continuity of Pressure: $p_f = p_p$.

Impose no slip for the Stokes velocity at the interface: $u_t = 0$, where $u_t$ is the tangential velocity on the fluid side.

**Boundary Conditions (External)**

Inlet ($x = 0$, on $\Omega_f$): Prescribed velocity: $u_x(y) = 6 \bar{U} y (H - y) / H^2$, $u_y = 0$.

Fluid Walls ($y = 0$ and $y = H$, on $\Omega_f$): No-slip and no penetration.

Porous Walls ($y = 0$ and $y = H$, on $\Omega_p$): No-flux (impermeable)

Outlet ($x = 1.0$, on $\Omega_p$): Fixed pressure, $p_p = 0$ (gauge pressure).

**Parameters:**

Dynamic Viscosity ($\mu$): $0.02$ Pa$\cdot$s

Permeability ($K$): $1.0 \times 10^{-6}$ m$^2$

Mean Inlet Speed ($\bar{U}$): $0.1$ m s$^{-1}$

Channel Height ($H$): $0.2$ m

**Output:**

Save the interface profiles for normal velocity ($u_x$) and tangential velocity ($u_y$) along $\Gamma$ (from $y=0$ to $y=H$) to `q15\_interface.csv`.

Save a color map of the pressure field $p$ (for both $p_f$ and $p_p$) as `q15\_p.png`.

Also, save the combined velocity field ($u$) and pressure field ($p$) for the entire domain to `q15\_solution.xdmf`.
