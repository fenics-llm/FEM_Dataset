# Problem_38: Multiphysics Problem 7 (Hard)

Source: `main-v4.tex`, Benchmark Problems, label `mf_q7`.

Original benchmark category: Multiphysics Q7.

Reference solution source: `Results_ALL-FEM-main/reference solutions/multiphysics - Copy/7/multi-7.py`.

**Geometry:**

Let $\Omega = [0, 1.0] \times [0, 0.20]$ m be a 2D rectangular channel.

The domain $\Omega$ is divided into two subdomains:

A porous filter: $\Pi = [0.4, 0.6] \times [0, 0.20]$ m.

A free-fluid region: $\Omega_f = \Omega \setminus \Pi$.

The interface between them is $\partial\Pi$.

**Model:**

Solve a coupled, steady, incompressible flow problem for velocity $u = (u_x, u_y)$ and pressure $p$. The governing equations differ by subdomain:

In $\Omega_f$ (Fluid): Steady Navier-Stokes equations.

In $\Pi$ (Porous): Steady incompressible Brinkman equations given by

$-\nabla p + \mu \nabla^2 u - (\mu/K) u = 0$

$\nabla \cdot u = 0$

At $\partial\Pi$ (Interface): Continuity of velocity and traction is enforced.

**Boundary Conditions:**

Inlet ($x = 0$): $u_x(y) = 6 \bar{U} y (H - y) / H^2$, with $\bar{U} = 1.0$ m s$^{-1}$. $u_y = 0$.

Walls ($y = 0$ and $y = H$): No-slip and no penetration.

Outlet ($x = 1.0$): traction-free, $(-pI + \mu(\nabla u + \nabla u^T)) n = 0$.

**Parameters:**

Density ($\rho$): $1.0$ kg m$^{-3}$

Dynamic viscosity ($\mu$): $0.01$ Pa$\cdot$s.

Permeability ($K$): $1.0 \times 10^{-6}$ m$^2$.

**Output:**

Save a color map of the velocity magnitude $|u|$ as `q14\_speed.png`.

Compute the pressure drop ($\Delta p$) across the porous block by sampling $p$ at the centerline just before and just after $\Pi$. Save this value to `q14\_dp.txt`.

Also, save the velocity field ($u$) and the pressure field ($p$) for the entire domain to `q14\_solution.xdmf`.
