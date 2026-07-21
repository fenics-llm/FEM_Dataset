# Problem_25: Fluid Mechanics Problem 9 (Medium)

Source: `main-v4.tex`, Benchmark Problems, label `fm_q9`.

Original benchmark category: Fluid Q9.

Reference solution source: `Results_ALL-FEM-main/reference solutions/fluid/9/fluid9.py`.

**Geometry:**

Let $\Omega = [0, 2.0] \times [0, 0.20]$ m be a rectangular channel.

**Mesh:**

Use a uniform mesh composed of $200 \times 20$ elements.

**Model:**

Fluid flow in the domain is governed by steady incompressible Navier-Stokes equations with velocity $u = (u_x, u_y)$ and pressure $p$.

After solving the fluid flow problem, solve the steady advection-diffusion for the transport of a chemical with concentration $c$, where the diffusion coefficient is $\kappa$ and the advective transport is driven by the velocity $u$.

**Boundary conditions:**

Flow ($u, p$):

Inlet ($x = 0$): $u_x(y) = 6 \bar{U} (y/H) (1 - y/H)$, with Mean inflow speed $\bar{U} = 0.1$ m s$^{-1}$ and $H = 0.20$ m, $u_y = 0$.

Walls ($y = 0$ and $y = 0.20$): no-slip and no penetration.

Outlet ($x = 2.0$): traction-free.

Concentration ($c$):

Inlet ($x = 0$): Dirichlet $c = 0$.

Walls ($y = 0$ and $y = 0.20$): impermeable, $-\kappa \nabla c \cdot n = 0$, where $n$ is unit normal vector.

Outlet ($x = 2.0$): Dirichlet $c = 1$.

**Parameters:**

Density $\rho = 1$ kg m$^{-3}$; Dynamic viscosity $\mu = 0.01$ Pa s.

Diffusivity $\kappa = 1.0 \times 10^{-3}$ m$^2$ s$^{-1}$.

**Output:**

Save the concentration field as a color map `q10\_conc.png`.

Also, save the velocity ($u$), pressure ($p$), and concentration ($c$) fields to `q10\_solution.xdmf`.
