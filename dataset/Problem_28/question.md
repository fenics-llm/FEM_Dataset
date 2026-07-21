# Problem_28: Fluid Mechanics Problem 12 (Hard)

Source: `main-v4.tex`, Benchmark Problems, label `fm_q12`.

Original benchmark category: Fluid Q12.

Reference solution source: `Results_ALL-FEM-main/reference solutions/fluid/12/fluid12.py`.

**Geometry:**

Let $\Omega = [0, 2.0] \times [0, 0.20]$ m be a 2-D rectangular channel of length $L = 2.0$ m and height $H = 0.20$ m.

**Model:**

Steady incompressible Navier-Stokes for velocity $u = (u_x, u_y)$ and pressure $p$ in $\Omega$, coupled to a steady advection-diffusion equation for temperature $T$:

Momentum: $\rho (u \cdot \nabla) u - \nabla \cdot [2 \mu(T) \varepsilon(u)] + \nabla p = 0$, where $\varepsilon(u) = (\nabla u + (\nabla u)^T)/2$.

Mass conservation: $\nabla \cdot u = 0$.

Energy: $u \cdot \nabla T - \kappa \nabla^2 T = 0$.

Temperature-dependent viscosity: $\mu(T) = \mu_{\text{ref}}  e^{-\beta (T - T_{\text{ref}})}$.

**Boundary conditions:**

Inlet ($x = 0$): $u_x(y) = 6 \bar{U} y (H - y) / H^2$, $u_y = 0$; $T = T_{\text{ref}}$.

Walls ($y = 0$ and $y = H$): no-slip and no penetration for the flow.

The boundary conditions at the walls for the temperature equation are:

Bottom wall ($y = 0$): Dirichlet $T = T_{\text{ref}} + 10$ K.

Top wall ($y = H$): $\partial T/\partial n = 0$.

Outlet ($x = L = 2.0$): traction-free for the flow.

For temperature: $-\kappa \partial T/\partial n = 0$.

where $n$ is outward normal vector.

**Parameters:**

Density: $\rho = 1.0$ kg m$^{-3}$.

Mean inlet speed: $\bar{U} = 1.0$ m s$^{-1}$.

Reference viscosity: $\mu_{\text{ref}} = 0.02$ Pa$\cdot$s.

$\beta = 0.05$ K$^{-1}$.

Reference temperature: $T_{\text{ref}} = 300$ K.

Thermal diffusivity: $\kappa = 1.0\times 10^{-3}$ m$^2$ s$^{-1}$.

**Output:**

Save $\mu(x, y)$ as a color map image `q13\_mu.png`.

Extract the streamwise velocity profile $u_x(y)$ along the mid-length line ($x = 1.0$) and save it to a CSV file named `q13\_profile.csv` with columns: y, ux.

Export solution fields ($u, p, T, \mu$) in XDMF format for post-processing as `q13\_solution.xdmf`.
