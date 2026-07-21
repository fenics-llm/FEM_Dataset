# Problem_26: Fluid Mechanics Problem 10 (Medium)

Source: `main-v4.tex`, Benchmark Problems, label `fm_q10`.

Original benchmark category: Fluid Q10.

Reference solution source: `Results_ALL-FEM-main/reference solutions/fluid/10/ref_fluid10.py`.

**Geometry:**

Let $\Omega = [0, 1] \times [0, 1]$ m be a unit square cavity.

**Model (Boussinesq, steady):**

Solve the steady incompressible Navier-Stokes with body force $f$ coupled to a steady advection-diffusion equation for temperature $T$.

The notation for velocity is $u = (u_x, u_y)$, pressure is $p$, and temperature is $T$.

The body force $f$ is modeled as a Boussinesq buoyancy term, $f = [0, \rho g \beta (T - T_{\text{ref}})]$, where $\rho$ is fluid density, $g$ is gravitational acceleration, $\beta$ is volumetric thermal expansion coefficient, $T_{\text{ref}}$ is reference temperature.

**Boundary conditions:**

Left wall ($x = 0$): $T = 1$ (hot), no-slip and no penetration.

Right wall ($x = 1$): $T = 0$ (cold), no-slip and no penetration.

Top/bottom ($y = 1$ and $y = 0$): adiabatic $\partial T/\partial n = 0$, no-slip and no penetration.

where $n$ is the outward unit normal vector.

**Parameters:**

Density $\rho = 1$ kg/m$^3$.

Dynamic viscosity ($\mu$): $1.5 \times 10^{-5}$ Pa$\cdot$s.

Thermal diffusivity $\alpha = 2.1 \times 10^{-5}$ m$^2$ s$^{-1}$.

$g \beta = 3.15 \times 10^{-5}$ m s$^{-2}$ K$^{-1}$.

$T_{\text{ref}} = 0.5$ K.

**Output:**

Temperature field: save a color map as `q11\_T.png`.

Report the average Nusselt number at the left wall.

Also, save the velocity field ($u$), pressure field ($p$), and temperature field ($T$) to `q11\_solution.xdmf`.
