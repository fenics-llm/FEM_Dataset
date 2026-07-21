# Problem_27: Fluid Mechanics Problem 11 (Medium)

Source: `main-v4.tex`, Benchmark Problems, label `fm_q11`.

Original benchmark category: Fluid Q11.

Reference solution source: `Results_ALL-FEM-main/reference solutions/fluid/11/fluid11.py`.

**Geometry:**

Let $\Omega = [0, 2.0] \times [0, 0.20]$ m be a 2D rectangular channel (Length $L = 2.0$ m, Height $H = 0.20$ m).

**Mesh:**

Structured mesh with $240 \times 24$ subdivisions.

**Model:**

Steady, incompressible, non-Newtonian flow.

The governing equations for velocity $u = (u_x, u_y)$, and pressure $p$ are:

$\rho (u \cdot \nabla)u = -\nabla p + \nabla \cdot \tau$

$\nabla \cdot u = 0$

The stress tensor $\tau$ is defined by a power-law model: $\tau = 2 \mu_{\text{eff}} D$, where $D = (\nabla u + (\nabla u)^T)/2$ is the strain-rate tensor and $\mu_{\text{eff}}$ is the effective viscosity.

**Material (Power-Law Fluid):**

The effective viscosity $\mu_{\text{eff}}$ is a function of the shear rate $|D| = (2 D:D)^{1/2} + 10^{-8}$.

$\mu_{\text{eff}}(|D|) = \mu_0 [|D|^{n-1}]$.

Density $\rho = 1.0$ kg m$^{-3}$.

Consistency index $\mu_0 = 0.5$ Pa$\cdot$s$^n$.

Flow behavior index $n = 0.5$.

**Boundary conditions:**

Inlet ($x = 0$): $u_x(y) = 6 \bar{U} y (H - y) / H^2$, with $\bar{U} = 1.0$ m s$^{-1}$. $u_y = 0$.

Walls ($y = 0$ and $y = H$): No-slip and no penetration.

Outlet ($x = 2.0$): traction-free, $(-pI + \tau) n = 0$.

**Output:**

Save a color map of the velocity magnitude $|u|$ as `q12\_speed.png`.

Extract the streamwise velocity profile $u_x(y)$ at mid-channel ($x = 1.0$) and save it as `q12\_profile.csv`.

Also, save the velocity field ($u$), pressure field ($p$), and effective viscosity ($\mu_{\text{eff}}$) to `q12\_solution.xdmf`.

Report the maximum velocity $u_x(y)$ at $x = L/2$.
