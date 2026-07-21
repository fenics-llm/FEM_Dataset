# Problem_23: Fluid Mechanics Problem 7 (Medium)

Source: `main-v4.tex`, Benchmark Problems, label `fm_q7`.

Original benchmark category: Fluid Q7.

Reference solution source: `Results_ALL-FEM-main/reference solutions/fluid/7/ref_fluid7.py`.

**Geometry:**

Let $\Omega = [0, 2.2] \times [0, 0.41]$ m be a rectangular channel with a circular hole of radius $R = 0.05$ m centered at $(0.20, 0.20)$ m.

**Model:**

Steady incompressible Navier-Stokes for velocity $u = (u_x, u_y)$ and pressure $p$ in $\Omega$.

**Boundary conditions:**

Inlet ($x = 0$): $u_x(y) = 6 \bar{U} y (H - y) / H^2$, $u_y = 0$, where $H = 0.41$.

Walls ($y = 0$ and $y = 0.41$) and circular boundary: no-slip and no penetration.

Outlet ($x = 2.2$): traction-free.

**Parameters:**

Dynamic viscosity $\mu = 0.001$ Pa s, Density $\rho = 1$ kg m$^{-3}$, and Mean inlet velocity $\bar{U} = 0.2$ m/s.

**Output:**

Drag coefficient $C_D$:

Compute $C_D$ from the drag force on the circle $F_D$ via $C_D = 2 F_D / (\rho \bar{U}^2 D)$.

Save a color map of speed $|u|$ over $\Omega$ as `q7\_speed.png`.

Also, save the velocity field ($u$) and pressure field ($p$) to `q7\_soln.xdmf`.
