# Problem_37: Multiphysics Problem 6 (Hard)

Source: `main-v4.tex`, Benchmark Problems, label `mf_q6`.

Original benchmark category: Multiphysics Q6.

Reference solution source: `None`.

Analyze the motion of a solid flag in a 2D flow channel using a fluid structure interaction model.

**Geometry:**

The fluid domain is a 2D channel of length $2.5$ m and height $0.41$ m. The flag is modeled as a rectangular elastic solid of length $0.35$ m and thickness $0.02$ m with the right-bottom corner placed at $(0.60 \text{ m}, 0.19 \text{ m})$ in the reference configuration.
[.2cm]
The pole of the flag is modeled as circular disk of radius $0.05$ m centered at $(0.20 \text{ m}, 0.20 \text{ m})$. The pole is assumed to be rigid and we remove this disk from our computational domain. Since this disk intersects the fluid and solid domain, portions of both the solid and fluid domain are removed.
[.2cm]
The boundary of the circular disk is denoted by $\Gamma = \Gamma_s \cup \Gamma_f$, where $\Gamma_f$ and $\Gamma_s$ are the partitions of $\Gamma$ that are in contact with the fluid and solid domain, respectively.

**Model:**

The fluid is incompressible and is solved in the moving domain using the Arbitrary Lagrangian Eulerian form of the Navier-Stokes equations. The solid flag is modeled using the St. Venant-Kirchhoff model.
[.2cm]
Enforce no-slip (fluid velocity equals solid velocity) and traction balance (fluid traction equals solid traction) at the fluid-solid interface.

**Boundary conditions:**

At the inlet ($x = 0$), prescribe a parabolic velocity profile with:

$u(t, 0, y) = \begin{cases} u_y(0,y)  \frac{1 - \cos(\pi t / 2)}{2}, & \text{if } t < 2.0 \text{ sec.}
u_y(0,y), & \text{if } t \ge 2.0 \text{ sec.} \end{cases}$

where $u_y(0,y) = 1.5 \bar{U}  y (H - y) / (H/2)^2$ with $H = 0.41$ and $\bar{U} = 1$ m/s.

At the outlet ($x = 2.5$ m), use a traction free boundary condition.

Impose no-slip and no penetration on the top and bottom channel walls.

At $\Gamma_f$: Impose no slip and no penetration for the fluid.

At $\Gamma_s$: Impose zero displacement for the solid flag.

**Parameters:**

Use fluid density $\rho_f = 1000$ kg/m$^3$ and kinematic viscosity $\nu_f = 1.0\times10^{-3}$ m$^2$/s. Use solid density $\rho_s = 10000$ kg/m$^3$, Poisson ratio $\nu_s = 0.4$, shear modulus $\mu_s = 0.5\times10^6$ Pa.

**Numerical outputs:**

Save the fluid velocity and pressure fields and the beam displacement in XDMF format.

Report the displacement components of point A with time, where the reference configuration of point A is given by $A(t=0) = (0.60 \text{ m}, 0.20 \text{ m})$.
