# Problem_30: Fluid Mechanics Problem 14 (Hard)

Source: `main-v4.tex`, Benchmark Problems, label `fm_q14`.

Original benchmark category: Fluid Q14.

Reference solution source: `None`.

**Geometry:**

Solve the turbulent flow over a cylinder in the square domain $\Omega = [-30 D, +30 D] \times [-30 D, +30 D]$. The circular cylinder of diameter $D$ is placed at $(0, 0)$.

**Model:**

Solve the incompressible, unsteady Navier-Stokes equations in $\Omega$ using a Variational Multiscale (residual-based VMS) formulation with streamline-upwind, pressure-stabilizing Petrov-Galerkin, and grad-div stabilizations.

We use $u$ to denote the velocity and $p$ to denote the pressure.

**Boundary conditions:**

Impose a uniform inflow $u = (U, 0)$ at $x = -30 D$. Impose traction-free outflow with reference pressure $p = 0$ at $x = +30 D$.

Impose no slip and no penetration on the top and bottom boundaries at $y = \pm 30 D$, and on the cylinder surface. Use the initial condition $u = (0,0)$ and $p = 0$, with an optional small perturbation to trigger vortex shedding.

**Parameters:**

Set $U = 1.0$ m/s, kinematic viscosity $\nu = 2.56 \times 10^{-5}$ m$^2$/s, and density $\rho = 1.0$ kg/m$^3$, diameter $D = 1$ m.

**Output:**

Report the mean drag coefficient computed over $t \in [8.0 \text{ s}, 10.0 \text{ s}]$.

Save the final velocity field ($u$) and pressure field ($p$) at $t=10.0$ s to `vms\_solution.xdmf`.
