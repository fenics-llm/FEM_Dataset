# Problem_31: Fluid Mechanics Problem 15 (Medium)

Source: `main-v4.tex`, Benchmark Problems, label `fm_q15`.

Original benchmark category: Fluid Q15.

Reference solution source: `Results_ALL-FEM-main/reference solutions/fluid/15/fluid15.py`.

**Geometry:**

Analyze the fluid flow in the unit square $[0, 1] \times [0, 1]$.

**Model:**

Solve the unsteady incompressible Navier-Stokes equations for velocity $U = (u, v)$ and pressure $p$ in the domain.

**Boundary conditions:**

Impose periodic boundary conditions on all four sides.

The initial condition for the velocity is:

$u(x, y, 0) = \sin(2\pi x) \cos(2\pi y)$

$v(x, y, 0) = -\cos(2\pi x) \sin(2\pi y)$

**Parameters:**

Set density $\rho = 1$ and kinematic viscosity $\nu = 1 \times 10^{-3}$. Simulate up to time $t = 1$.

**Numerical formulation:**

Use the supplied mesh and a periodic Taylor--Hood $[P_2]^2\times P_1$ space. Fix the pressure at $(0,0)$ to remove its constant nullspace. Advance with $\Delta t=0.0025$ using backward Euler, with convection linearized as $(U^n\cdot\nabla)U^{n+1}$.

Every velocity field after $t=0$ must be obtained by assembling and solving this finite-element time-stepping problem through all intervening steps. Direct evaluation, interpolation, or projection of an analytical or manufactured velocity field at the requested output times is not permitted.

**Numerical outputs:**

Save the velocity field at $t = 0$, $0.25$, $0.5$, and $1.0$ in `solution.xdmf` and `solution.h5`.
