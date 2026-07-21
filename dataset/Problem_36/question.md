# Problem_36: Multiphysics Problem 5 (Hard)

Source: `main-v4.tex`, Benchmark Problems, label `mf_q5`.

Original benchmark category: Multiphysics Q5.

Reference solution source: `None`.

Analyze the fluid flow in a 2D elastic tube using a fluid structure interaction model.

**Geometry:**

The fluid domain is a rectangle of length $6$ cm and height $1$ cm. The upper and lower walls are modeled as 2-D linear-elastic solids of uniform thickness $0.1$ cm, attached along the top and bottom boundaries of the fluid domain (so the outer faces of the upper and lower walls lie at $y = 1.1$ cm and $y = -0.1$ cm, respectively).

**Model:**

The fluid is incompressible and is solved in the moving domain using the Arbitrary Lagrangian Eulerian (ALE) form of the Navier-Stokes equations. The solid walls are small-strain, 2-D (plane-strain) linear-elastic bodies.

Enforce no-slip (fluid velocity equals wall velocity) and traction balance (fluid traction matches solid traction) at the fluid-solid interface.

**Boundary conditions:**

At the inlet ($x = 0$), apply the normal-traction condition on the fluid:

$\sigma_f n_f = [ -(2\cdot10^4)/2 \cdot (1 - \cos(\pi t / 2.5\cdot10^{-3})), 0 ]^T$ for $t < 0.005$ s, and $\sigma_f n_f = 0$ thereafter.

At the outlet ($x = 6$ cm), impose zero traction, $\sigma_f n_f = 0$.

Here, $\sigma_f$ is the Cauchy stress tensor of the fluid and $n_f$ is the outward unit normal vector to the fluid domain.

The outer faces of the lower and upper walls are traction free.

The initial conditions are zero displacement and zero velocity for the solid and fluid domain respectively.

**Parameters:**

Use a fluid viscosity $\mu_f = 0.003$ poise and a fluid density $\rho_f = 1$ g$\cdot$cm$^{-3}$.

Solid density $\rho_s = 1.1$ g$\cdot$cm$^{-3}$, and Poisson ratio $\nu_s = 0.49$, Young's modulus $E_s = 3.0\times10^5$ Pa. Use a mixed displacement-pressure formulation for the discretization of the solids mechanics equations.

The simulation time step is $\Delta t = 1.0\times10^{-4}$ s.

**Numerical outputs:**

Save the output velocity and displacement at time $0.005$s and $0.1$s in XDMF format.
