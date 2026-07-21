# Problem_34: Multiphysics Problem 3 (Hard)

Source: `main-v4.tex`, Benchmark Problems, label `mf_q3`.

Original benchmark category: Multiphysics Q3.

Reference solution source: `None`.

**Geometry:**

Solve the spinodal decomposition of a chemical in the unit square $\Omega = (0, 1) \times (0, 1)$.

**Model:**

Solve the non-dimensional Cahn-Hilliard equations in mixed form for the concentration $c$ and chemical potential $\mu$:

$\partial c/\partial t = \nabla \cdot ( M(c) \nabla \mu )$, with $\mu = 3\alpha \mu_c - \nabla^2 c$ and $\mu_c = (0.5/ \theta) \ln(c/(1 - c)) + 1 - 2c$.

Use the degenerate mobility $M(c) = c(1 - c)$.

Initialize with $c(x, y, 0) = \bar{c} + r(x, y)$, where $\bar{c} = 0.63$ and $r$ is a zero-mean uniform perturbation in $[-0.05, 0.05]$.

Advance in time using a backward Euler scheme.

**Boundary conditions:**

Impose periodic boundary conditions for both $c$ and $\mu$ on $\partial\Omega$.

**Parameters:**

Set $\theta = 1.5$, $\alpha = 3000$, and final time $T = 0.04$.

Report the fields at $t = 0, 3\text{e}-6, 1\text{e}-4, 1\text{e}-3$, and $4\text{e}-2$ in XDMF format.

**Output:**

Save the concentration and chemical potential fields at $t = 0, 3\text{e}-6, 1\text{e}-4, 1\text{e}-3$, and $4\text{e}-2$ to a time-series file named `cahn\_hilliard.xdmf`

**Hints:**

Use adaptive time stepping: start with $\Delta t$ in the range $1\text{e}-7$ to $5\text{e}-7$ and then increase or reduce the time step based on the nonlinear iterations required for convergence. Make sure the simulation can recover if the time step selected is too large.

Discretize by splitting the fourth-order equation into two coupled second-order equations and solve them with linear finite elements.
