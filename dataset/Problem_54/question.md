# Problem_54: Goal-Oriented Adaptive Poisson Problem

**Geometry and initial mesh:**

Let $\Omega=(0,1)^2$ with the supplied $8\times8$ triangular mesh. All
quantities are nondimensional.

**Model:**

Find the scalar field $u$ satisfying

$$
-\nabla^2u=f
\quad\text{in }\Omega,
$$

where

$$
f(x,y)=10\exp\left(
-\frac{(x-0.5)^2+(y-0.5)^2}{0.02}
\right).
$$

**Boundary conditions:**

$$
u=0 \quad\text{on }x=0\text{ and }x=1,
$$

and

$$
\nabla u\cdot n=\sin(5x)
\quad\text{on }y=0\text{ and }y=1.
$$

**Numerical method:**

Use continuous piecewise-linear elements and goal-oriented adaptive
refinement for

$$
M(u)=\int_\Omega u\,dx.
$$

Refine until the estimated error in $M$ is below $10^{-5}$.

**Output:**

Save the final-mesh scalar field $u$ to `solution.xdmf` and `solution.h5`.
