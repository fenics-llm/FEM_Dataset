# Problem_57: Poisson Equation with Mixed Boundary Conditions

**Geometry and mesh:**

Let $\Omega=(0,1)^2$ with the supplied $32\times32$ triangular mesh. All
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

Use continuous piecewise-linear finite elements.

**Output:**

Save $u$ to `solution.xdmf` and `solution.h5`.
