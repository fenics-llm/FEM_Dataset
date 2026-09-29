# Problem_56: Mean-Constrained Neumann Poisson Problem

**Geometry and mesh:**

Let $\Omega=(0,1)^2$ with the supplied $64\times64$ triangular mesh. All
quantities are nondimensional.

**Model:**

Find a scalar field $u$ and a spatially constant scalar $c$ such that

$$
-\nabla^2u+c=f
\quad\text{in }\Omega,
$$

where

$$
f(x,y)=10\exp\left(
-\frac{(x-0.5)^2+(y-0.5)^2}{0.02}
\right).
$$

**Boundary condition and normalization:**

On the entire boundary, impose

$$
\nabla u\cdot n=-\sin(5x).
$$

Fix the additive nullspace by requiring

$$
\int_\Omega u\,dx=0.
$$

**Numerical method:**

Use continuous piecewise-linear elements for $u$ and a global real element for
$c$.

**Output:**

Save $u$ and a spatially constant representation of $c$ to `solution.xdmf`
and `solution.h5`.
