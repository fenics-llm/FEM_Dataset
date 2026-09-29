# Problem_63: Steady Heat Conduction with Variable Diffusivity

**Geometry and mesh:**

Let $\Omega=(0,5)^2$ with a uniform $50\times50$ triangular subdivision.

**Model:**

Find the steady temperature $u$ satisfying

$$
-\nabla\cdot\left(k(y)\nabla u\right)=0
\quad\text{in }\Omega,
$$

where

$$
k(y)=1+100y(5-y).
$$

**Boundary conditions:**

$$
u(0,y)=100y(5-y)
\quad\text{on }x=0,
$$

$$
u=0
\quad\text{on }y=0\text{ and }y=5,
$$

and the natural insulated condition

$$
k\nabla u\cdot n=0
\quad\text{on }x=5.
$$

**Numerical method:**

Use continuous piecewise-linear Lagrange elements.

**Output:**

Save $u$ to `solution.xdmf` and `solution.h5`.
