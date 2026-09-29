# Problem_61: Steady Heat Conduction with Piecewise Boundary Temperature

**Geometry and mesh:**

Let $\Omega=(0,5)\times(0,1)$ with a uniform $50\times10$ triangular
subdivision.

**Model:**

Find the steady temperature $u$ satisfying

$$
\nabla\cdot(k\nabla u)=0
\quad\text{in }\Omega,
\qquad k=1,
$$

**Boundary conditions:**

$$
u=10
\quad\text{on }y=1,
$$

and

$$
u=100
\quad\text{on }y=0,\ x=0,\text{ and }x=5,
$$

away from the two top corners. Assign the top value $u=10$ at those corner
nodes.

**Numerical method:**

Use continuous piecewise-linear Lagrange elements.

**Output:**

Save $u$ to `solution.xdmf` and `solution.h5`.
