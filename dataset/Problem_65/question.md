# Problem_65: Poisson Equation with a High-Degree Manufactured Solution

**Geometry and mesh:**

Let $\Omega=(0,1)^2$ with a uniform $10\times10$ triangular mesh.

**Model:**

Find the scalar field $u$ satisfying

$$
-\nabla^2u=f
\quad\text{in }\Omega,
$$

where, with $p=10$,

$$
u_0(x,y)
=2^{4p}x^p(1-x)^p y^p(1-y)^p,
\qquad
f=-\nabla^2u_0.
$$

**Boundary conditions:**

$$
u=u_0
\quad\text{on }\partial\Omega.
$$

**Numerical method:**

Use continuous piecewise-linear Lagrange elements and sufficiently accurate
quadrature for the high-degree manufactured load.

**Output:**

Save $u$ to `solution.xdmf` and `solution.h5`.
