# Problem_62: Poisson Equation with a Point Heat Source

**Geometry and mesh:**

Let $\Omega=(-1,1)^2$ with a uniform $100\times100$ triangular subdivision
containing the origin as a mesh vertex.

**Model:**

Find the temperature $u$ satisfying, in the distributional sense,

$$
-\nabla^2u=2+100\,\delta_{(0,0)}
\quad\text{in }\Omega.
$$

**Boundary conditions:**

$$
u=0
\quad\text{on }\partial\Omega.
$$

**Numerical method:**

Use continuous piecewise-linear Lagrange finite elements and apply the
concentrated load with a point-source operator at the origin.

**Output:**

Save $u$ to `solution.xdmf` and `solution.h5`.
