# Problem_58: Three-Dimensional Stokes Flow with Taylor-Hood Elements

**Geometry and mesh:**

Let

$$
\Omega=(0,1)^3.
$$

Use the supplied unit-cube mesh with 16 subdivisions in each coordinate
direction.

**Model:**

Find the velocity $\mathbf{u}:\Omega\to\mathbb{R}^3$ and pressure
$p:\Omega\to\mathbb{R}$ satisfying

$$
-\nabla^2\mathbf{u}+\nabla p=\mathbf{0},
\qquad
\nabla\cdot\mathbf{u}=0
\quad\text{in }\Omega.
$$

**Boundary conditions:**

Impose no slip on the two $y$-normal walls,

$$
\mathbf{u}=\mathbf{0}
\quad\text{on }y=0\text{ and }y=1,
$$

the inflow profile

$$
\mathbf{u}=(-\sin(\pi y),0,0)
\quad\text{on }x=1,
$$

and the pressure condition

$$
p=0
\quad\text{on }x=0.
$$

Impose a zero-traction boundary condition on the boundary portions where
velocity is not prescribed, namely $x=0$ and $z=0,1$.

**Numerical method:**

Use Taylor--Hood elements: continuous quadratic velocity and continuous
linear pressure.

**Output:**

Save $\mathbf{u}$ and $p$ to `solution.xdmf` and `solution.h5`.
