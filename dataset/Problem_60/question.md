# Problem_60: Transient Heat Conduction on a Perforated Square

**Geometry and mesh:**

Let

$$
\Omega=(-1,1)^2\setminus
\left\{(x,y):(x-0.5)^2+(y-0.5)^2<0.25^2\right\}.
$$

Use the supplied conforming triangular CSG mesh generated with resolution 40.

**Model:**

Find the temperature $u$ satisfying

$$
\frac{\partial u}{\partial t}-\nabla^2u=0
\quad\text{in }\Omega,\qquad 0<t\le5.
$$

**Initial condition:**

$$
u(\mathbf{x},0)=40.
$$

**Boundary conditions:**

$$
u=10
\quad\text{on the outer square boundary},
\qquad
u=100
\quad\text{on the circular-hole boundary}.
$$

**Numerical method:**

Use continuous piecewise-linear elements and backward Euler with
$\Delta t=0.25$ for exactly 20 steps.

**Output:**

Save the temperature $u$ at $t=5$ to `solution.xdmf` and `solution.h5`.
