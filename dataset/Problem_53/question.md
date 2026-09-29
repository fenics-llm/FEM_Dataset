# Problem_53: Diffusion Through Two Material Layers

**Geometry and mesh:**

Let

$$
\Omega=(0,1)\times(0,1)
$$

be a nondimensional square composed of

$$
\Omega_0=\{(x,y)\in\Omega:y\leq 0.5\},
\qquad
\Omega_1=\{(x,y)\in\Omega:y>0.5\}.
$$

Use the supplied $20\times20$ triangular mesh.

**Model:**

Find the scalar field $u$ satisfying

$$
-\nabla\cdot\bigl(k\nabla u\bigr)=0
\quad\text{in }\Omega,
$$

where $k=1.5$ in $\Omega_0$ and $k=50$ in $\Omega_1$. All quantities are
nondimensional.

**Boundary and interface conditions:**

$$
u(x,0)=0,\qquad u(x,1)=1,
$$

and

$$
k\nabla u\cdot n=0
\quad\text{on }x=0\text{ and }x=1.
$$

Across $y=0.5$, $u$ and the normal flux $k\nabla u\cdot n$ are continuous.

**Numerical method:**

Use continuous piecewise-linear finite elements.

**Output:**

Save the scalar field $u$ to `solution.xdmf` and `solution.h5`.
