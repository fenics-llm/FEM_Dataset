# Problem_55: Biharmonic Equation with a C0 Interior-Penalty Method

**Geometry and mesh:**

Let $\Omega=(0,1)^2$ with the supplied $32\times32$ triangular mesh. All
quantities are nondimensional.

**Model:**

Find the scalar field $u$ satisfying

$$
\nabla^4u=4\pi^4\sin(\pi x)\sin(\pi y)
\quad\text{in }\Omega.
$$

**Boundary conditions:**

Impose the simply supported conditions

$$
u=0,\qquad \nabla^2u=0
\quad\text{on }\partial\Omega.
$$

**Numerical method:**

Use continuous quadratic Lagrange elements and the symmetric C0
interior-penalty formulation with penalty parameter $\alpha=8$.

**Output:**

Save the scalar field $u$ to `solution.xdmf` and `solution.h5`.
