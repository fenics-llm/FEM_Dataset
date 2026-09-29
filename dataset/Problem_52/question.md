# Problem_52: Manufactured Poisson Equation

**Geometry and mesh:**

Let

$$
\Omega=(0,1)\times(0,1)
$$

with the supplied $6\times6$ triangular mesh. All quantities are
nondimensional.

**Model:**

Find the scalar field $u$ satisfying

$$
-\nabla^2u=-6 \quad\text{in }\Omega.
$$

**Boundary conditions:**

$$
u=1+x^2+2y^2 \quad\text{on }\partial\Omega.
$$

**Numerical method:**

Use continuous piecewise-linear finite elements. Solve the linear system with
conjugate gradients and an ILU preconditioner.

**Output:**

Save the scalar field $u$ to `solution.xdmf` and `solution.h5`.
