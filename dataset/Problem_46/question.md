# Problem_46: Three-dimensional cantilever modal analysis

**Geometry and mesh:**

Let $\Omega=(0,20)\times(0,0.5)\times(0,1)$ and use the
accompanying `mesh.xdmf` tetrahedral mesh.

**Model:**

Find elastic eigenpairs satisfying
$$-\nabla\cdot\sigma(u)=\rho\omega^2u$$
for isotropic elasticity with $E=10^5$, $\nu=0$, and $\rho=10^{-3}$.

**Boundary conditions:**

Set $u=0$ on $x=0$ and impose zero traction elsewhere.

**Numerical method:**

Use vector CG1 elements, a consistent mass matrix, and compute the six lowest
positive modes.

**Output:**

Save displacement eigenmodes $u_1,\ldots,u_6$ to `solution.xdmf` and
`solution.h5`. Eigenvector amplitude and sign are not prescribed; comparisons
must be invariant to multiplication of each mode by a nonzero scalar.
