# Problem_42: Self-weight cantilever in plane stress

**Geometry and mesh:**

Let $\Omega=(0,25)\times(0,1)$. Use the accompanying `mesh.xdmf`, a crossed
$250\times10$ rectangular mesh.

**Model:**

Find $u=(u_x,u_y)$ such that
$$
-\nabla\cdot\sigma(u)=f,\qquad
\sigma(u)=\lambda_*\operatorname{tr}(\varepsilon(u))I+2\mu\varepsilon(u),
$$
where $\varepsilon(u)$ is the infinitesimal (linearized) strain tensor,
$E=10^5$, $\nu=0.3$, $\mu=E/[2(1+\nu)]$,
$\lambda_*=E\nu/(1-\nu^2)$, and $f=(0,-10^{-3})$.

**Boundary conditions:**

Impose $u=0$ on $x=0$ and zero traction on $\partial\Omega\setminus\{x=0\}$.

**Numerical method:**

Use continuous piecewise-quadratic vector Lagrange elements.

**Output:**

Save the displacement $u$ to `solution.xdmf` and `solution.h5`.
