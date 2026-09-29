# Problem_47: Generalized-alpha elastodynamics of a cantilever

**Geometry and mesh:**

Let $\Omega=(0,1)\times(0,0.1)\times(0,0.04)$ and use the accompanying
`mesh.xdmf` tetrahedral mesh.

**Model:**

Find the time-dependent displacement
$u(\boldsymbol{x},t)=[u_x,u_y,u_z]^T$ satisfying
$$
\rho\frac{\partial^2u}{\partial t^2}
-\nabla\cdot\sigma(u)=0.
$$
The isotropic linear-elastic stress is
$$
\sigma(u)=\lambda\operatorname{tr}(\varepsilon(u))I+2\mu\varepsilon(u),
\qquad
\lambda=\frac{E\nu}{(1+\nu)(1-2\nu)},
\qquad
\mu=\frac{E}{2(1+\nu)},
$$
where $\varepsilon(u)$ is the infinitesimal (linearized) strain tensor and
$I$ is the identity tensor. Use $E=1000$, $\nu=0.3$, and density $\rho=1$.

**Initial conditions:**

Set $u=0$ and $\partial u/\partial t=0$ at $t=0$.

**Mechanical boundary conditions:**

Let $n$ be the outward unit normal. Clamp $u=0$ on $x=0$. On $x=1$,
impose $\sigma(u)n=(0,p(t),0)$, where $p(t)=t/0.8$ for
$0\le t\le0.8$ and $p(t)=0$ for $t>0.8$. Impose $\sigma(u)n=0$ on the
four lateral faces $y=0$, $y=0.1$, $z=0$, and $z=0.04$.

**Numerical method:**

Use vector-valued CG1 finite elements for the spatial discretization. Advance
the solution with a generalized-alpha time-stepping scheme using
$\alpha_m=0.2$, $\alpha_f=0.4$, $\gamma=0.7$, and $\beta=0.36$. Integrate
from $t=0$ to $T=4$ using 50 uniform time steps of size $\Delta t=0.08$.

**Output:**

Save the final displacement $u(T)$ to `solution.xdmf` and `solution.h5`.
