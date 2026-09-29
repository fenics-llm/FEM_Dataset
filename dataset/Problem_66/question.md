# Problem_66: Poisson Equation with an Interior Line of Limited Regularity

**Geometry and mesh:**

Let $\Omega=(-1,1)^2$ with a uniform $4\times4$ triangular mesh aligned with
$x=0$.

**Model:**

With $\alpha=2.5$, find the scalar field $u$ satisfying

$$
-\nabla^2u=f
\quad\text{in }\Omega,
$$

where

$$
u_0(x,y)=
\begin{cases}
\cos(\pi y/2), & x\le0,\\
\cos(\pi y/2)+x^\alpha, & x>0,
\end{cases}
$$

and $f=-\nabla^2u_0$, namely

$$
f(x,y)=
\begin{cases}
\dfrac{\pi^2}{4}\cos(\pi y/2), & x\le0,\\
\dfrac{\pi^2}{4}\cos(\pi y/2)
-\alpha(\alpha-1)x^{\alpha-2}, & x>0.
\end{cases}
$$

**Boundary conditions:**

$$
u=u_0
\quad\text{on }\partial\Omega.
$$

**Numerical method:**

Use continuous piecewise-linear Lagrange elements.

**Output:**

Save $u$ to `solution.xdmf` and `solution.h5`.
