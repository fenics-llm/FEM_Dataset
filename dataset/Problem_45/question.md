# Problem_45: Fully coupled transient thermoelasticity

**Geometry and mesh:**

Let
$\Omega=(0,1\,\mathrm m)^2\setminus\{x^2+y^2<(0.1\,\mathrm m)^2\}$ in
the first quadrant and use the accompanying `mesh.xdmf` triangular mesh.

**Model:**

Let $\Theta(\boldsymbol{x},t)$ be the temperature increment from the reference
temperature $T_0$, given by $\Theta=T-T_0$. For
$t\in[10,10^4]\,\mathrm s$, find the time-dependent plane-strain displacement
$u(\boldsymbol{x},t)=[u_x,u_y]^T$ and temperature increment
$\Theta(\boldsymbol{x},t)$
satisfying the coupled quasistatic thermoelastic system
$$
-\nabla\cdot\sigma(u,\Theta)=0,\qquad
\sigma(u,\Theta)=\lambda\operatorname{tr}(\varepsilon(u))I
+2\mu\varepsilon(u)-\kappa\Theta I,
$$
where $\varepsilon(u)$ is the infinitesimal (linearized) strain tensor. The
coupled thermal equation is
$$
c_V\frac{\partial\Theta}{\partial t}
+\kappa T_0\operatorname{tr}\!\left(
\frac{\partial\varepsilon(u)}{\partial t}\right)
-k\nabla^2\Theta=0.
$$
Use $T_0=293\,\mathrm K$, $E=70000\,\mathrm{MPa}$, $\nu=0.3$,
$\rho=2700\,\mathrm{kg/m^3}$, $\alpha=2.31\times10^{-5}\,\mathrm K^{-1}$,
$c_V=\rho(910\,\mathrm{J/(kg\,K)})=2.457\,\mathrm{MPa/K}$, and
$k=237\times10^{-6}\,\mathrm{MPa\,m^2/(s\,K)}$,
$\mu=E/[2(1+\nu)]$,
$\lambda=E\nu/[(1+\nu)(1-2\nu)]$, and
$\kappa=\alpha(3\lambda+2\mu)$.

**Initial conditions:**

Set $u=0$ and $\Theta=0$ at $t_0=10\,\mathrm s$.

**Mechanical boundary conditions:**

On the symmetry boundary $x=0$, impose $u_x=0$ and zero tangential traction
$(\sigma n)\cdot\boldsymbol e_y=0$. On the symmetry boundary $y=0$, impose
$u_y=0$ and $(\sigma n)\cdot\boldsymbol e_x=0$. Impose zero traction
$\sigma n=0$ on the circular-hole boundary
$x^2+y^2=(0.1\,\mathrm m)^2$ and on the exterior edges $x=1\,\mathrm m$ and
$y=1\,\mathrm m$, where $n$ is the outward unit normal of $\Omega$.

**Thermal boundary conditions:**

For $t>t_0$, impose $\Theta=10\,\mathrm K$ on the circular-hole boundary.
Insulate every other boundary by imposing $q\cdot n=0$, where
$q=-k\nabla\Theta$ is the heat flux.

**Numerical method:**

Use a monolithic mixed finite-element formulation with vector CG2 displacement
$u$ and scalar CG1 temperature increment $\Theta$. Advance both fields together
with implicit Euler over 100 nonuniform time steps at
$$
t_i=10^{1+3i/100}\,\mathrm s,\qquad i=0,\ldots,100,
$$
so that $t_0=10\,\mathrm s$, $t_{100}=10^4\,\mathrm s$, and the step from
$t_{i-1}$ to $t_i$ has size
$$
\Delta t_i=t_i-t_{i-1}
=\left(10^{1+3i/100}-10^{1+3(i-1)/100}\right)\,\mathrm s,\qquad
i=1,\ldots,100.
$$

**Output:**

Save final displacement $u$ and temperature increment $\Theta$ to
`solution.xdmf` and `solution.h5`.
