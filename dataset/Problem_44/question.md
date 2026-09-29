# Problem_44: Sequential steady thermoelasticity of a heated beam

**Geometry and mesh:**

Let $\Omega=(0,5\,\mathrm m)\times(0,0.3\,\mathrm m)$ and use the accompanying
`mesh.xdmf` crossed triangular mesh.

**Model:**

Let $\Theta=T-T_{\mathrm{ref}}$ be the temperature increment from the
stress-free reference temperature $T_{\mathrm{ref}}$, not an absolute
temperature. This is one-way sequential coupling: $\Theta$ affects $u$ through
the thermoelastic stress, while $u$ does not enter the thermal equation.
First find the temperature increment $\Theta$ from $-\nabla^2\Theta=0$.
Then find the plane-strain displacement $u=[u_x,u_y]^T$ from
$-\nabla\cdot\sigma(u,\Theta)=f$, with density
$\rho=2400\,\mathrm{kg/m^3}$, gravitational acceleration
$g=9.81\,\mathrm{m/s^2}$, body-force density
$f=(0,-\rho g\,10^{-6})\,\mathrm{MPa/m}$, and
$$
\sigma(u,\Theta)=\lambda\operatorname{tr}(\varepsilon(u))I
+2\mu\varepsilon(u)-\alpha(3\lambda+2\mu)\Theta I,
$$
where $\varepsilon(u)$ is the infinitesimal (linearized) strain tensor,
$E=5\times10^4\,\mathrm{MPa}$, $\nu=0.2$,
$\alpha=10^{-5}\,\mathrm K^{-1}$,
$\mu=E/[2(1+\nu)]$, and $\lambda=E\nu/[(1+\nu)(1-2\nu)]$.

**Thermal boundary conditions:**

Set $\Theta=50\,\mathrm K$ on $y=0$ and $\Theta=0$ on
$y=0.3\,\mathrm m$, $x=0$, and $x=5\,\mathrm m$.

**Mechanical boundary conditions:**

Clamp $u=0$ on $x=0$ and $x=5\,\mathrm m$. Impose zero traction
$\sigma n=0$ on $y=0$ and $y=0.3\,\mathrm m$, where $n$ is the outward unit
normal.

**Numerical method:**

Solve the thermal problem with CG1, then the mechanical problem with vector
CG2.

**Output:**

Save displacement $u$ and temperature increment $\Theta$ to `solution.xdmf`
and `solution.h5`.
