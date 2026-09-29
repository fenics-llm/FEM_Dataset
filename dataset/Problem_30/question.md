# Problem_30: Fluid Mechanics Problem 14 (Hard)

Source: Schäfer--Turek/DFG benchmark 2D-3 (Re = 100, fixed time interval).

Reference source: `DFG_2D-3_benchmark.html`.

**Geometry:**

Consider the channel

$$
\Omega=([0,2.2]\times[0,0.41])\setminus\overline{B_{0.05}(0.2,0.2)},
$$

containing a circular cylinder of diameter $D=0.1$ centered at $(0.2,0.2)$. Use a conforming body-fitted mesh refined near the cylinder and in its wake.

**Model:**

Solve the two-dimensional unsteady incompressible Navier--Stokes equations

$$
\frac{\partial u}{\partial t}-\nu\Delta u+(u\cdot\nabla)u+\nabla p=0,
\qquad \nabla\cdot u=0,
$$

using a residual-based Variational Multiscale formulation with SUPG, PSPG, and grad-div stabilization. State the finite-element spaces and stabilization parameters used.

**Boundary and initial conditions:**

At the inlet $x=0$, prescribe

$$
u(0,y,t)=\left(\frac{4U(t)y(0.41-y)}{0.41^2},0\right),
\qquad U(t)=1.5\sin\left(\frac{\pi t}{8}\right).
$$

Impose $u=0$ on the upper and lower channel walls and on the cylinder. At the outlet $x=2.2$, impose

$$
\nu\frac{\partial u}{\partial n}-pn=0.
$$

Use the initial condition $u(x,0)=(0,0)$.

**Parameters:**

Set density $\rho=1$ and kinematic viscosity $\nu=0.001$, giving a maximum Reynolds number $\mathrm{Re}=100$ based on the mean inlet velocity and cylinder diameter. Simulate over $0\leq t\leq8$. Use Crank--Nicolson time integration with $\Delta t=1/1600$ and verify that the reported quantities are sufficiently mesh-converged.

**Benchmark quantities:**

With $\sigma=\nu\nabla u-pI$ and $\eta$ pointing outward from the cylinder into the fluid, compute

$$
\begin{pmatrix}F_D\\F_L\end{pmatrix}=\int_{\partial B_{0.05}}\sigma\eta\,ds,
\qquad
C_D=\frac{2F_D}{U_{\mathrm{mean}}^2D},
\qquad
C_L=\frac{2F_L}{U_{\mathrm{mean}}^2D},
$$

where $U_{\mathrm{mean}}=1$ and $D=0.1$. Also compute

$$
\Delta p(t)=p(0.15,0.2,t)-p(0.25,0.2,t).
$$

Report $\max C_D(t)$, the time at which it occurs, $\max C_L(t)$, the time at which it occurs, and $\Delta p(8)$. Compare them with

$$
C_{D,\max}=2.950921575,
\qquad C_{L,\max}\approx0.47795,
\qquad \Delta p(8)\approx0.1116.
$$

**Output:**

Save the velocity and pressure at $t=8$ to `vms_solution.xdmf`. Save the time histories of $C_D$, $C_L$, and $\Delta p$ in a machine-readable text file.
