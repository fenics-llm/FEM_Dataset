# Problem_41: High-Womersley Oscillatory Channel Flow (Hard)

**Geometry and mesh:**

$$
\Omega=(0,L)\times(-a,a),
\qquad
L=2.0\ \mathrm{m},
\qquad
a=0.50\ \mathrm{m}.
$$

Use the supplied uniform $64\times160$ crossed triangular mesh.

**Model:**

Solve for velocity $u=(u_x,u_y)$ and pressure $p$:

$$
\rho\frac{\partial u}{\partial t}
-\mu\nabla^2u
+\nabla p
=G_0\cos(\omega t)e_x,
\qquad
\nabla\cdot u=0
\quad\text{in }\Omega.
$$

Use

$$
\rho=1.0\ \mathrm{kg\,m^{-3}},
\qquad
\mu=1.0\ \mathrm{Pa\,s},
\qquad
G_0=1.0\ \mathrm{N\,m^{-3}},
\qquad
\omega=400\ \mathrm{s^{-1}}.
$$

**Boundary and initial conditions:**

$$
u(0,y,t)=u(L,y,t),
\qquad
p(0,y,t)=p(L,y,t),
$$

$$
u(x,\pm a,t)=0,
\qquad
\int_\Omega p\,dx=0,
\qquad
u(x,y,0)=0.
$$

**Time discretization:**

Let $P=2\pi/\omega$. Use Crank--Nicolson with $\Delta t=P/100$, evaluate
the body force at the temporal midpoint, and solve until $T=80P$.

**Output:**

Save $u$ and $p$ to `solution.xdmf` at

$$
\frac{t}{P}\in
\left\{
79,79+\frac14,79+\frac12,79+\frac34,80
\right\}.
$$
