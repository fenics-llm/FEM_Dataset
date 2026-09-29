# Problem_40: Low-Womersley Oscillatory Channel Flow (Medium)

**Geometry and mesh:**

Let

$$
\Omega=(0,L)\times(-a,a),
\qquad
L=2.0\ \mathrm{m},
\qquad
a=0.50\ \mathrm{m}.
$$

Use the supplied uniform triangular mesh corresponding to a $128\times64$
crossed subdivision of $\Omega$.

**Model:**

Solve for the velocity $u=(u_x,u_y)$ and pressure $p$:

$$
\rho\frac{\partial u}{\partial t}
-\mu\nabla^2u
+\nabla p
=G_0\cos(\omega t)e_x,
\qquad
\nabla\cdot u=0
\quad\text{in }\Omega\times(0,6\pi].
$$

Use

$$
\rho=1.0\ \mathrm{kg\,m^{-3}},
\qquad
\mu=1.0\ \mathrm{Pa\,s},
\qquad
G_0=1.0\ \mathrm{N\,m^{-3}},
\qquad
\omega=1.0\ \mathrm{s^{-1}}.
$$

**Boundary and initial conditions:**

Impose streamwise periodicity,

$$
u(0,y,t)=u(L,y,t),
\qquad
p(0,y,t)=p(L,y,t),
$$

no slip at $y=\pm a$,

$$
u(x,\pm a,t)=0,
$$

and the pressure normalization

$$
\int_\Omega p\,dx=0.
$$

Initially,

$$
u(x,y,0)=0.
$$

**Time discretization:**

Use Crank--Nicolson with

$$
\Delta t=\frac{\pi}{100}\ \mathrm{s},
$$

evaluating the harmonic body force at the temporal midpoint of each step.

**Output:**

Save $u$ and $p$ to `solution.xdmf` at

$$
t\in
\left\{
4\pi,\frac{9\pi}{2},5\pi,\frac{11\pi}{2},6\pi
\right\}\ \mathrm{s}.
$$
