# Problem_43: Axisymmetric linear elasticity of a pressurized hollow sphere

**Geometry and mesh:**

Let
$$
\omega=\{(r,z):r\geq0,\ z\geq0,\ R_i^2<r^2+z^2<R_e^2\},
\qquad R_i=9,\quad R_e=11,
$$
be the quarter-annular meridional computational section. Revolving $\omega$
about the $z$-axis gives the upper half of a hollow spherical shell; symmetry
across $z=0$ represents the complete shell. Use the accompanying `mesh.xdmf`
triangular mesh.

**Model:**

Find the axisymmetric displacement
$\boldsymbol u=u_r(r,z)\boldsymbol e_r+u_z(r,z)\boldsymbol e_z$ satisfying the
small-strain isotropic linear-elasticity equations
$$
\frac{\partial\sigma_{rr}}{\partial r}
+\frac{\partial\sigma_{rz}}{\partial z}
+\frac{\sigma_{rr}-\sigma_{\theta\theta}}{r}=0,
\qquad
\frac{\partial\sigma_{rz}}{\partial r}
+\frac{\partial\sigma_{zz}}{\partial z}
+\frac{\sigma_{rz}}{r}=0
\quad\text{in }\omega,
$$
with the axisymmetric strain tensor and the three-dimensional isotropic
constitutive law
$$
\varepsilon(u)=
\begin{bmatrix}
u_{r,r}&0&\tfrac12(u_{r,z}+u_{z,r})\\
0&u_r/r&0\\
\tfrac12(u_{r,z}+u_{z,r})&0&u_{z,z}
\end{bmatrix}_{(\boldsymbol e_r,\boldsymbol e_\theta,\boldsymbol e_z)},
\qquad
\sigma(u)=\lambda\operatorname{tr}(\varepsilon(u))I+2\mu\varepsilon(u),
$$
Here $\varepsilon(u)$ is the axisymmetric infinitesimal (linearized) strain
tensor. Use $E=10^5$, $\nu=0.3$, $\mu=E/[2(1+\nu)]$, and
$\lambda=E\nu/[(1+\nu)(1-2\nu)]$. There is no body force.

**Boundary conditions:**

Impose $u_z=0$ on $z=0$, $u_r=0$ on $r=0$, traction
$\sigma(u)n=-pn$ with $p=10$ on $r^2+z^2=R_e^2$, and zero traction on
$r^2+z^2=R_i^2$, where $n$ is the outward unit normal of the elastic body.

**Numerical method:**

Use continuous quadratic vector Lagrange elements and the axisymmetric weak
form with radial integration weight $r$.

**Output:**

Save the displacement $u$ to `solution.xdmf` and `solution.h5`.
