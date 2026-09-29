# Problem_50: Hertz contact with a rigid spherical indenter

**Geometry and mesh:**

Use the accompanying `mesh.xdmf` graded hexahedral mesh of the quarter block
$\Omega=[0,1]^2\times[-1,0]$. The undeformed top surface is the flat plane
$z=0$; the mesh is graded toward $x=y=z=0$ but contains no indenter-shaped
depression.

**Model:**

Find the displacement $u=[u_x,u_y,u_z]^T$ satisfying
$$
-\nabla\!\cdot\sigma(u)=0\qquad\text{in }\Omega,
$$
where $\varepsilon(u)$ is the infinitesimal strain tensor and
$$
\sigma(u)=\lambda\operatorname{tr}(\varepsilon(u))I+2\mu\varepsilon(u),
\qquad
\mu=\frac{E}{2(1+\nu)},\qquad
\lambda=\frac{E\nu}{(1+\nu)(1-2\nu)},
$$
with $E=10$ and $\nu=0.3$. A rigid sphere of radius $R=0.5$ is represented,
without meshing it, by the paraboloidal obstacle height
$$h(x,y)=-d+\frac{x^2+y^2}{2R}$$
for indentation depth $d=0.02$. On the top surface define the gap
$g=h-u_z$ and frictionless penalty pressure
$$p=k\langle-g\rangle_+=k\langle u_z-h\rangle_+,\qquad k=10^4,$$
where $\langle s\rangle_+=\max(s,0)$.

**Boundary conditions:**

Clamp $u=0$ on $z=-1$. On the symmetry face $x=0$, impose $u_x=0$ and zero
tangential traction components. On the symmetry face $y=0$, impose $u_y=0$
and zero tangential traction components. On the top face $z=0$, impose the
frictionless penalty-contact traction $\sigma n=-p\boldsymbol e_z$, so its
tangential components vanish. Impose zero traction on the exterior faces
$x=1$ and $y=1$.

**Numerical method:**

Use vector CG1 elements. Find $u_h$ satisfying the essential conditions and
$$
\int_\Omega\sigma(u_h):\varepsilon(v_h)\,dx
+\int_{z=0}k\langle u_{h,z}-h\rangle_+v_{h,z}\,ds=0
$$
for every admissible vector-CG1 test function $v_h$, using Newton's method.

**Output:**

Save displacement $u$ to `solution.xdmf` and `solution.h5`.
