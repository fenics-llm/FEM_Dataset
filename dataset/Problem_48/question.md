# Problem_48: Thermally degraded beam cross-section stiffness

**Geometry and mesh:**

Use the accompanying `box_girder_bridge.xdmf` mesh of the cross-section
$A\subset\mathbb{R}^2$ in the $(y,z)$ plane, with coordinates measured in
metres. Let $n$ be the outward unit normal. The top and bottom boundaries are
the exterior facets at $z=\max_A z$ and $z=\min_A z$, respectively.

**Temperature problem:**

Find the steady temperature $T$ satisfying
$$
-\nabla^2T=0\qquad\text{in }A.
$$
Set $T=20\,^{\circ}\mathrm C$ on the top boundary and
$T=800\,^{\circ}\mathrm C$ on the bottom boundary. Impose
$\nabla T\cdot n=0$ on every other exterior facet.

Define the temperature-dependent stiffness-reduction factor and elastic
moduli by
$$
r(T)=\exp\left[-\frac{T-20\,^{\circ}\mathrm C}
{211\,^{\circ}\mathrm C}\right],\qquad
E(T)=50r(T)\,\mathrm{GPa},\qquad
\mu(T)=\frac{E(T)}{2(1+\nu)},\qquad \nu=0.2.
$$

**Section quantities and shear-warping problem:**

Fix the beam axis at the geometric centroid
$$
|A|=\int_A dA,
\qquad
z_0=\frac{1}{|A|}\int_A z\,dA.
$$
Using the temperature-dependent Young modulus, define
$$
H_N=\int_A E(T)\,dA,
\qquad
H_{NM}=\int_A E(T)(z-z_0)\,dA,
$$
$$
H_M=\int_A E(T)(z-z_0)^2\,dA,
\qquad
\Delta_H=H_N H_M-H_{NM}^2.
$$
For the unit vertical shear resultant $V_z=1\,\mathrm{GN}$, define
$$
f(y,z)=
\frac{E(T)}{\Delta_H}
\left[H_N(z-z_0)-H_{NM}\right]V_z.
$$
Let $e_x$ be the unit vector along the beam axis. Find the scalar axial
warping amplitude $u(y,z)$, defined by the cross-section warping displacement
$$
\boldsymbol{u}^{\mathrm{warp}}(y,z)=u(y,z)e_x.
$$
Only this axial displacement component is unknown in the cross-section
problem. It satisfies
$$
-\nabla\cdot\bigl(\mu(T)\nabla u\bigr)=f
\qquad\text{in }A,
$$
with
$$
\mu(T)\nabla u\cdot n=0
\qquad\text{on }\partial A,
\qquad
\int_Au\,dA=0.
$$

After solving for $T$ and $u$, derive the following postprocessing quantities;
they are not additional unknown fields:
$$
\sigma_{xy}=\mu(T)\frac{\partial u}{\partial y},
\qquad
\sigma_{xz}=\mu(T)\frac{\partial u}{\partial z},
$$
$$
U_1=\frac12\int_A\mu(T)|\nabla u|^2\,dA,
\qquad
\overline H_V=\int_A\mu(T)\,dA,
$$
$$
\kappa_z=\frac{V_z^2}{2U_1\overline H_V},
\qquad
H_V=\kappa_z\overline H_V.
$$

**Numerical method:**

Use scalar CG1 elements for $T$. Solve the warping problem using a mixed
CG1/Real formulation for $(u,\lambda)$, where the Real-valued Lagrange
multiplier $\lambda$ enforces the zero-mean constraint.

**Output:**

Save only the primitive temperature $T$ and scalar axial warping amplitude $u$
to `solution.xdmf` and `solution.h5`. The shear stresses and section-stiffness
quantities are derived postprocessing outputs, not canonical solution fields.
