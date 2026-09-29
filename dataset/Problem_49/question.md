# Problem_49: Linear buckling of a three-dimensional solid beam

**Geometry and mesh:**

Let $\Omega=(0,1)\times(0,0.01)\times(0,0.03)$ and use the
accompanying `mesh.xdmf` tetrahedral mesh.

**Model:**

Use small-strain isotropic linear elasticity with $E=1000$, $\nu=0$,
$$
\varepsilon(w)=\tfrac12(\nabla w+\nabla w^T),\qquad
\sigma(w)=\lambda_L\operatorname{tr}(\varepsilon(w))I
            +2\mu\varepsilon(w),
$$
where $\mu=E/[2(1+\nu)]$ and
$\lambda_L=E\nu/[(1+\nu)(1-2\nu)]$.

First find the reference pre-stress displacement $u_0$ from
$$
-\nabla\!\cdot\sigma(u_0)=0\qquad\text{in }\Omega,
$$
under the unit compressive traction specified below, and set
$\sigma_0=\sigma(u_0)$.

Then find the three smallest positive load factors $\lambda_i$ and nonzero
buckling modes $\phi_i$ satisfying
$$
-\nabla\!\cdot\left[\sigma(\phi_i)
  +\lambda_i(\nabla\phi_i)\sigma_0\right]=0
\qquad\text{in }\Omega.
$$
Equivalently, for every admissible test function $v$,
$$
\int_\Omega \sigma(\phi_i):\varepsilon(v)\,dx
+\lambda_i\int_\Omega
\sigma_0:\left[(\nabla\phi_i)^T\nabla v\right]dx=0.
$$
Here $\lambda_i$ scales the unit reference compression, while $\phi_i$ is
the infinitesimal displacement direction in which the pre-stressed beam first
loses stability.

**Boundary conditions for the pre-stress problem:**

On the left end $\Gamma_L=\{x=0\}$, impose $u_0=0$. On the right end
$\Gamma_R=\{x=1\}$, impose $(u_{0y},u_{0z})=(0,0)$ and the axial traction
$$
e_x\cdot\sigma(u_0)n=-1.
$$
The constrained transverse traction components on $\Gamma_R$ are reactions.
On the four lateral faces $\Gamma_{\rm lat}$, impose
$\sigma(u_0)n=0$.

**Boundary conditions for each buckling mode:**

Use the homogeneous counterparts of the displacement constraints:
$$
\phi_i=0\quad\text{on }\Gamma_L,
\qquad
(\phi_{iy},\phi_{iz})=(0,0)\quad\text{on }\Gamma_R.
$$
On unconstrained components, impose the natural incremental-traction
conditions
$$
\left[\sigma(\phi_i)+\lambda_i(\nabla\phi_i)\sigma_0\right]n=0
\quad\text{on }\Gamma_{\rm lat},
$$
and
$$
e_x\cdot\left[\sigma(\phi_i)
+\lambda_i(\nabla\phi_i)\sigma_0\right]n=0
\quad\text{on }\Gamma_R.
$$

**Numerical method:**

Use vector CG2 elements for $u_0$ and the buckling modes. Assemble the elastic
stiffness $K$ and $K_G=-G(\sigma_0)$, where $G(\sigma_0)$ discretizes the
second integral above, and solve
$$K\Phi_i=\lambda_iK_G\Phi_i.$$

**Output:**

Save buckling modes $\phi_1,\phi_2,\phi_3$ to `solution.xdmf` and
`solution.h5`. Eigenvector amplitude and sign are not prescribed; comparisons
must be invariant to multiplication of each mode by a nonzero scalar.
