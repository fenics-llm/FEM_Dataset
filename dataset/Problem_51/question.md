# Problem_51: Clamped Reissner-Mindlin plate with selective reduced integration

**Geometry and mesh:**

Let $\Omega=(0,1)^2$ and use the accompanying `mesh.xdmf` structured
$50\times50$ bilinear quadrilateral mesh.

**Model:**

Find the transverse deflection $w:\Omega\to\mathbb R$ and rotation
$\theta:\Omega\to\mathbb R^2$. Define
$$
\chi(\theta)=\operatorname{sym}\nabla\theta,
\qquad
M(\theta)=D\left[(1-\nu)\chi(\theta)
+\nu\operatorname{tr}(\chi(\theta))I\right],
$$
and
$$
Q(w,\theta)=F(\theta-\nabla w),
\qquad
D=\frac{Et^3}{12(1-\nu^2)},
\qquad
F=\frac{5Et}{12(1+\nu)}.
$$
Use
$$
E=1000,\qquad \nu=0.3,\qquad t=10^{-3},\qquad f=-t^3=-10^{-9}.
$$
Solve
$$
-\nabla\cdot M(\theta)+Q(w,\theta)=0
\qquad\text{in }\Omega,
$$
and
$$
\nabla\cdot Q(w,\theta)=f
\qquad\text{in }\Omega.
$$

**Boundary conditions:**

Fully clamp the plate:
$$
w=0,\qquad \theta=0
\qquad\text{on }\partial\Omega.
$$

**Numerical method:**

Use continuous bilinear quadrilateral elements for both $w$ and $\theta$.
Use full quadrature for the bending and load terms and quadrature degree zero,
corresponding to one Gauss point per cell, for the shear term.

**Output:**

Save the transverse deflection $w$ and rotation $\theta$ to `solution.xdmf`
and `solution.h5`.
