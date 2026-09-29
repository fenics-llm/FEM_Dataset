# Problem_59: Diffusion-Reaction Equation with Mixed Boundary Conditions

**Geometry and mesh:**

Let $\Omega=(0,1)$ with a uniform mesh of 10 linear elements.

**Model:**

Find $u:\Omega\to\mathbb{R}$ satisfying

$$
-u''+u=x
\quad\text{in }\Omega.
$$

**Boundary conditions:**

$$
u(0)=0,
\qquad
u'(1)=1-\frac{\cosh(1)}{\sinh(1)}.
$$

**Numerical method:**

Use continuous piecewise-linear Lagrange finite elements and include the
Neumann datum through the boundary term in the weak formulation.

**Output:**

Save $u$ to `solution.xdmf` and `solution.h5`.
