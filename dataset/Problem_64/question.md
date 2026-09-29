# Problem_64: Tensioned Circular Membrane under Gaussian Pressure

**Geometry and mesh:**

Let

$$
\Omega=\{(x,y):x^2+y^2<R^2\},
\qquad R=0.3.
$$

Use the supplied circular CSG mesh generated with resolution 40.

**Model:**

Find the transverse membrane deflection $w$ satisfying

$$
-T\nabla^2w=p(x,y)
\quad\text{in }\Omega,
$$

where $T=10$ and

$$
p(x,y)=4\exp\left[
-\frac12\left(
\frac{(x-x_0)^2+(y-y_0)^2}{\sigma^2}
\right)
\right],
$$

with

$$
\theta=0.2,\qquad
(x_0,y_0)=0.6R(\cos\theta,\sin\theta),\qquad
\sigma=0.025.
$$

**Boundary conditions:**

$$
w=0
\quad\text{on }\partial\Omega.
$$

**Numerical method:**

Use continuous piecewise-linear Lagrange elements and solve with conjugate
gradients and an ILU preconditioner.

**Output:**

Save $w$ to `solution.xdmf` and `solution.h5`.
