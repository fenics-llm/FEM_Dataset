# Problem 34: Mixed Cahn--Hilliard Formulation and Proposed Implementation

## Purpose and scope

This document specifies the weak form and proposed numerical implementation
for Problem 34 before a reference solver is written.  It follows the benchmark
statement in `question.md`.  The formulation uses two second-order equations,
so continuous Lagrange finite elements can be used instead of directly
discretising the fourth-order Cahn--Hilliard equation.

The approach is based on the mixed-field strategy in the
[legacy DOLFIN Cahn--Hilliard demo](https://fenics.readthedocs.io/projects/dolfin/en/2017.2.0/demos/cahn-hilliard/python/demo_cahn-hilliard.py.html),
but uses this problem's logarithmic chemical free energy and degenerate
mobility.

## Strong form

Let \(\Omega=(0,1)^2\).  The unknown concentration \(c\) and chemical
potential \(\mu\) satisfy

\[
\frac{\partial c}{\partial t}=\nabla\!\cdot\left(M(c)\nabla\mu\right),
\qquad
\mu=3\alpha\,\mu_c(c)-\nabla^2c,
\tag{1}
\]

where

\[
M(c)=c(1-c),
\qquad
\mu_c(c)=\frac{1}{2\theta}\ln\!\left(\frac{c}{1-c}\right)+1-2c,
\tag{2}
\]

with \(\theta=1.5\) and \(\alpha=3000\).  Both \(c\) and \(\mu\) are
periodic across opposite sides of the unit square.  The initial concentration
is

\[
c^0(x,y)=0.63+r(x,y),\qquad r\in[-0.05,0.05],\qquad \int_\Omega r\,dx=0.
\tag{3}
\]

Thus the initial concentration lies strictly inside \((0,1)\), as required
by the logarithm, and has mean \(\bar c=0.63\).

## Weak form

Let \(V_h\subset H^1_{\mathrm{per}}(\Omega)\) be a periodic continuous
piecewise-linear (P1) scalar space, and let

\[
W_h=V_h\times V_h.
\]

At time \(t_{n+1}=t_n+\Delta t\), given \(c_h^n\), find
\((c_h^{n+1},\mu_h^{n+1})\in W_h\) such that for every
\((q_h,v_h)\in W_h\),

\[
\begin{aligned}
R_c(c_h^{n+1},\mu_h^{n+1};q_h)
={}&
\int_\Omega \frac{c_h^{n+1}-c_h^n}{\Delta t}q_h\,dx
+\int_\Omega M(c_h^{n+1})\nabla\mu_h^{n+1}\cdot\nabla q_h\,dx=0,
\\
R_\mu(c_h^{n+1},\mu_h^{n+1};v_h)
={}&
\int_\Omega \mu_h^{n+1}v_h\,dx
-\int_\Omega 3\alpha\,\mu_c(c_h^{n+1})v_h\,dx
-\int_\Omega \nabla c_h^{n+1}\cdot\nabla v_h\,dx=0.
\end{aligned}
\tag{4}
\]

The residual is \(R=R_c+R_\mu\).  Periodicity cancels the boundary terms in
both integrations by parts.  Equation (4) is backward Euler: all nonlinear
terms, including \(M(c)\), are evaluated at \(t_{n+1}\).

Choosing \(q_h=1\) in the first equation shows that the discrete total mass
\(\int_\Omega c_h\,dx\) is conserved (up to nonlinear-solver tolerance).
This is a required verification metric.

## Newton linearisation

For Newton's method, define the increment
\((\delta c,\delta\mu)\in W_h\).  The Jacobian action of (4) is

\[
\begin{aligned}
J_c[\delta c,\delta\mu;q]
={}&
\int_\Omega \frac{\delta c}{\Delta t}q\,dx
+\int_\Omega \left[M'(c)\delta c\,\nabla\mu
+M(c)\nabla\delta\mu\right]\cdot\nabla q\,dx,
\\
J_\mu[\delta c,\delta\mu;v]
={}&
\int_\Omega \delta\mu\,v\,dx
-\int_\Omega 3\alpha\,\mu_c'(c)\delta c\,v\,dx
-\int_\Omega \nabla\delta c\cdot\nabla v\,dx,
\end{aligned}
\tag{5}
\]

where

\[
M'(c)=1-2c,
\qquad
\mu_c'(c)=\frac{1}{2\theta}\left(\frac{1}{c}+\frac{1}{1-c}\right)-2.
\tag{6}
\]

The implementation should construct (4) directly in UFL and use
`derivative(residual, w, dw)` for (5), rather than manually coding a separate
Jacobian.  This retains the exact derivative of the mobility--potential
coupling.

## Proposed FEniCS implementation

The future `solution.py` should keep the following objects visible in its main
solve section:

```python
# Periodic scalar P1 space and mixed concentration/potential space.
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, MixedElement([P1, P1]), constrained_domain=periodic_bc)

w = Function(W)       # current Newton iterate: (c, mu)
w_previous = Function(W)
dw = TrialFunction(W)
q, v = TestFunctions(W)
c, mu = split(w)
c_previous, _ = split(w_previous)

mobility = c*(1.0 - c)
mu_c = (1.0/(2.0*theta))*ln(c/(1.0 - c)) + 1.0 - 2.0*c

residual = (
    ((c - c_previous)/dt)*q*dx
    + dot(mobility*grad(mu), grad(q))*dx
    + mu*v*dx
    - 3.0*alpha*mu_c*v*dx
    - dot(grad(c), grad(v))*dx
)
jacobian = derivative(residual, w, dw)
```

`periodic_bc` must map the right boundary to the left boundary and the top
boundary to the bottom boundary, while identifying the upper-right corner
only once.  No Dirichlet boundary condition is applied.

Each accepted time step will:

1. Set `w_previous` to the prior converged mixed solution.
2. Use that same solution as the initial Newton iterate in `w`.
3. Solve `residual == 0` with `jacobian` using a damped Newton method and a
   visible direct linear solver (`mumps` for the initial 2-D reference run).
4. Accept the step only after Newton converges and the nodal concentration
   remains strictly in \((0,1)\).
5. Save the requested fields at \(t=0\), \(3\times10^{-6}\),
   \(10^{-4}\), \(10^{-3}\), and \(4\times10^{-2}\).

## Adaptive time stepping and logarithm admissibility

The logarithmic term is defined only for \(0<c<1\).  The solver must not
silently clip \(c\), because clipping changes the benchmark PDE and can break
mass conservation.  Instead:

- Begin with \(\Delta t=10^{-7}\) (inside the benchmark's suggested range).
- Use damped Newton updates.  A candidate update that makes any concentration
  degree of freedom nonpositive or at least one is rejected and retried with a
  smaller damping factor.
- If Newton fails or no admissible damped update is found, reject the entire
  time step, restore the previous converged state, and reduce \(\Delta t\),
  for example by a factor of two.
- After a comfortably converged step (six Newton iterations or fewer in the
  initial implementation), increase \(\Delta t\) by 50%.  For moderate
  convergence (seven through eleven iterations), increase it conservatively
  by 10%; decrease it when Newton takes many iterations (twelve or more
  initially).  Exact output times take precedence: cap a proposed step so it
  lands exactly on the next requested output time.

The first implementation should record accepted/rejected steps, Newton
iteration counts, minimum/maximum concentration, and total mass.  These
diagnostics are necessary to distinguish a physical spinodal pattern from a
time-step or nonlinear-solver artifact.

## Verification plan before acceptance

1. **Periodic field check:** inspect that opposite-side traces of both \(c\)
   and \(\mu\) agree.
2. **Mass conservation:** compare \(\int_\Omega c_h\,dx\) at every accepted
   step with its initial value.  The relative drift should be consistent with
   nonlinear-solver tolerance.
3. **Admissibility:** verify \(0<c_h<1\) throughout the simulation; report
   the global minimum and maximum at every saved time.
4. **Energy trend:** evaluate the corresponding discrete free-energy trend as
   a diagnostic.  The continuous model dissipates free energy.  A fully
   implicit solve is a useful baseline, but this non-convex logarithmic energy
   does not by itself give an unconditional discrete energy-dissipation proof;
   a material energy increase must therefore be investigated and may motivate
   a convex-splitting or other energy-stable time discretisation.
5. **Resolution/time-step sensitivity:** repeat at a finer mesh and/or tighter
   time-step adaptation.  Compare mass, concentration range, and qualitative
   domain morphology at the required output times.

The mesh resolution, nonlinear tolerances, damping implementation, maximum
time step, and output/checkpoint naming will be fixed only after this
formulation is reviewed.
