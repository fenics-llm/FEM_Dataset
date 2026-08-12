# Problem 32: Conservative ALE Weak Form

Let \(\dot c = \partial c/\partial t\vert_X\) denote the time derivative at
fixed ALE (mesh) coordinate \(X\), and let \(w\) be the mesh velocity.  Thus,
at a fixed physical position \(x\),

\[
\partial_t c = \dot c - w\cdot\nabla c.
\]

The physical conservative balance and the benchmark's relative-flux boundary
condition are

\[
\partial_t c - \nabla\cdot(D\nabla c) + \kappa c = 0
\quad\text{in }\Omega(t),
\qquad
(-D\nabla c - wc)\cdot n = 0
\quad\text{on }\Gamma(t).
\]

Equivalently, the normal diffusive flux is

\[
D\nabla c\cdot n = -wc\cdot n.
\]

This is important: homogeneous diffusive Neumann data
\(D\nabla c\cdot n=0\) would be a different problem and would not account
for the expanding boundary.

## Weak form with the `wc` term explicit

For every \(v\in H^1(\Omega(t))\), multiply the strong balance by \(v\) and
integrate the diffusion term by parts.  Substitution of the relative-flux
condition gives: find \(c(t)\in H^1(\Omega(t))\) such that

\[
\boxed{
\int_{\Omega(t)} v\dot c\,dx
- \int_{\Omega(t)} v\,w\cdot\nabla c\,dx
+ \int_{\Omega(t)} D\nabla c\cdot\nabla v\,dx
+ \int_{\Omega(t)} \kappa c v\,dx
+ \int_{\Gamma(t)} v\,c\,w\cdot n\,ds = 0.
}
\tag{WF-32a}
\]

The final boundary integral is the required \(wc\) contribution.  Its plus
sign follows from

\[
-\int_{\Gamma(t)}vD\nabla c\cdot n\,ds
= +\int_{\Gamma(t)}v\,wc\cdot n\,ds.
\]

No Dirichlet boundary condition is imposed.

## Equivalent volume-only form

Integrate the mesh-advection term in (WF-32a) by parts:

\[
-\int_{\Omega}v\,w\cdot\nabla c\,dx
+\int_{\Gamma}vcw\cdot n\,ds
= \int_{\Omega}c\,w\cdot\nabla v\,dx
+\int_{\Omega}cv\,\nabla\cdot w\,dx.
\]

Therefore the equivalent, boundary-term-free form is

\[
\boxed{
\int_{\Omega(t)} v\dot c\,dx
+\int_{\Omega(t)} D\nabla c\cdot\nabla v\,dx
+\int_{\Omega(t)} \kappa cv\,dx
+\int_{\Omega(t)} c\,w\cdot\nabla v\,dx
+\int_{\Omega(t)} cv\,\nabla\cdot w\,dx = 0.
}
\tag{WF-32b}
\]

Use exactly one of (WF-32a) and (WF-32b), never both: adding the boundary
term to (WF-32b) double-counts the mesh-flux contribution.

## Backward-Euler residual

After moving the mesh to \(\Omega^{n+1}\), let \(c^n_{\rm ALE}\) be the
previous concentration transported with its mesh nodes onto that new mesh.
With \(c\equiv c^{n+1}\), a direct FEniCS representation of (WF-32a) is

```python
# WF-32a: c_old is c^n_ALE represented on the current, moved mesh.
F = (
    v * (c - c_old) / dt * dx
    - v * dot(w, grad(c)) * dx
    + D * dot(grad(c), grad(v)) * dx
    + kappa * c * v * dx
    + v * c * dot(w, n) * ds
)
```

The corresponding volume-only residual from (WF-32b) is

```python
# WF-32b: mathematically equivalent to WF-32a.
F = (
    v * (c - c_old) / dt * dx
    + D * dot(grad(c), grad(v)) * dx
    + kappa * c * v * dx
    + c * dot(w, grad(v)) * dx
    + c * v * div(w) * dx
)
```

For this prescribed radial mesh motion, the explicit-boundary form is the
safer initial implementation: it does not require evaluating
\(\nabla\cdot w=s/\|x\|\), which is singular in its pointwise expression at
the disk centre even though the weak integral is well-defined.

## Conservation check

Set \(v=1\) in (WF-32b).  The result is the ALE transport theorem applied to
the total chemical amount \(M(t)=\int_{\Omega(t)}c\,dx\):

\[
\frac{dM}{dt}=-\kappa M,
\qquad
M(t)=M(0)e^{-\kappa t}.
\]

This is the primary validation quantity.  Domain expansion lowers a uniform
concentration while preserving total chemical amount when \(\kappa=0\); it
must not create a spurious gain or loss of chemical.
