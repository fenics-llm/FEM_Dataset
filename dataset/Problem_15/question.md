# Problem_15: Solid Mechanics Problem 15 (Hard)

Source: `main-v4.tex`, Benchmark Problems, label `sm_q15`.

Original benchmark category: Solid Q15.

Reference solution source: `Results_ALL-FEM-main/reference solutions/solid/15/solid15.py`.

**Geometry:**

Let $\Omega = (0, 1.20) \times (0, 0.20)$ m be a rectangular strip with three circular holes of radius $a = 0.03$ m, centered on the midline $y = 0.10$ m at $x = 0.30$ m, $0.60$ m, and $0.90$ m.

**Model:**

Finite-strain Saint-Venant-Kirchhoff elasticity. Assume plane strain. The unknown is the displacement field $u = (u_x, u_y)$.

Define:

$F = I + \nabla u$ (deformation gradient)

$E = 0.5 (F^T F - I)$ (Green-Lagrange strain)

$S = \lambda \text{tr}(E) I + 2 \mu E$ (second Piola-Kirchhoff stress)

$J = \det(F)$, and Cauchy stress $\sigma = (1/J) F S F^T$

**Material:**

Use Lamé parameters: $\lambda = 5.769$ MPa, $\mu = 3.846$ MPa.

**Boundary conditions and loading:**

Left edge ($x = 0$): fixed, $u_x = 0, u_y = 0$.

Right edge ($x = 1.20$): prescribed displacement, $u_x = +0.012$ m, $u_y = 0$.

Top ($y = 0.20$) and bottom ($y = 0$): traction-free.

Hole boundaries: traction-free.

**Output:**

Save a plot of the deformed configuration as `q15\_def.png`.

Compute principal Green-Lagrange strains; save a color map of the maximum principal value $E_{\max}$ as `q15\_Emax.png`.

Define $s = S - (1/3) \text{tr}(S) I$ and $\sigma_{\text{vm}}(s) = \sqrt{1.5 (s:s)}$; save a color map as `q15\_vmS.png`.

Export the final displacement $u$ and $E_{\max}$ in XDMF format (`q15\_u.xdmf`, `q15\_Emax.xdmf`).

**Notes:**

Use load stepping with Newton iterations. Use load stepping and stop when max principal Green-Lagrange strain $E_{\max} \le 0.03$; if $E_{\max}$ exceeds $0.03$, do not advance the load.
