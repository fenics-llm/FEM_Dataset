# Problem_6: Solid Mechanics Problem 6 (Easy)

Source: `main-v4.tex`, Benchmark Problems, label `sm_q6`.

Original benchmark category: Solid Q6.

Reference solution source: `Results_ALL-FEM-main/reference solutions/solid/6/oss-multi-solid6.py`.

**Geometry:**

Let $\Omega = (0, 1.0) \times (0, 0.20)$ meter be a rectangular plate with a semicircular notch of radius $a = 0.05$ m removed from the top edge. The center of this semicircular notch is $[0.5, 0.2]$.

**Model:**

Plane-stress linear elasticity for displacement $u = (u_x, u_y)$ in $\Omega$.

$\sigma$ is the Cauchy stress tensor and $n$ the outward unit normal.

**Material:**

Young's modulus $E = 200$ GPa, Poisson's ratio $\nu = 0.30$.

**Boundary conditions:**

Bottom edge ($y = 0$): fixed, $u_x = 0, u_y = 0$.

Top edge ($y = 0.20$, on $x$ from $[0, 0.45) \cup (0.55, 1.0]$): $\sigma n = (0, -10 \text{ MPa} \cdot \text{m})$.

Both vertical sides ($x = 0$ and $x = 1$) and the notch arc: Traction-free, $\sigma n = (0, 0)$.

**Output:**

Compute the von Mises equivalent stress (plane-stress) and save a color map as `q6\_vm.png`.

Save the resulting displacement field in XDMF format.
