# Problem_5: Solid Mechanics Problem 5 (Easy)

Source: `main-v4.tex`, Benchmark Problems, label `sm_q5`.

Original benchmark category: Solid Q5.

Reference solution source: `Results_ALL-FEM-main/reference solutions/solid/5/gpt5-solid5.py`.

**Geometry:**

Let $\Omega = (0, 1.0) \times (0, 0.20)$ meter be a rectangular plate with a rectangular notch cut from the left edge defined by the region $(0, 0.06) \times (0.08, 0.12)$.

**Model:**

Plane-stress linear elasticity for displacement $u = (u_x, u_y)$ in $\Omega$.

$\sigma$ is the Cauchy stress tensor and $n$ is the outward unit normal.

**Material:**

Young's modulus $E = 200$ GPa, Poisson's ratio $\nu = 0.30$.

**Boundary conditions and loads:**

Left outer boundary ($x = 0$): fixed, $u_x = 0, u_y = 0$.

Right edge ($x = 1.0$): $\sigma n = (2 \text{ MPa} \cdot \text{m}, 0)$.

Top and bottom edges, and notch boundaries: traction-free, $\sigma n = 0$.

**Output:**

Compute the von Mises equivalent stress and save a color map as `q5\_vm.png`.

Save the resulting displacement field in XDMF format.
