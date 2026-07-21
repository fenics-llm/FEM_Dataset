# Problem_33: Multiphysics Problem 2 (Hard)

Source: `main-v4.tex`, Benchmark Problems, label `mf_q2`.

Original benchmark category: Multiphysics Q2.

Reference solution source: `Results_ALL-FEM-main/reference solutions/multiphysics - Copy/2/multi-2.py`.

**Geometry:**

Unit square $\Omega = (0, 1) \times (0, 1)$.

**Mesh:**

Use a structured mesh with $200 \times 200$ elements.

**Model:**

Solve in $\Omega$ the Allen-Cahn curvature-flow equation for the phase field $\phi$, given by $\partial\phi/\partial t = -M ( (1/\varepsilon) W'(\phi) - \varepsilon \nabla^2\phi )$, where the double-well potential is $W(\phi) = 0.25 (\phi^2 - 1)^2$. Solve this equation using the finite element method. Use an implicit time stepping algorithm.

**Boundary conditions:**

Impose homogeneous Neumann conditions on the entire boundary, namely $\nabla\phi \cdot n = 0$ on $\partial\Omega$.

Initialize the phase field using $\phi(x, y, 0) = \tanh(d_{\text{rect}}(x, y) / (\sqrt{2} \varepsilon))$, where $d_{\text{rect}}$ is a signed distance to the surface of a square that is centered at $(0.5, 0.5)$ with side length equal to $0.5$. The signed distance function is negative inside the square and positive outside.

**Parameters:**

Simulate up to the final time $T = 0.20$. Use an interface thickness $\varepsilon = 0.01$ and a mobility $M = 1.0$. Choose a time step $\Delta t = 1.0\text{e}-3$.

**Output:**

Save the phase field, $\phi$, at each of the specified time steps: $t = 0.00, 0.05, 0.10$, and $0.20$. Store these files in XDMF format.
