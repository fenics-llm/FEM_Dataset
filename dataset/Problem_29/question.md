# Problem_29: Fluid Mechanics Problem 13 (Medium)

Source: `main-v4.tex`, Benchmark Problems, label `fm_q13`.

Original benchmark category: Fluid Q13.

Reference solution source: `Results_ALL-FEM-main/reference solutions/fluid/13/fluid13.py`.

**Geometry:**

Solve for the transport of a chemical in a rectangular channel $\Omega = (0, L) \times (0, H)$ with length $L = 1.0$ m and height $H = 0.10$ m.
Boundary partition: $\Gamma_{\text{in}} = \{0\}\times[0, H]$, $\Gamma_{\text{out}} = \{L\}\times[0, H]$, $\Gamma_w = (0, L)\times\{0\} \cup (0, L)\times\{H\}$.

**Mesh:**

Use a uniform mesh composed of $100 \times 10$ elements.

**Model:**

Steady advection-diffusion for a scalar concentration $c$ in $\Omega$ with a given velocity field $u(x,y) = (u_x(y), 0)$, \; $u_x(y) = U_{\text{max}} \cdot [4 y (H - y) / H^2]$.

**Boundary conditions:**

Inlet ($\Gamma_{\text{in}}$): $c = 0$.

Outlet ($\Gamma_{\text{out}}$): $c = 1$.

Walls ($\Gamma_w$): zero diffusive normal flux.

**Parameters:**

$L = 1.0$ m; $H = 0.10$ m; $U_{\text{max}} = 0.75$ m s$^{-1}$.

Diffusivity of the chemical: $D = 1.0\times10^{-5}$ m$^2$ s$^{-1}$.

**Output:**

Save the resulting concentration field output in xdmf format.

**Hint:**

Check the Peclet number (using the mesh size as length scale) to determine whether a stabilization scheme is required. If stabilization is required, include it in your formulation.
