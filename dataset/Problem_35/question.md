# Problem_35: Multiphysics Problem 4 (Hard)

Source: `main-v4.tex`, Benchmark Problems, label `mf_q4`.

Original benchmark category: Multiphysics Q4.

Reference solution source: `None`.

**Geometry:**

Solve the coupled Stokes-Darcy problem on the rectangle $\Omega = (0, \pi) \times (-1, 1)$. The upper subdomain is the Stokes region $\Omega_S = (0, \pi) \times (0, 1)$ and the lower subdomain is the Darcy region $\Omega_D = (0, \pi) \times (-1, 0)$. The interface is $\Gamma = (0, \pi) \times \{0\}$ with unit normal $n = (0, -1)$ pointing from $\Omega_S$ into $\Omega_D$ and unit tangent $t = (1, 0)$.

**Model:**

In $\Omega_S$, the flow is incompressible Stokes with a body force: $\nabla \cdot u_S = 0$ and $-\nabla \cdot \sigma(u_S, p_S) = b$, with $\sigma = -p_S I + 2 \nu \text{sym}(\nabla u_S)$. Here, $\text{sym}(\nabla u_S)$ denotes the symmetric part of $\nabla u_S$.

In $\Omega_D$, the flow is incompressible Darcy: $\nabla \cdot u_D = 0$ and $u_D = -(k/\mu) \nabla p_D$.

On $\Gamma$, impose mass continuity $u_S \cdot n = u_D \cdot n$, normal traction balance $n \cdot \sigma n = -p_D/\rho$, and the Saffman tangential condition $(\alpha/\sqrt{k}) u_S\cdot t = - t \cdot \sigma n$.

Here, $u_S$ and $u_D$ are the fluid velocity in the Stokes and Darcy subdomains respectively, while $p_S$ and $p_D$ are the fluid pressure in the Stokes and Darcy subdomains, respectively.

**Boundary conditions:**

Dirichlet boundary conditions for both Stokes and Darcy boundaries.

Stokes (on $\partial\Omega_S = \{y=1\} \cup \{x=0\}\times\{y\in[0,1]\} \cup \{x=\pi\}\times\{y\in[0,1]\}$):

$u_S(x,y) = [ w'(y) \cos x , w(y) \sin x ]$.

Darcy ($\partial\Omega_D = \{y=-1\} \cup \{x=0\}\times\{y\in[-1,0]\} \cup \{x=\pi\}\times\{y\in[-1,0]\}$):

$p_D(x,y) = \rho g \exp(y) \sin x$.

$w(y) = -K - (g y)/(2\nu) + (K/2 - \alpha g/(4\nu^2)) y^2$ and $K = k \rho g / \mu$.

**Parameters:**

Take $g = 1, \rho = 1, \nu = 1, k = 1, \mu = 1, K = 1, \alpha = 1$.

Body force $b$:

$b_x(x,y) = [ (\nu K - (\alpha g)/(2\nu)) y - g/2 ] \cos(x)$

$b_y(x,y) = [ ( (\nu K)/2 - (\alpha g)/(4\nu) ) y^2 - (g/2) y + ( (\alpha g)/(2\nu) - 2\nu K ) ] \sin(x)$

**Numerical formulation:**

Use the supplied mesh. Assemble and solve one coupled finite-element variational system for $(u_S,p_S,p_D)$, using $[P_2]^2$ for $u_S$ and $P_1$ for both pressures. Represent these fields on the full mesh with zero extension outside their corresponding subdomains, without constraining the interface trace, and impose $p_S(0,0.5)=0$. The volume equations must be integrated only over their corresponding subdomains, and the three interface conditions must enter the coupled weak form. Compute $u_D=-(k/\mu)\nabla p_D$ from the solved Darcy pressure in a vector $DG_0$ space, extended by zero outside $\Omega_D$.

The manufactured expressions in this problem may be used only as boundary data and for validation. Directly interpolating or projecting them as the returned solution fields is not permitted.

**Output:**

Save $u_S$, $p_S$, $p_D$, and $u_D$ in `solution.xdmf` and `solution.h5`.
