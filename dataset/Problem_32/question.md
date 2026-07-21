# Problem_32: Multiphysics Problem 1 (Hard)

Source: `main-v4.tex`, Benchmark Problems, label `mf_q1`.

Original benchmark category: Multiphysics Q1.

Reference solution source: `None`.

**Geometry:**

Solve the transport of a chemical inside an expanding circular disk $\Omega(t)$ with radius $R(t)$, where $R(t) = R_0 + s \cdot t$ with constant rate $s = 1.0 \times 10^{-4}$ m s$^{-1}$. Here, $R_0$ is the initial radius and is equal to $0.05$ m. The boundary of this circular disk is denoted by $\Gamma(t)$.

**Mesh:**

Unstructured triangular mesh of $\Omega(0)$ with characteristic size $h_0 \approx 1.0 \times 10^{-3}$ m.

The mesh motion is given by $w(x,t) = s x / \|x\|$ for $x \neq (0,0)$, $w(x,t) = 0$ for $x = 0$.

**Model:**

Solve the diffusion-reaction equation for concentration $c$ on a moving domain using an Arbitrary Lagrangian-Eulerian (ALE) description, with constant diffusivity $D$ and first-order decay equal to $\kappa c$, where $\kappa$ is the decay rate and is positive.

**Boundary conditions:**

Zero total flux at the moving boundary: $(-D\nabla c - w c) \cdot n = 0$ on $\Gamma(t)$ for all $t \ge 0$, where $D$ is the diffusion coefficient.

Initial condition: $c(x,0) = 1$ for $x \in \Omega(0)$.

**Parameters:**

Diffusivity: $D = 1.0 \times 10^{-5}$ m$^2$ s$^{-1}$.

Decay rate: $\kappa = 1.0 \times 10^{-4}$ s$^{-1}$.

Use the time step $\Delta t = 0.01$ s. Simulate up to 10 seconds.

**Output:**

Save the output concentration as function of time in an XDMF file.

Report the total concentration in the domain after every 100 time steps.
