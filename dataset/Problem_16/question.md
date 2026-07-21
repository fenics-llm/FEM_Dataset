# Problem_16: Solid Mechanics Problem 16 (Hard)

Source: `main-v4.tex`, Benchmark Problems, label `sm_q16`.

Original benchmark category: Solid Q16.

Reference solution source: `None`.

**Geometry:**

Analyze a plate in plane strain using the symmetric quarter domain $[0, 100] \times [0, 180]$ mm with a quarter hole of radius $50$ mm at the origin.

**Model:**

Use small-strain von Mises elastoplasticity with associative flow and perfect plasticity.

Let $u = (u_x, u_y)$, $\varepsilon = 1/2 (\nabla u + \nabla u^T)$.

The strain tensor is split into two elastic and plastic strain as $\varepsilon_{ij} = \varepsilon^e_{ij} + \varepsilon^p_{ij}$, where $\varepsilon^e$ represents the elastic part and $\varepsilon^p$ the plastic part. The plastic part verifies $\varepsilon^p_{kk}=0$. Here, and in what follows, repeated indices indicate summation.

The quasi-static momentum balance is $\sigma_{ij,j} = 0$, where $\sigma_{ij} = (\lambda + 2/3 \mu) \varepsilon_{kk} \delta_{ij} + s_{ij}$.

Here, $s_{ij} = 2 \mu ( e_{ij} - e^p_{ij} )$, where $e_{ij} = \varepsilon_{ij} - (1/3) \varepsilon_{kk} \delta_{ij}$ and $e^p_{ij}$ is the plastic deviatoric strain.

The plastic flow is given by $\Delta e^p_{ij} = \Delta \varepsilon^p  (3/2)  s_{ij} / q$, where $q = \sqrt{ (3/2) s_{ij} s_{ij} }$ is the equivalent stress, and $\Delta \varepsilon^p$ is the increment of equivalent plastic strain (the plastic multiplier), defined as $\Delta \varepsilon^p = \sqrt{2/3 \Delta \varepsilon^p_{ij} \Delta \varepsilon^p_{ij}}$.

The yield condition is $F = q - \sigma_Y \le 0$ and the plastic loading follow the Kuhn-Tucker conditions, i.e., $\Delta \varepsilon^p \ge 0$, $F \le 0$, $\Delta \varepsilon^p F = 0$. The initial plastic strain is zero.

**Boundary conditions:**

Apply symmetry on $x = 0$ using $u_x = 0$ and zero traction in a direction tangent to $x=0$.

Apply symmetry on $y = 0$ using $u_y = 0$ and zero traction in a direction tangent to $y=0$.

At the top edge ($y = 180$), prescribe $u = (0, 1)$ mm.

Set the edge $x = 100$ mm and the hole boundary traction-free.

**Parameters:**

Use $\lambda = 19.44$ GPa and $\mu = 29.17$ GPa. Set $\sigma_Y = 243$ MPa.

**Output:**

Save the resulting displacement in an xdmf file.
