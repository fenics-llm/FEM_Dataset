# Problem_29: Fluid Mechanics Problem 13 (Medium)

Source: `main-v4.tex`, Benchmark Problems, label `fm_q13`, revised to make the
outlet concentration layer a controlled, marginally resolved SUPG test.

Original benchmark category: Fluid Q13.

Reference solution source: `Results_ALL-FEM-main/reference solutions/fluid/13/fluid13.py`, subsequently rewritten for the revised parameters.

**Geometry:**

Solve for the transport of a chemical in the rectangular channel
$\Omega=(0,L)\times(0,H)$, with length $L=1.0$ m and height $H=0.10$ m.
The boundary partition is $\Gamma_{\mathrm{in}}=\{0\}\times[0,H]$,
$\Gamma_{\mathrm{out}}=\{L\}\times[0,H]$, and
$\Gamma_w=(0,L)\times\{0\}\cup(0,L)\times\{H\}$.

**Mesh:**

Use $440\times44$ uniform rectangular subdivisions, with each rectangle split
by a crossed, reflection-symmetric triangulation. The streamwise spacing is

$$h_x=\frac{L}{440}=2.272727\times10^{-3}\ \mathrm{m}.$$

**Model:**

Solve the steady advection--diffusion equation

$$u\cdot\nabla c-D\nabla^2c=0$$

for the scalar concentration $c$, with prescribed Poiseuille velocity

$$u(x,y)=(u_x(y),0),\qquad
u_x(y)=U_{\max}\frac{4y(H-y)}{H^2}.$$

**Boundary conditions:**

- Inlet $\Gamma_{\mathrm{in}}$: $c=0$.
- Outlet $\Gamma_{\mathrm{out}}$: $c=1$.
- Walls $\Gamma_w$: zero diffusive normal flux,
  $D\nabla c\cdot n=0$.

**Parameters:**

$$L=1.0\ \mathrm{m},\qquad H=0.10\ \mathrm{m},$$

$$U_{\max}=1.0\times10^{-2}\ \mathrm{m\,s^{-1}},\qquad
D=1.0\times10^{-5}\ \mathrm{m^2\,s^{-1}}.$$

The maximum-velocity 1D outlet-layer scale is

$$\delta=\frac{D}{U_{\max}}=1.0\times10^{-3}\ \mathrm{m}
=0.1\%\,L.$$

For the constant-velocity 1D guide

$$c_{1D}(x)=
\frac{\exp[-Pe_L(1-x/L)]-\exp(-Pe_L)}{1-\exp(-Pe_L)},
\qquad Pe_L=\frac{U_{\max}L}{D}=1000,$$

the interval from $c=0.01$ to the outlet value $c=1$ has width
$\ln(100)\delta=4.60517\times10^{-3}$ m, corresponding to approximately
$2.03$ streamwise mesh intervals. This 1D expression is a centerline guide,
not an exact solution of the 2D Poiseuille problem.

The maximum streamwise cell Peclet number is

$$Pe_h=\frac{U_{\max}h_x}{2D}=1.13636,$$

so include streamline-upwind/Petrov--Galerkin (SUPG) stabilization.

**Output and checks:**

- Save the concentration field in XDMF and checkpoint-readable form.
- Check the concentration bounds $0\le c\le1$.
- Check symmetry about $y=H/2$ and downstream monotonicity of the
  cross-sectional mean.
- Compare the centreline profile with the constant-$U_{\max}$ 1D guide above,
  while treating that comparison as diagnostic rather than exact validation.
