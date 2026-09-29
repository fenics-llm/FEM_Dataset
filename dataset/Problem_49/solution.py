# Extracted from buckling_3d_solid.ipynb.
# Markdown cells are retained as comments; notebook-only magics are disabled.

# # Linear Buckling Analysis of a 3D solid
#
# This demo has been written in collaboration with [Eder Medina (Harvard University)](mailto:e_medina@g.harvard.edu).
#
# In this numerical tour,  we will demonstrate how to compute the buckling modes of a three-dimensional elastic solid under a conservative loading.
# The critical buckling loads and buckling modes are computed by solving a generalized non-hermitian eigenvalue problem using the `SLEPcEigensolver`. This example is closely related to the [Modal analysis of an elastic structure](https://comet-fenics.readthedocs.io/en/latest/demo/modal_analysis_dynamics/cantilever_modal.py.html) and the [Eulerian buckling of a beam](../beam_buckling/beam_buckling.ipynb) demos. For more details on the theoretical formulation, the reader can refer to [[NGU00]](#References).
#
# $$\newcommand{\bq}{\boldsymbol{q}}
# \newcommand{\bu}{\boldsymbol{u}}
# \newcommand{\bv}{\boldsymbol{v}}
# \newcommand{\bt}{\boldsymbol{t}}
# \newcommand{\be}{\boldsymbol{e}}
# \newcommand{\beps}{\boldsymbol{\varepsilon}}
# \newcommand{\bsig}{\boldsymbol{\sigma}}
# \newcommand{\T}{^{\text{T}}}$$
# ## Geometrically Nonlinear Elasticity
#
# To formulate the buckling problem, we must consider a general geometrically nonlinear framework. More precisely, we will consider the case of small-strain/large displacements that are expected in typical buckling analysis.
#
# We consider the Green-Lagrange strain tensor $\be$ defined as:
#
# $$e_{ij}  = \frac{1}{2} (u_{i,j} + u_{j,i}) + \frac{1}{2} u_{k,i}u_{k,j}$$
#
# Note that this can be rewritten in tensor format as
#
# $$ \be = \beps(\bu) + \dfrac{1}{2}\bq(\bu, \bu) $$
#
# where we isolate the linear $\beps(\bu)$ and nonlinear quadratic part $\bq(\bu, \bu)=\nabla\bu\T\nabla\bu$.
#
# Since we restrict to a small strain assumption, we adopt a standard linear Hookean behaviour for the material material model which can be described by the following quadratic elastic energy density:
#
# $$ \psi(\be) = \dfrac{1}{2}\be:\mathbb{C}: \be $$
#
# where $\mathbb{C}$ is the elastic moduli tensor.
#
# We seek to find the underlying equilibrium displacement field $\bu$ as a function of a scalar loading parameter $\lambda$ under a set of prescribed external surface forces $\bt(\lambda)$.
#
# ## Equilibrium and stability conditions
# Let us consider the total potential energy functional:
#
# $$ \mathcal{E}_\text{pot}(\bu;\lambda) = \int_\Omega \psi(\be(\bu)) d\Omega - \int_{\Gamma_N} \bt(\lambda) \cdot \bu dS $$
#
# The equilibrium equations are obtained from the first variation which must vanish:
#
# $$ \partial_u\mathcal{E}_\text{pot}[\delta\bu] = \int_\Omega \be(\bu):\mathbb{C}:\delta\be[\delta\bu] d\Omega - \int_{\Gamma_N} \bt(\lambda) \cdot \delta\bu dS $$
#
# where $\delta\bu$ is any kinematically admissible perturbation direction and where:
#
# $$
# \delta\be[\delta\bu] = \beps(\delta\bu) + \nabla\bu\T\nabla\delta\bu = \beps(\delta\bu) + \bq(\bu,\delta\bu)
# $$
#
# The stability conditions of an equilibrium are obtained from the Lejeune-Dirichlet theorem stating that the energy second variation must be positive. Here the second variation bilinear form is given by:
#
# $$ \partial_{uu}\mathcal{E}_\text{pot}[\delta\bu,\delta \bv] = \int_\Omega (\delta \be[\delta\bv]:\mathbb{C}:\delta\be[\delta\bu] + \be:\mathbb{C}:\delta^2\be[\delta\bu,\delta\bv] )d\Omega $$
#
# where:
# $$
# \delta^2\be[\delta\bu,\delta\bv] = \delta\bq(\bu,\delta \bu)[\delta\bv] = \nabla \delta\bv\T\nabla\delta\bu = \bq(\delta\bv,\delta\bu)
# $$
#
# ## Bifurcation point
#
# A given equilibrium solution $\lambda^c, \bu^c=\bu(\lambda^c)$ is a bifurcation point if the second variation vanishes for some $\delta\bu$. $\lambda^c$ defines the bifurcation load and $\delta\bu$ the bifurcation mode. Hence, we have:
# $$
# \int_\Omega ((\beps(\delta\bv) + \bq(\bu^c,\delta\bv)):\mathbb{C}:(\beps(\delta\bu) + \bq(\bu^c,\delta\bu)) + \be(\bu^c):\mathbb{C}:\bq(\delta\bv,\delta\bu))d\Omega = 0
# $$
#
# ## Linear buckling analysis
#
# When performing a linear buckling analysis, we look for bifurcation points on the fundamental equilibrium branch obtained by a small displacement assumption. As a result, in the above bifurcation equation the quadratic contribution in the first term can be neglected and the Green-Lagrange strain $\be(\bu^c)$ can be replaced by the linearized strain $\beps(\bu^c)$. Besides, in the small-displacement assumption, if the loading depends linearly on the load parameter $\lambda$, the small displacement solution also depends linearly on it. We therefore have $\bu^c= \lambda^c\bu_0$ and the bifurcation condition becomes:
# $$
# \int_\Omega (\beps(\delta\bv):\mathbb{C}:\beps(\delta\bu)  + \lambda^c\beps(\bu_0):\mathbb{C}:\bq(\delta\bv,\delta\bu))d\Omega = 0
# $$
#
# After, introducing the inital pre-stress $\bsig_0 = \mathbb{C}:\beps(\bu_0)$, the last term can be rewritten as $\bsig_0:\bq(\delta\bv,\delta\bu) = (\nabla \delta \bv) \bsig_0 (\nabla \delta\bu)\T$ so that :
# $$
# \int_\Omega (\beps(\delta\bv):\mathbb{C}:\beps(\delta\bu)  + \lambda^c(\nabla \delta \bv) \bsig_0 (\nabla \delta\bu)\T)d\Omega = 0
# $$
#
# We recognize here a linear eigenvalue problem where the first term corresponds to the classical linear elasticity bilinear form and the second term is an additional contribution depending on the pre-stressed state $\bu_0$. Transforming this variational problem into matrix notations yields:
# $$
# [\mathbf{K}]+\lambda[\mathbf{K}_G(\bu_0)]=0
# $$
#
# where $[\mathbf{K}]$ is the classical linear elastic stiffness matrix and $[\mathbf{K}_G(\bu_0)]$ the so-called *geometrical stiffness* matrix.
# ## FEniCS implementation
#
# ### Defining the Geometry
#
# Here we will consider a fully three dimensional beam-like structure with length $L = 1$ and of rectangular cross-section with height $h = 0.01$, and width $b = 0.03$.
from dolfin import *
import json
import os
import numpy as np

problem_dir = os.path.dirname(os.path.abspath(__file__))

L, b, h = 1, 0.01, 0.03
Nx, Ny, Nz = 51, 5, 5

mesh = BoxMesh(Point(0, 0, 0), Point(L, b, h), Nx, Ny, Nz)

# ### Solving for the pre-stressed state $\bu_0$
#
# We first need to compute the pre-stressed linear state $\bu_0$. It is obtained by solving a simple linear elasticity problem. The material is assumed here isotropic and we consider a pre-stressed state obtained from the application of a unit compression applied at the beam right end in the $X$ direction while the $Y$ and $Z$ displacement are fixed, mimicking a simple support condition. The beam is fully clamped on its left end.
E , nu = 1e3, 0.

mu = Constant(E/2/(1+nu))
lmbda = Constant(E*nu/(1+nu)/(1-2*nu))

def eps(v):
    return sym(grad(v))
def sigma(v):
    return lmbda*tr(eps(v))*Identity(v.geometric_dimension())+2*mu*eps(v)

# Compute the linearized unit preload
N0 = 1
T = Constant((-N0, 0, 0))
V = VectorFunctionSpace(mesh, 'Lagrange', degree = 2)
v = TestFunction(V)
du = TrialFunction(V)

# Two different ways to define how to impose boundary conditions
def left(x,on_boundary):
    return near(x[0], 0.)

def right(x, on_boundary):
    return near(x[0], L)

boundary_markers = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)
AutoSubDomain(right).mark(boundary_markers, 1)

# Clamped Boundary and Simply Supported Boundary Conditions
bcs = [DirichletBC(V, Constant((0, 0, 0)), left),
       DirichletBC(V.sub(1), Constant(0), right),
       DirichletBC(V.sub(2), Constant(0), right)]

# Linear Elasticity Bilinear Form
a = inner(sigma(du), eps(v))*dx

# External loading on the right end
l = dot(T, v)*ds(1)

# We will reuse our stiffness matrix in the linear buckling computation and hence we supply a `PETScMatrix` and `PETScVector` to the system assembler. We then use the `solve` function on these objects and specify the linear solver method.

K = PETScMatrix()
f = PETScVector()
assemble_system(a, l, bcs, A_tensor=K, b_tensor=f)

u = Function(V,name = "Displacement")
solve(K, u.vector(), f, "mumps")

# Output the trivial solution
artifact_dir = os.path.join(problem_dir, "reference solution artifacts")
os.makedirs(artifact_dir, exist_ok=True)
ffile = XDMFFile(os.path.join(artifact_dir, "buckling_results.xdmf"))
ffile.parameters["functions_share_mesh"] = True
ffile.parameters["flush_output"] = True
ffile.write(u, 0)

# ### Forming the geometric stiffness matrix and the buckling eigenvalue problem
#
# From the previous solution, the prestressed state $\bsig_0$ entering the geometric stiffness matrix expression is simply given by `sigma(u)`. After having formed the geometric stiffness matrix, we will call the `SLEPcEigenSolver` for solving a generalized eigenvalue problem $\mathbf{Ax}=\lambda \mathbf{Bx}$ with here $\mathbf{A}=\mathbf{K}$ and $\mathbf{B}=-\mathbf{K}_G$. We therefore include directly the negative sign in the definition of the geometric stiffness form:
kgform = -inner(sigma(u), grad(du).T*grad(v))*dx

KG = PETScMatrix()
assemble(kgform, KG)

# Zero out the rows to enforce boundary condition
for bc in bcs:
    bc.zero(KG)

# The buckling points of interest are the smallest critical loads, i.e. the points corresponding to the smallest $\lambda^c$.
#
# We now pass some parameters to the SLEPc solver. In particular, we peform a shift-invert transform as discussed in the [Eulerian buckling of an elastic beam](https://comet-fenics.readthedocs.io/en/latest/demo/beam_buckling/beam_buckling.html) tour.
# What Solver Configurations exist?
eigensolver = SLEPcEigenSolver(K, KG)
print(eigensolver.parameters.str(True))

eigensolver.parameters['problem_type'] = 'gen_hermitian'
eigensolver.parameters['solver'] = "krylov-schur"
eigensolver.parameters['tolerance'] = 1.e-12
eigensolver.parameters['spectral_transform']  = 'shift-and-invert'
eigensolver.parameters['spectral_shift'] = 1e-2

# Request the smallest 3 eigenvalues
N_eig = 3
print("Computing {} first eigenvalues...".format(N_eig))
eigensolver.solve(N_eig)
print("Number of converged eigenvalues:", eigensolver.get_number_converged())

eigenmode = Function(V)

# Analytical beam theory solution
# F_cr,n = alpha_n^2*EI/L^2
# where alpha_n are solutions to tan(alpha) = alpha
alpha = np.array([1.4303, 2.4509, 3.4709] + [(n+1/2) for n in range(4, 10)])*pi
I = b**3*h/12
S = b*h

critical_loads = []
beam_loads = []
buckling_modes = []
for i in range(eigensolver.get_number_converged()):
    # Extract eigenpair
    r, c, rx, cx = eigensolver.get_eigenpair(i)

    critical_load_an = alpha[i]**2*float(E*I/N0/S)/L**2
    print("Critical load factor {}: {:.5f} FE | {:.5f} Beam".format(i+1, r, critical_load_an))
    critical_loads.append(float(r))
    beam_loads.append(float(critical_load_an))

    # Initialize function and assign eigenvector
    eigenmode.vector()[:] = rx
    eigenmode.rename("Eigenmode {}".format(i+1), "")

    ffile.write(eigenmode, 0)
    mode = Function(V)
    mode.assign(eigenmode)
    mode.rename("phi_{}".format(i+1), "buckling mode {}".format(i+1))
    buckling_modes.append(mode)

ffile.close()
with XDMFFile(mesh.mpi_comm(), os.path.join(problem_dir, "mesh.xdmf")) as out:
    out.write(mesh)
with XDMFFile(mesh.mpi_comm(), os.path.join(problem_dir, "solution.xdmf")) as out:
    out.parameters["functions_share_mesh"] = True
    for mode in buckling_modes:
        out.write(mode, 0.0)
with XDMFFile(mesh.mpi_comm(), os.path.join(problem_dir, "read_checkpoint.xdmf")) as out:
    for i, mode in enumerate(buckling_modes, start=1):
        out.write_checkpoint(mode, "phi_{}".format(i), 0, XDMFFile.Encoding.HDF5, i > 1)
with open(os.path.join(problem_dir, "read_checkpoint.json"), "w") as out:
    json.dump({"fields":[
        {"name":"phi_{}".format(i),"family":"CG","degree":2,
         "value_shape":[3],"mode":i,"critical_load_factor":critical_loads[i-1]}
        for i in range(1,len(buckling_modes)+1)],"critical_loads_fe":critical_loads,
               "critical_loads_beam":beam_loads},out,indent=2)

# > Note that the analysis above is not limited to problems subject to Neumann boundary conditions.
# One can equivalenty compute the prestress resulting from a displacement control simulation and the critical buckling load will correspond to the critical buckling displacement of the structure.
# ## References
#
# [NGU00]: Nguyen, Q. S. (2000). *Stability and nonlinear solid mechanics*. Wiley.
