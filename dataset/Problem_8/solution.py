# Copied from Results_ALL-FEM-main/reference solutions/solid/8/solid8.py for dataset Problem_8.
# Original benchmark: Solid Mechanics Problem 8 (Medium).

from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Geometry & mesh: Omega = (0,1.0) x (0,0.20), structured 50 x 25
# ------------------------------------------------------------
Lx, Ly = 1.0, 0.20
nx, ny = 50, 25
mesh = RectangleMesh(Point(0.0, 0.0), Point(Lx, Ly), nx, ny)

# ------------------------------------------------------------
# Function spaces
# ------------------------------------------------------------
V = VectorFunctionSpace(mesh, "CG", 1)        # displacement
Tspace = TensorFunctionSpace(mesh, "CG", 1)   # stress tensor
Sspace = FunctionSpace(mesh, "CG", 1)         # scalar fields (e.g., von Mises)

# ------------------------------------------------------------
# Material data (orthotropic lamina, local axes 1-2), plane-stress
# ------------------------------------------------------------
E1   = 40e9     # Pa
E2   = 10e9     # Pa
G12  = 5e9      # Pa
nu12 = 0.25
nu21 = nu12*E2/E1

# Plane-stress reduced stiffness in local 1-2 (Voigt with engineering shear gamma_12)
den = 1.0 - nu12*nu21
Q11 = E1/den
Q22 = E2/den
Q12 = nu12*E2/den
Q66 = G12

# Rotation (local 1-2 -> global x-y) by theta (anti-clockwise)
theta_deg = 30.0
theta = np.deg2rad(theta_deg)
c = np.cos(theta)
s = np.sin(theta)

# Transformed reduced stiffness Qbar (engineering shear form)
Q11b = Q11*c**4 + 2.0*(Q12 + 2.0*Q66)*s**2*c**2 + Q22*s**4
Q22b = Q11*s**4 + 2.0*(Q12 + 2.0*Q66)*s**2*c**2 + Q22*c**4
Q12b = (Q11 + Q22 - 4.0*Q66)*s**2*c**2 + Q12*(s**4 + c**4)
Q16b = (Q11 - Q12 - 2.0*Q66)*c**3*s - (Q22 - Q12 - 2.0*Q66)*c*s**3
Q26b = (Q11 - Q12 - 2.0*Q66)*c*s**3 - (Q22 - Q12 - 2.0*Q66)*c**3*s
Q66b = (Q11 + Q22 - 2.0*Q12 - 2.0*Q66)*s**2*c**2 + Q66*(s**4 + c**4)

# Assemble 3x3 stiffness in global x-y for [exx, eyy, gamma_xy]^T -> [sxx, syy, txy]^T
D_np = np.array([[Q11b, Q12b, Q16b],
                 [Q12b, Q22b, Q26b],
                 [Q16b, Q26b, Q66b]], dtype=float)

# Wrap as Constant matrix for UFL
D = as_matrix(((D_np[0,0], D_np[0,1], D_np[0,2]),
               (D_np[1,0], D_np[1,1], D_np[1,2]),
               (D_np[2,0], D_np[2,1], D_np[2,2])))

# ------------------------------------------------------------
# Kinematics & constitutive law (plane-stress, engineering shear)
# ------------------------------------------------------------
u = TrialFunction(V)
v = TestFunction(V)

def strain_vec(w):
    eps = sym(grad(w))
    exx = eps[0,0]
    eyy = eps[1,1]
    gxy = 2.0*eps[0,1]  # engineering shear gamma_xy
    return as_vector((exx, eyy, gxy))

def stress_vec(w):
    return D*strain_vec(w)  # [sxx, syy, txy]

def stress_tensor_from_vec(sv):
    # Convert [sxx, syy, txy] to 2x2 tensor
    return as_tensor(((sv[0], sv[2]),
                      (sv[2], sv[1])))

# Bilinear and linear forms
a = inner(stress_vec(u), strain_vec(v))*dx

# ------------------------------------------------------------
# Boundary conditions
# ------------------------------------------------------------
# Bottom edge: y = 0 fixed (ux=0, uy=0)
class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.0, DOLFIN_EPS)

bottom = Bottom()
bc_bot = DirichletBC(V, Constant((0.0, 0.0)), bottom)

# Neumann traction on top edge: t = (0, -10 MPa)
tmag = -10e6  # Pa (N/m^2); in 2D (unit thickness) this acts as N/m
t = Constant((0.0, tmag))

# Mark boundaries to apply traction only on top
class Top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], Ly, DOLFIN_EPS)

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
Top().mark(boundaries, 1)
dsN = Measure("ds", domain=mesh, subdomain_data=boundaries)

L = dot(t, v)*dsN(1)  # only on top

# ------------------------------------------------------------
# Solve
# ------------------------------------------------------------
u_sol = Function(V, name="displacement")
solve(a == L, u_sol, bc_bot,
      solver_parameters={"linear_solver": "mumps"})

# ------------------------------------------------------------
# Post-processing: stress tensor and von Mises (plane-stress)
# ------------------------------------------------------------
sv = stress_vec(u_sol)
sigma = stress_tensor_from_vec(sv)                  # 2x2 tensor
sigma_fun = project(sigma, Tspace, solver_type="mumps")
sigma_fun.rename("stress", "stress")

# von Mises for plane-stress: sqrt(sxx^2 - sxx*syy + syy^2 + 3*txy^2)
sxx, syy, txy = sv[0], sv[1], sv[2]
vm_expr = sqrt(sxx**2 - sxx*syy + syy**2 + 3.0*txy**2)
vm_fun = project(vm_expr, Sspace, solver_type="mumps")
vm_fun.rename("von_mises", "von_mises")

# ------------------------------------------------------------
# Save XDMF outputs (u, sigma, von Mises)
# ------------------------------------------------------------
with XDMFFile(mesh.mpi_comm(), "q8_solution.xdmf") as xdmf:
    xdmf.parameters["functions_share_mesh"] = True
    xdmf.parameters["rewrite_function_mesh"] = False
    xdmf.write(u_sol, 0.0)
    xdmf.write(sigma_fun, 0.0)
    xdmf.write(vm_fun, 0.0)

# ------------------------------------------------------------
# Plots: colour maps of ux and von Mises
# ------------------------------------------------------------
plt.figure()
p0 = plot(u_sol.sub(0), title="Horizontal displacement, u_x", mode="color")
plt.colorbar(p0)
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.tight_layout()
plt.savefig("q8_ux.png", dpi=300)

plt.figure()
p1 = plot(vm_fun, title="von Mises stress", mode="color")
plt.colorbar(p1)
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.tight_layout()
plt.savefig("q8_vm.png", dpi=300)

print("Done. Files written: q8_solution.xdmf, q8_ux.png, q8_vm.png")

# -----------------------------------------------------------------------------
# Canonical dataset outputs
# -----------------------------------------------------------------------------
import json

mesh_file = XDMFFile(mesh.mpi_comm(), "mesh.xdmf")
mesh_file.write(mesh)
mesh_file.close()

u_sol.rename("u", "displacement field")

solution_file = XDMFFile(mesh.mpi_comm(), "solution.xdmf")
solution_file.parameters["flush_output"] = True
solution_file.parameters["functions_share_mesh"] = True
solution_file.write(u_sol, 0.0)
solution_file.close()

checkpoint_file = XDMFFile(mesh.mpi_comm(), "read_checkpoint.xdmf")
checkpoint_file.write_checkpoint(u_sol, "u", 0.0, XDMFFile.Encoding.HDF5, False)
checkpoint_file.close()

read_checkpoint_metadata = {
    "format": "FEniCS XDMFFile.write_checkpoint",
    "time": 0.0,
    "fields": [
        {
            "symbol": "u_sol",
            "checkpoint_name": "u",
            "meaning": "displacement field",
            "read_example": "XDMFFile(mesh.mpi_comm(), 'read_checkpoint.xdmf').read_checkpoint(function, 'u', -1)",
        },
    ],
}
with open("read_checkpoint.json", "w", encoding="utf-8") as metadata_file:
    json.dump(read_checkpoint_metadata, metadata_file, indent=2)
    metadata_file.write("\n")
