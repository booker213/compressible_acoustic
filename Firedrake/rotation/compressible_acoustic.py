from firedrake import *

# Script to solve compressible acoustic waves in upto three dimensions
# c^2_0 and \rho_0 will be considered to be scaled to be equal to 1.

# Test problem will be the one-dimensional waves considered in the MATLAB scripts.

# Create mesh
# Current mesh is a unit line  with Nx elements.

Nx = 16
Ny = 16
Nz = 16
mesh = UnitCubeMesh(Nx,Ny, Nz)

# Declare timestep
dt = 1./16.

# Declare initial and end time
# Period of waves considered in test is 1s
# We will consider 3 periods initially
t = 0.
end = 1.

# Declare order of the basis in the elements
# Test problem will consider order 0 - a finite volume scheme
order_basis = 0 

# Declare flux indicator function
theta = Constant(0.5)


# Declare function spaces on the mesh
# and create mixed function space
# for coupled Poisson Bracket.

V = VectorFunctionSpace(mesh, "DG", order_basis)
R = FunctionSpace(mesh, "DG", order_basis)

W = V*R

# Define trial and test functions on the mixed space

(u, rho) = TrialFunctions(W)
(dFdu_vec, dFdrho) = TestFunctions(W)


# Initial conditions

# Function space
w0 = Function(W)
(u0,rho0) = split(w0)

# Interpolate expressions
u0,rho0 = w0.split()

u0.interpolate(Expression(["sin(2*pi*x[0])*sin(2*pi*0.125)", "sin(2*pi*x[1])*sin(2*pi*0.125)","sin(2*pi*x[2])*sin(2*pi*0.125)"] ))

rho0.interpolate(Expression("cos(2*pi*x[0])*cos(2*pi*0.125) + cos(2*pi*x[1])*cos(2*pi*0.125)+ cos(2*pi*x[2])*cos(2*pi*0.125) "))

# Assemble initial energy
E0 = assemble ( (0.5*(inner(u0,u0) + rho0**2))*dx)

#print E0

# Declare angular velocity
Omega = Function(V)
Omega.interpolate(Expression(["0.0", "0.0", "0.0"]))


# Set up internal boundary normal on mesh 
# Needed for numerical fluxes in bilinear form

n = FacetNormal(mesh)


# Create bilinear form for linear solver
# Bilinear problem is of the form 
# a(u,v) = L(v)
# Using our Poisson bracket and an implicit midpoint rule
# we see that
# a(u,v) = u^{n+1}*dFdu_vec + rho^{n+1)*dFdrho - 0.5*dt*PB(u^{n+1}, rho^{n+1})
# L(v) = u^{n}*dFdu_vec + rho^{n)*dFdrho + 0.5*dt*PB(u^{n}, rho^{n})


# We note that there are no boundary surface integrals ds, as we require
# the normal of the variational derivative and test function to vanish 
# at the boundary.

#Define varitional derivatives
dHdu = u
dHdrho = rho

dHdu0 = u0
dHdrho0 = rho0

a0 = (dot(u, dFdu_vec) + rho*dFdrho )*dx
a1 = (dot(-grad(dHdrho), dFdu_vec) + dot(dHdu, grad(dFdrho)) - dot(cross(2*Omega, dHdu),dFdu_vec))*dx
a2 = (jump( dFdrho)*dot((dHdu('-')*(1-theta)+dHdu('+')*theta), n('-')))*dS
a3 = (-jump( dHdrho)*dot((dFdu_vec('-')*(1-theta)+dFdu_vec('+')*theta), n('-')))*dS

a = a0 - 0.5*dt*(a1+a2+a3)


L0 = (dot(u0, dFdu_vec) + rho0*dFdrho )*dx
L1 = (dot(-grad(dHdrho0), dFdu_vec) + dot(dHdu0, grad(dFdrho))- dot(cross(2*Omega, dHdu0),dFdu_vec))*dx
L2 = (jump( dFdrho)*dot((dHdu0('-')*(1-theta)+dHdu0('+')*theta), n('-')))*dS
L3 = (-jump( dHdrho0)*dot((dFdu_vec('-')*(1-theta)+dFdu_vec('+')*theta), n('-')))*dS

L = L0 + 0.5*dt*(L1+L2+L3)

# Storage for visualisation
outfile = File('./Results/compressible_acoustic_results.pvd')

u0.rename("Velocity")
rho0.rename("Density")


# Output initial conditions
outfile.write(u0,rho0, time = t)


out = Function(W)
# File for energy output
E_file = open('./Results/energy.txt', 'w')


# Solve loop

while (t < end):
 # Update time
 t+= dt
 
 solve(a == L, out)
 u, rho = out.split()

 # Assign appropriate name in results file
 u.rename("Velocity")
 rho.rename("Density")
 
 # Output results
 outfile.write(u, rho, time =t)
 
 # Assign output as previous timestep for next time update
 u0.assign(u)
 rho0.assign(rho)

 # Assemble initial energy
 E = assemble ( (0.5*(inner(u,u) + rho**2))*dx)
 
 E_file.write('%-10s %-10s\n' % (t,abs((E-E0)/E0)))
 # Print time and energy drift, drift should be around machine precision.
 print t, abs((E-E0)/E0)

# Create analytic solutions for error analysis
exact_rho= Function(R)
exact_rho.interpolate(Expression("cos(2*pi*x[0])*cos(2*pi*(t+0.125)) + cos(2*pi*x[1])*cos(2*pi*(t+0.125))+ cos(2*pi*x[2])*cos(2*pi*(t+0.125)) ", t = t))


# Print error for density
error_rho = errornorm(rho, exact_rho,  norm_type='L2')
print error_rho


# Close energy write
E_file.close()
