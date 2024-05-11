#Nina Martinez Diers
#2D Arterial fluid dynamics adapted using 12 steps to Navier–Stokes
#Implementing Vessel Flow

'''To adapt this to a model of an arterial segment, we break the periodic boundary condition and use an initial condition of laminar flow, with constant x-velocity $u$ running through the artery ($v=0$) and a constant pressure through the artery. For our boundary conditions, we fix the velocity $(u,v=0)$ entering the artery and the pressure at the exit of the artery. Using pseudotime in the pressure equation makes sure our pressure field evolves approriately.

bernoulli's eq: $$P_{in} = P_{out} + \rho(\vec{v_{out}}^2/2 - \vec{v_{in}}^2/2 + g(h_{out}-h_{in}))$$

Initial conditions:

hout = hin #horizontal artery

flow rate in: 0.8->4 cubic cm/s. We will use (0.8->4)^2/3 for sq cm/s flow for planar approximation of vessel.
start with v = 0 at in for laminar flow, u = one of above quantities.

Arterial pressure usually ranges from 9.332-17.332 kPa, so values within that range is what we will use for BCs

Viscosity ν = 0.003528 Pa.s and density ρ = 1060 kg.m-3.'''
import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D

def build_up_b(rho, dt, dx, dy, u, v):
    '''The function `build_up_b` below represents the contents of the square brackets 
    in the pressure Poisson equation (PPE), so that the entirety of the PPE is slightly more manageable. '''

    b = numpy.zeros_like(u)
    b[1:-1, 1:-1] = (rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) +
                                      (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -
                            ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 -
                            2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *
                                 (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx))-
                            ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2))
    return b

def vessel_pressure_poisson(p, dx, dy, b):
    '''The function `vessel_pressure_poisson` is also defined to help segregate 
    the different rounds of calculations.  The pseudo-time variable `nit` is a 
    sub-iteration in the Poisson calculation to help ensure a divergence-free field.  '''
    pn = numpy.empty_like(p)
    
    for q in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 + 
                          (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) /
                          (2 * (dx**2 + dy**2)) -
                          dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * 
                          b[1:-1,1:-1])

        
        p[:, 0] = p[:, 1]   # dp/dx = 0 at x = 0 #Pressure is not changing coming in. Change this for pulsatile
        
        # Wall boundary conditions (No-slip)
        p[-1, :] = p[-2, :]  # dp/dy = 0 at y = 2
        p[0, :] = p[1, :]  # dp/dy = 0 at y = 0
        
        #pressure BC at x = 2. Pressure is always going to be the same going out. Change this for pulsatile.
        p[:,-1] = P
        
    return p

def vessel_flow(nt, u, v, dt, dx, dy, p, rho, nu, F):
    '''`vessel_flow` finishes the calculation of the PDE system.'''
    un = numpy.empty_like(u)
    vn = numpy.empty_like(v)
    b = numpy.zeros((ny, nx))
    
    for n in range(nt):
    
        un = u.copy()
        vn = v.copy()

        b = build_up_b(rho, dt, dx, dy, u, v)
        p = vessel_pressure_poisson(p, dx, dy, b)

        u[1:-1, 1:-1] = (un[1:-1, 1:-1] -
                         un[1:-1, 1:-1] * dt / dx * 
                        (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                         vn[1:-1, 1:-1] * dt / dy * 
                        (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                         dt / (2 * rho * dx) * 
                        (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                         nu * (dt / dx**2 * 
                        (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                         dt / dy**2 * 
                        (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])) + 
                         F * dt)

        v[1:-1, 1:-1] = (vn[1:-1, 1:-1] -
                         un[1:-1, 1:-1] * dt / dx * 
                        (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                         vn[1:-1, 1:-1] * dt / dy * 
                        (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                         dt / (2 * rho * dy) * 
                        (p[2:, 1:-1] - p[0:-2, 1:-1]) +
                         nu * (dt / dx**2 *
                        (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                         dt / dy**2 * 
                        (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))
        
        #BC: constant u @ x = 0 (except at walls)
        u[1:-2, 0] = a
        
        #BC: v=0 @ x = 0
        v[:,0] = 0
        
        
       # Wall BC: u,v = 0 @ y = 0,2
        u[0, :] = 0
        u[-1, :] = 0
        v[0, :] = 0
        v[-1, :]=0
        
    return u, v, p

#Huge artery size (unrealistic), slow movement, high viscous (out of normal range), normal density
L = 2
w = 1
dt = 0.00001

nx = 200
ny = 100
nit = 50
dx = L / (nx - 1)
dy = w / (ny - 1)
x = numpy.linspace(0, L, nx)
y = numpy.linspace(0, w, ny)
X, Y = numpy.meshgrid(x, y)

rho = 1060
nu = 1
a = 8.6e-5 #u_in
P = 10

#arbitrary
F = 1
#change F = 10 to see very high turbulence

Re = rho*a*w/nu
print('Reynolds number = ', Re)

for nt in [100,500,1000,1500,2000]:
    #Call function
    #initial condition is laminar flow in x direction
    u = a* numpy.ones((ny, nx))#x-velocity is some constant (starting with all ones)
    v = numpy.zeros((ny, nx))#y-velocity is zero
    p = P* numpy.ones((ny, nx))#for constant velocity and a horizontal artery, pressure is zero everywhere 
    b = numpy.zeros((ny, nx))#this is an empty matrix that will be populated with correct values in build-up-b function
    u, v, p = vessel_flow(nt, u, v, dt, dx, dy, p, rho, nu, F)
    
    #Plot results in velocity field graph with pressure field as contour:
    fig = pyplot.figure(figsize=(13,4), dpi=100)
    # plotting the pressure field as a contour
    pyplot.subplot(1, 2, 1)
    pyplot.contourf(X, Y, p, alpha=0.5, cmap=cm.viridis)  
    pyplot.colorbar()
    # plotting the pressure field outlines
    pyplot.contour(X, Y, p, cmap=cm.viridis)  
    # plotting velocity field
    pyplot.quiver(X[::7, ::7], Y[::7, ::7], u[::7, ::7], v[::7, ::7]) #only shows every 7th point so we can see the data better
    pyplot.xlabel('X (m)')
    pyplot.ylabel('Y (m)')
    
    #Plot results in streamline field graph also with pressure field as contour:
    pyplot.subplot(1, 2, 2)
    pyplot.contourf(X, Y, p, alpha=0.5, cmap=cm.viridis)
    pyplot.colorbar(label='Pressure (Pa)')
    pyplot.contour(X, Y, p, cmap=cm.viridis)
    pyplot.streamplot(X, Y, u, v)
    pyplot.xlabel('X (m)')