#Nina Martinez Diers

from math import sqrt
from numpy import array,arange
import matplotlib.pyplot as plt

#latex fonts
plt.rc('text',usetex=True)
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

#constants and initial vectors
G = 6.67408e-11 #gravitational constant with units m3 kg-1 s-2
M = 1.989e30 #mass of the sun w units kg

#fixed-step runge-kutta, adapted from odesim.py and rk4.py algorithms
#consulted pendulum_adaptive.py to help fix bugs
def f(r):
    x = r[0]
    y = r[1]
    R = sqrt(x**2 + y**2)
    vx = r[2]
    vy = r[3]
    #set of 1st order equations
    fx = vx
    fy = vy
    fvx = -G*M*x/R**3
    fvy = -G*M*y/R**3
    return array([fx,fy,fvx,fvy],float)

delta = 1000/1.63e9 #1km/yr w units m/s

a = 0.0
b = 1.63e9
N = 250000
h = (b-a)/N #we know this is a good step-size for h because it worked for the runge-kutta without the adaptive timestep
rho = 1 #because we know h is ok, we don't want to update it for the first time-step
t = a

tpoints = []
x_points = [] #in degrees, for plot in 8.4(a)
y_points = []

#initial conditions
initial_x = 4e12 #4 billion kilometers in meters
initial_y = 0 #m
initial_vx = 0 #m/s
initial_vy = 500 #m/s

#vectorize the differential equations
r = [initial_x,initial_y,initial_vx,initial_vy]

#initial values
tpoints.append(t)
x_points.append(r[0])
y_points.append(r[1])

while t < b:
    #calculate points for t+2h in 1 timestep
    k1 = 2*h*f(r)
    k2 = 2*h*f(r+0.5*k1)
    k3 = 2*h*f(r+0.5*k2)
    k4 = 2*h*f(r+k3)
    r1 = r + (k1+2*k2+2*k3+k4)/6
    #calculate points for t+2h in 2 timesteps
    r2 = r
    for i in range(2):
        k1 = h*f(r2)
        k2 = h*f(r2+0.5*k1)
        k3 = h*f(r2+0.5*k2)
        k4 = h*f(r2+k3)
        r2 += (k1+2*k2+2*k3+k4)/6
    #calculate rho for two estimates of points for t+2h
    diff = abs(r2[0]-r1[0])
    if diff == 0:
        rho = 1.1
    else:
        rho = 30*h*delta / abs(diff)
        
    if rho >= 1.0:
        t+=2*h
        h*=min(rho**0.25,1+10*delta) #insurance to make sure not to scale h too large
        r=r2
        #update lists with most current values
        tpoints.append(t)
        x_points.append(r[0])
        y_points.append(r[1])
    else:
        h*=rho**0.25

#Make plot
    
plt.plot(x_points,y_points,'k')
plt.xlabel(r'$x(t)$', fontsize = 24)
plt.ylabel(r'$y(t)$', fontsize = 24)
plt.show()
plt.clf()

#Make plot of time-points
plt.plot(x_points,y_points,'k.', markersize =0.1, alpha=0.5)
plt.xlabel(r'$x(t)$', fontsize = 24)
plt.ylabel(r'$y(t)$', fontsize = 24)
plt.show()