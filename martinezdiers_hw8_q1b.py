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



a = 0.0
b = 165*365*24*3600
N = 250000
h = (b-a)/N

tpoints = arange(a,b,h)
x_points = [] #in degrees, for plot in 8.4(a)
y_points = []

#initial conditions
initial_x = 4e12 #4 billion kilometers in meters
initial_y = 0 #m
initial_vx = 0 #m/s
initial_vy = 500 #m/s

#vectorize the differential equations
r = [initial_x,initial_y,initial_vx,initial_vy]

for t in tpoints:
    #solve x and y equations:
    x_points.append(r[0])
    y_points.append(r[1])
    k1 = h*f(r)
    k2 = h*f(r+0.5*k1)
    k3 = h*f(r+0.5*k2)
    k4 = h*f(r+k3)
    r += (k1+2*k2+2*k3+k4)/6

print(max(x_points))
#Make plot
    
plt.plot(x_points,y_points,'k')
plt.xlabel(r'$x(t)$', fontsize = 24)
plt.ylabel(r'$y(t)$', fontsize = 24)
plt.show()
plt.clf()
