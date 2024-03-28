#Newman 8.4 and 8.5
#Nina Martinez Diers
#adapted from rk4.py, an implementation of the 4th order runge-kutta method
from math import sin, cos
from numpy import arange, array, pi
import matplotlib.pyplot as plt
from vpython import *

#latex fonts
plt.rc('text',usetex=True)
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

#constants among all problems
g = 9.81 #m/s^2
l = 0.10 #m

#Newman 8.4(a)
def f(r):
    theta = r[0]
    omega = r[1]
    ftheta = omega
    fomega = -g/l*sin(theta)
    return array([ftheta, fomega], float)

a = 0.0
b = 10.0
N = 100000
h = (b-a)/N

tpoints = arange(a,b,h)
theta_points = [] #in degrees, for plot in 8.4(a)
theta_radians = [] #in radians, for animation in 8.4(b)
initial_theta = 179 #in degrees
initial_omega = 0.0

#vectorize the differential equation, with theta in radians
r = [initial_theta*pi/180, initial_omega]

for t in tpoints:
    theta_points.append(r[0]*180/pi)
    theta_radians.append(r[0])
    k1 = h*f(r)
    k2 = h*f(r+0.5*k1)
    k3 = h*f(r+0.5*k2)
    k4 = h*f(r+k3)
    r += (k1+2*k2+2*k3+k4)/6

#Make plot
    
plt.plot(tpoints,theta_points,'k')
plt.xlabel(r'Time(s)', fontsize = 24)
plt.ylabel(r'$\theta(t)$', fontsize = 24)
plt.show()
plt.clf()

#Newman 8.4(b)

#8.4(b) --> Animations extra-credit
#adapted revolve.py from Computational Physics by Mark Newman

scene = canvas() #this will show the animation
#make pendulum out of sphere and cylinder
#initial position of mass at the end of the pendulum is at the origin
s = sphere(pos = vec(0,0,0), radius=0.01)

#the hinge of the cylinder is at position (0,l,0) with an axis of (0,-1,0). 
#This makes a pendulum of length l that starts hanging straight down, terminating at the origin.
c = cylinder(pos = vec(0,l,0), axis = vec(0,-1,0), radius=0.002)
for theta in theta_radians: #need theta vals in radians for trig in calculations of x and y
    rate(3000) #this rate makes the program run at a rate that is comfortably observed
    x = l*sin(theta)
    y = l*(1-cos(theta))
    s.pos = vec(x,y,0)
    c.axis = vec(x,y-l,0) #this axis makes the end of the cylinder terminate 
                          #at the location of the end of the pendulum


#8.5(a) - the driven pendulum
C = 2 #s^-2
Omega = 5 #s^-2

def f_driven(r,t):
    theta = r[0]
    omega = r[1]
    ftheta = omega
    fomega = -g/l*sin(theta) + C*cos(theta)*sin(Omega*t)
    return array([ftheta, fomega], float)

a = 0.0
b = 100.0
N = 1000000
h = (b-a)/N

tpoints = arange(a,b,h)
theta_points = []
initial_theta = 0.0 #in degrees. start with the pendulum handing straight down.
initial_omega = 0.0 #m/s. start at rest
#convert theta to radians
initial_theta = initial_theta *pi/180

#vectorize the differential equation
r = [initial_theta, initial_omega]

for t in tpoints:
    theta_points.append(r[0]*180/pi) #convert theta back to degrees and store value in list
    k1 = h*f_driven(r,t)
    k2 = h*f_driven(r+0.5*k1,t+0.5*h)
    k3 = h*f_driven(r+0.5*k2,t+0.5*h)
    k4 = h*f_driven(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6

plt.plot(tpoints,theta_points,'k')
plt.xlabel(r'Time(s)', fontsize = 24)
plt.ylabel(r'$\theta(t)$', fontsize = 24)
plt.show()
plt.clf()


#8.5(b)
Omega = 10 #s^-2

#in order to find resonant Omega value:
#for W in [0.1,1,10,35,50]:
#   Omega = W
#indent the rest of the code to include execute it for all values of W to get a sense of a good value for Omega
theta_points = []
initial_theta = 0.0 #in degrees. start with the pendulum handing straight down.
initial_omega = 0.0 #m/s. start at rest
#convert theta to radians
initial_theta = initial_theta *pi/180

#vectorize the differential equation
r = [initial_theta, initial_omega]

for t in tpoints:
    theta_points.append(r[0]*180/pi) #in degrees
    k1 = h*f_driven(r,t)
    k2 = h*f_driven(r+0.5*k1,t+0.5*h)
    k3 = h*f_driven(r+0.5*k2,t+0.5*h)
    k4 = h*f_driven(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6

plt.plot(tpoints,theta_points,'k')
plt.xlabel(r'Time(s)', fontsize = 24)
plt.ylabel(r'$\theta(t)$', fontsize = 24)
plt.show()
plt.clf()