#Exercise 8.2: The Lotka-Volterra Equations (Predator-Prey Interactions)
#Nina Martinez Diers
#adapt rk4.py to lotka-volterra equations
from math import sin, cos
from numpy import arange, array, pi
import matplotlib.pyplot as plt

#lotka-volterra constants
alpha = 1
beta = 0.5
gamma = beta
delta = 2

#initial conditions
x = 2 #starting population of prey
y = 2 #starting population of predator
#vectorize the differential equation
r = [x,y]

def f(r): #no explicit time-dependence
    x = r[0]
    y = r[1]
    fx = alpha*x - beta*x*y
    fy = gamma*x*y-delta*y
    return array([fx, fy], float)

a = 0.0
b = 30.0
N = 100000
h = (b-a)/N

tpoints = arange(a,b,h)
xpoints = []
ypoints = []

for t in tpoints:
    xpoints.append(r[0])
    ypoints.append(r[1])
    k1 = h*f(r)
    k2 = h*f(r+0.5*k1)
    k3 = h*f(r+0.5*k2)
    k4 = h*f(r+k3)
    r += (k1+2*k2+2*k3+k4)/6

#Make plot
#latex fonts
plt.rc('text',usetex=True)
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

plt.plot(tpoints,xpoints,'k', label=r'prey population')
plt.plot(tpoints,ypoints,'r',label=r'predator population')
plt.xlabel(r"Time", fontsize = 24)
plt.ylabel(r"Population", fontsize  = 24)
plt.legend()
plt.show()