#rumke-kutta 4th order
from math import sin
from numpy import arange
from pylab import plot,xlabel,ylabel,show

#acceleration due to air resistance function
#acceleration = g - bv
def f(v,t):
    g = 10 #m s^-2
    b = 1 #s^-1
    return g - b*v

a = 0.0
b = 10.0
N = 100000
h = (b-a)/N

tpoints = arange(a,b,h)
vpoints = []
v = 0.0 #start at v = 0 but accelerate up to terminal v

for t in tpoints:
    vpoints.append(v)
    k1 = h*f(v,t)
    k2 = h*f(v+0.5*k1,t+0.5*h)
    k3 = h*f(v+0.5*k2,t+0.5*h)
    k4 = h*f(v+k3,t+h)
    v += (k1+2*k2+2*k3+k4)/6

plot(tpoints,vpoints)
xlabel("t")
ylabel("v(t)")
show()

#asymptotically, a = 0
#results with g = 10 and b = 1: terminal velocity is 10m/s