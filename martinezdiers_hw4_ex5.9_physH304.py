#Nina Martinez Diers
#Calculate Heat capacity as a function of Temperature.
from numpy import ones,copy,cos,tan,pi,linspace
from math import exp
import matplotlib.pyplot as plt

#Functions to calculate points and weights of gaussian quadrature. Will be called to calculate uncertainty in the wave function
#gaussxw and gaussxwab functions written by Mark Newman, June 4, 2011.
def gaussxw(N):

    # Initial approximation to roots of the Legendre polynomial
    a = linspace(3,4*N-1,N)/(4*N+2)
    x = cos(pi*a+1/(8*N*N*tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = ones(N,float)
        p1 = copy(x)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x,w

def gaussxwab(N,a,b):
    x,w = gaussxw(N)
    return 0.5*(b-a)*x+0.5*(b+a),0.5*(b-a)*w

#calculate heat capacity of 1000 cm^3 of solid aluminum
def cv(T):
    V = 1000 / (100)**(-3) #in meters^3
    numberDensity = 6.022e28 #in meters^-3
    kB = 1.3806503e-23 #boltzmann constant in m^2*kg*s^-2*K^-1
    DeByeTemperature = 428 #in Kelvin
    N = 50 #sample points for gaussian quadrature
    x,w = gaussxwab(N,0,DeByeTemperature/T)
    s = 0.0
    for k in range(N):
        s+= w[k] * x[k]**4 * exp(x[k]) / (exp(x[k])-1)**2
    return (9*V*numberDensity*DeByeTemperature*kB*s)


#generate data for graph of heat capacity as a function of temperature
T = range(5,501)
Cv = []
for t in T:
    Cv.append(cv(t))

#Make plot
#latex fonts
plt.rc('text',usetex=True)
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

#plot data
plt.plot(T,Cv,'k')

#define plot specifics
plt.xlabel(r'Temperature (K)', fontsize = 16)
plt.ylabel(r'Heat Capacity ($m^2*kg*s^{-2}*K^{-1}$)', fontsize = 24)
plt.show()