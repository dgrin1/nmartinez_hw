#Hermite polynomials
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, exp, factorial

#Functions to calculate points and weights of gaussian quadrature. Will be called to calculate uncertainty in the wave function
#gaussxw and gaussxwab functions written by Mark Newman, June 4, 2011.
def gaussxw(N):

    # Initial approximation to roots of the Legendre polynomial
    a = np.linspace(3,4*N-1,N)/(4*N+2)
    x = np.cos(np.pi*a+1/(8*N*N*np.tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = np.ones(N,float)
        p1 = np.copy(x)
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

#Define function for hermite polynomial
#H(n+1,x) = 2xH(n,x) - 2nH(n-1,x) is equivalent to H(n,x) = 2xH(n-1,x) - 2(n-1)H(n-2,x)

#First try at writing Hermite polynomial function. It works, but is inefficient and can be very slow.
#def H(n,x):
#    if n == 0:
#        Hermite = 1
#    elif n ==1:
#        Hermite = 2*x
#    else:
#        Hermite = 2*x*H(n-1,x) - 2*(n-1)*H(n-2,x)
#    return Hermite

#Better Hermite polynomial Function:
def H(n,x):
    #Make list of Hermite polynomials with values for n = 0 and n = 1
    Hermite = [1,2*x]
    #Expand list calculating all values from n = 0 through desired n
    while len(Hermite) <= n:
        index = len(Hermite)
        Hermite.append(2*x*Hermite[index-1] - 2*(index-1)*Hermite[index-2])
    #return last element in the list, which is nth Hermite polynomial for given value of x
    return Hermite[n]

#use latex fonts
plt.rc('text',usetex=True)
plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})

x = np.linspace(-4,4,100,endpoint = True)

#check that Hermite polynomial algorithm works:
#for n in range(0,4): 
#    Hermite = []
#    for i in x:
#        Hermite.append(H(n,i))
#    plt.plot(x,Hermite,label = n)
#plt.legend()
#plt.show()

#Calculate HO wavefuction for 0≤n≤3 from -4≤x≤4
for n in range(0,4):
    psi = []
    for i in x:
        psi.append((2**n * factorial(n) * np.pi**(1/2))**(-1/2)*exp(-i**2 / 2)*H(n,i))
    plt.plot(x,psi,label = r'n = '+str(n))

#make plot
plt.xlabel(r'$x$', fontsize = 16)
plt.ylabel(r'$\psi_n$', fontsize = 24)
plt.legend()
plt.show()
plt.clf()

#Calculate wavefunction for n = 30 from -10≤x≤10
n = 30
psi = []
x = np.linspace(-10,10,500, endpoint = True)
for i in x:
    psi.append((2**n * factorial(n) * np.pi**(1/2))**(-1/2)*exp(-i**2 / 2)*H(n,i))
plt.plot(x,psi,'k')
#make plot
plt.xlabel(r'$x$', fontsize = 16)
plt.ylabel(r'$\psi_{30}$', fontsize = 24)
plt.show()

#Calculate quantum uncertainty using gaussian quadrature
#Calculate the uncertainty in psi_5(x)
n = 5
#perform change of variables for psi(x)--> psi(tan(z))
#Use gaussian quadrature to calculate <x^2> of the wavefunction psi over 100 points from -infinity to infinity
#This is the same as calculating <tan^2(z)> over 100 points from -pi/2 to pi/2
N = 100 
a = -np.pi/2
b = np.pi/2
z,w = gaussxwab(N,a,b)
s = 0.0 
for k in range(N):
    # s+=w[k] * f(z[k])
    #f(z) = tan(z)**2 * (psi(tan(z)))**2 / cos(z)**2
    s+=w[k] * np.tan(z[k])**2 * ((2**n * factorial(n) * np.pi**(1/2))**(-1/2)*exp(-np.tan(z[k])**2 / 2)*H(n,np.tan(z[k])))**2 / (np.cos(z[k])**2)


rms_uncertainty = s**(1/2)
print(rms_uncertainty)