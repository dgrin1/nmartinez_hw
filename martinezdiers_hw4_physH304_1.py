#Simpson's rule
from math import exp
import numpy as np
import matplotlib.pyplot as plt

#define error function integrand
def erf_integrand(t):
    fx = exp(-t**2)
    return fx
#3rd derivative
def erf_integrand_3deriv(t):
    return 4*t*erf_integrand(t)*(3-2*t**2)

#determine number of bins

#(h/3)[f(a)+f(b)+4sum(f(a+kh))+2sum(f(a+kh))]

#first step: for 
#terms are function, lower and upper limits of integration, number of bins
def simpsons(f,a,b,N):
    h = (b-a)/N
    epsilon = (h**4/180)*(erf_integrand_3deriv(a)-erf_integrand_3deriv(b))
    oddsum = 0
    evensum = 0
    for k in range(1,N,2):
        oddsum += f(a+k*h)
    for k in range(2,N,2):
        evensum += f(a+k*h)
    I = (h/3)*(f(a)+f(b)+4*oddsum +2*evensum)
    return I, epsilon

#for erf, a = x0 = 0 and b = xf = x
#calculate erf as a function of x
x0 = 0
xf = 3
N = 30 #this will make h=0.1
x = np.linspace(x0,xf,N,endpoint = True)
erf = []
err = []
for i in x:
    I,epsilon = simpsons(erf_integrand,x0,i,N)
    erf.append(I) #use simpsons to calculate the error function integral
    err.append(epsilon)

print(err)
#use latex fonts
plt.rc('text',usetex=True)
plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})

#make plot
plt.plot(x,erf,'k', label=r'Error Function')
plt.plot(x,err,'r', label=r'Approximation Error')
plt.xlabel(r'$x$', fontsize = 16)
plt.ylabel(r'$E(x)$', fontsize = 24)
plt.legend()
plt.show()