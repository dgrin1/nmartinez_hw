#Nina Martinez Diers
#Plot Error of an integral as a function of N, identify scaling with N .....(log /log) — internal or scipy comparison
import numpy as np
from math import log10
import matplotlib.pyplot as plt
#use latex fonts
plt.rc('text',usetex=True)
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

#Using Simpson's rule to calculate an integral and evaluate the error of that integral
#integral we are going to be evaluating: sine(x)

#limits of integration:
a = 0
b = np.pi*6

#define integrand
def integrand(x):
    fx = np.sin(x)
    return fx
#3rd derivative
def integrand_3deriv(x):
    return -integrand(x)


#(h/3)[f(a)+f(b)+4sum(f(a+kh))+2sum(f(a+kh))]

#define simpson's rule
#terms are function, lower and upper limits of integration, number of bins
def simpsons(f,err,a,b,N):
    h = (b-a)/N
    epsilon = abs((h**4/180)*(err(a)-err(b))) #abs value of error of the simpson's integration
    oddsum = 0
    evensum = 0
    for k in range(1,N,2):
        oddsum += f(a+k*h)
    for k in range(2,N,2):
        evensum += f(a+k*h)
    I = (h/3)*(f(a)+f(b)+4*oddsum +2*evensum)
    return I, epsilon


#calculate error of the integral as a function of N
N = np.arange(1,1001)
error = []
for i in N:
    I,epsilon = simpsons(integrand,integrand_3deriv,a,b,i)
    error.append(epsilon) #use simpsons to calculate the error of the integral

plt.plot(N,error,'k')
plt.xlabel(r'Log10 of number of bins',fontsize=16)
plt.ylabel(r"Log10 of Simpson's rule error for sin($x$)",fontsize=16)
plt.xscale('log')
plt.yscale('log')
plt.text(N[0],error[0],'({}, {})'.format(N[0],error[0]))
plt.text(N[10],error[10],"slope = "+str(log10(error[900])-log10(error[10]) / (log10(N[900])-log10(N[10]))))
plt.show()