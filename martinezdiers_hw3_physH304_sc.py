#Nina Martinez Diers
#Madelung constant problem (NaCl)
from math import sqrt, cosh
from numpy import pi
import matplotlib.pyplot as plt

x = range(11)
Mlist = []

#method 1 to calculate madelung constant with sum
def Odd(x):
    if x % 2 ==0: #if x is divisible by 2
        return False #x is not even
    else:
        return True #x is odd

def Madelung(N):
    M = 0
    for i in range(-N,N+1):
        for j in range(-N,N+1):
            for k in range(-N,N+1):
                if i == 0 and j == 0 and k == 0:
                    m = 0
                else:
                    m = 1 / sqrt(i**2 + j**2 + k**2)
                if Odd(i+j+k) == True:
                    m = -m
                M += m
    return M

xlist = []
for i in x:
    xlist.append(2*i)
    Mlist.append(Madelung(i))

#latex fonts


plt.plot(xlist,Mlist,'-k')
plt.xlabel(r'Atoms on Crystal Edge')
plt.ylabel(r'Madelung Constant for Crystal')
plt.show()

#alternative method for calculating madelung constant with hyperbolic trig
Mlist2 = []
def Madelung2(N):
    M = 0
    for m in range(1,N+1,2):
        for n in range(1,N+1,2):
            M += 1 / (cosh(pi * sqrt(m**2 + n**2) / 2))**2
    return 12*pi*M

for i in x:
    Mlist2.append(Madelung2(i))

plt.plot(xlist,Mlist2,'-b')
plt.xlabel(r'Atoms on Crystal Edge')
plt.ylabel(r'Madelung Constant for Crystal')
plt.show()