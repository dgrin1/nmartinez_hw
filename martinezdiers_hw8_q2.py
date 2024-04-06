#Nina Martinez Diers
#Calculating the energy levels and wavefunctions for quadratic and quartic potential wells
from math import sqrt
from numpy import array,arange,linspace
import matplotlib.pyplot as plt
#adapted from Mark Newman's squarewell.py in Computational Physics

#use latex fonts
plt.rc('text',usetex=True)
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

# Constants
m = 9.1094e-31     # Mass of electron
hbar = 1.0546e-34  # Reduced Planck's
e = 1.6022e-19     # Electron charge

V0 = 50*1.602177e-19 #50eV converted to Joules
a = 1e-11 #in meters

L = 20*a #width of the quadratic well for x = -10a -> 10a
N = 1000
h = L/N

# Quadratic potential function
def V(x):
    return V0*x**2/a**2

# Wavefunction function
def f(r,x,E):
    psi = r[0]
    phi = r[1]
    fpsi = phi
    fphi = (2*m/hbar**2)*(V(x)-E)*psi
    return array([fpsi,fphi],float)

# Calculate the wavefunction for a particular energy
def solve(E):
    psi = 0.0 #to satisfy boundary condition psi(-10a)=0
    phi = 1.0
    r = array([psi,phi],float)

    for x in arange(-10*a,10*a,h):
        k1 = h*f(r,x,E)
        k2 = h*f(r+0.5*k1,x+0.5*h,E)
        k3 = h*f(r+0.5*k2,x+0.5*h,E)
        k4 = h*f(r+k3,x+h,E)
        r += (k1+2*k2+2*k3+k4)/6

    return r[0]

# Main program to find the energy using the secant method

print("Quantum energy levels for the quadratic well")
#Ground state
E1 = 0
E2 = e
psi2 = solve(E1)

target = e/1000
while abs(E1-E2)>target:
    psi1,psi2 = psi2,solve(E2)
    E1,E2 = E2,E2-psi2*(E2-E1)/(psi2-psi1)

print("E =",E2/e,"eV")
E0=E2

#1st excited state
E1=2*E0
E2=E1+e
psi2 = solve(E1)

while abs(E1-E2)>target:
    psi1,psi2 = psi2,solve(E2)
    E1,E2 = E2,E2-psi2*(E2-E1)/(psi2-psi1)
print("E_1 =",E2/e,"eV")

#2nd excited state
E1=4*E0
E2=E1+e
psi2 = solve(E1)

while abs(E1-E2)>target:
    psi1,psi2 = psi2,solve(E2)
    E1,E2 = E2,E2-psi2*(E2-E1)/(psi2-psi1)
print("E_2 =",E2/e,"eV")

#Now change the potential function:
# Quartic Potential function
def V(x):
    return V0*x**4/a**4

# Main program to find the energy using the secant method
#List of energies in Joules for first 3 energy levels
E_levels = []

print("Quantum energy levels for the quartic well")
#Ground state
E1 = 0
E2 = e
psi2 = solve(E1)

target = e/1000
while abs(E1-E2)>target:
    psi1,psi2 = psi2,solve(E2)
    E1,E2 = E2,E2-psi2*(E2-E1)/(psi2-psi1)

print("E =",E2/e,"eV")
E_levels.append(E2)
E0=E2

#1st excited state
E1=2*E0
E2=E1+e
psi2 = solve(E1)

target = e/1000
while abs(E1-E2)>target:
    psi1,psi2 = psi2,solve(E2)
    E1,E2 = E2,E2-psi2*(E2-E1)/(psi2-psi1)
print("E_1 =",E2/e,"eV")
E_levels.append(E2)
    
#2nd excited state
E1=1.5*E2
E2=E1+e
psi2 = solve(E1)

target = e/1000
while abs(E1-E2)>target:
    psi1,psi2 = psi2,solve(E2)
    E1,E2 = E2,E2-psi2*(E2-E1)/(psi2-psi1)

print("E_2 =",E2/e,"eV")
E_levels.append(E2)

#Change size of well and plot NORMALIZED wavefunction for first 3 energy levels
labels = ['Ground state', 'First excited state', 'Second excited state']
x_vals = arange(-5*a,5*a,h)

# Calculate the wavefunction for each energy state
for E in E_levels:
    psi = 0.0 #to satisfy boundary condition psi(-10a)=0
    phi = 1.0
    r = array([psi,phi],float)
    
    psi_x = [] #normalized psi(x) for box
    for x in x_vals:
        psi_x.append(r[0])#/s
        k1 = h*f(r,x,E)
        k2 = h*f(r+0.5*k1,x+0.5*h,E)
        k3 = h*f(r+0.5*k2,x+0.5*h,E)
        k4 = h*f(r+k3,x+h,E)
        r += (k1+2*k2+2*k3+k4)/6

    #Riemann integration to find the square of the integration constant
    s=0.0
    for i in psi_x:
        s+= i**2*h
    print(s)

    #normalized psi
    psi_normal = []
    for i in range(len(psi_x)):
        psi_normal.append(psi_x[i]/sqrt(s))

    plt.plot(x_vals,psi_normal, label = labels[E_levels.index(E)])

plt.xlabel(r'$x$', fontsize = 16)
plt.ylabel(r'$\psi(x)$', fontsize = 24)
plt.legend()
plt.show()