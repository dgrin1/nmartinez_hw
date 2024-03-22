#(b) evaluate H_mn for an electon in a 1D potential well with a width of Å and a=10eV
from math import sqrt
from numpy import pi, zeros, linalg, linspace, sin
import matplotlib.pyplot as plt

#arbitrary numbers to evaluate H_mn
m = 1
n = 3

#define constants
M = 9.1094e-31 #kg - mass of electron
L = 5e-10 #m
a = 10 / 6.2415e18 #converted eV to J
hbar = 1.0545e-34 #J⋅s

if m == n:
    H = (1/2)*((hbar*pi*n/L)**2/(M) + a)
    #print('m==n')
elif m%2 == 0 and n%2 == 1:
    H = -2*a*m*n/(m**2 - n**2)**2 * (2/pi)**2
    #print('m%2!=n%2')
elif m%2 == 1 and n%2 == 0:
    H = -2*a*m*n/(m**2 - n**2)**2 * (2/pi)**2
    #print('m%2!=n%2')
elif m%2 == 0 and n%2 == 0:
    H = 0
    #print('m%2==n%2')
elif m%2 == 1 and n%2 == 1:
    H = 0
    #print('m%2==n%2')
else:
    print('error')

print(H)

#(c) create 10x10 matrix for m,n 1->10 
def calcH(m,n):
    if m == n:
        H = (1/2)*((hbar*pi*n/L)**2/(M) + a)
    elif m%2 == 0 and n%2 == 1:
        H = -2*a*m*n/(m**2 - n**2)**2 * (2/pi)**2
    elif m%2 == 1 and n%2 == 0:
        H = -2*a*m*n/(m**2 - n**2)**2 * (2/pi)**2
    elif m%2 == 0 and n%2 == 0:
        H = 0
    elif m%2 == 1 and n%2 == 1:
        H = 0
    else:
        print('error')
    return H

#create 10x10 matrix of zeroes and populate with values of H_mn
matrixH = zeros([10,10])
for m in range(1,11):
    for n in range(1,11):
        matrixH[m-1,n-1] = calcH(m,n)
#the resulting matrixH is symmetric

#print first 10 energy levels of quantum well with this approximation
#calculate eigenvalues using numpy.linalg. These eigenvalues are the energies of each energy level
eigenvalues = linalg.eigvalsh(matrixH)
#print first 10 energy levels of quantum well with this approximation
for i in range(len(eigenvalues)):
    print('Energy of n =', i+1, 'is', eigenvalues[i]*6.2415e18, 'eV')


#(d) create 100x100 matrix for m,n 1->100 and calculate first 10 eigenvalues 
#create 10x10 matrix of zeroes and populate with values of H_mn
matrixH = zeros([100,100])
for m in range(1,101):
    for n in range(1,101):
        matrixH[m-1,n-1] = calcH(m,n)
#the resulting matrixH is symmetric

#print first 10 energy levels of quantum well with this approximation
#calculate eigenvalues using numpy.linalg. These eigenvalues are the energies of each energy level
eigenvalues = linalg.eigvalsh(matrixH)
#print first 10 energy levels of quantum well with this approximation
for i in range(10):
    print('Energy of n =', i+1, 'is', eigenvalues[i]*6.2415e18, 'eV')

#(e) calculate wavefunction for n = 1,2,3
#calculate eigenvectors of matrixH. We use these to calculate the wavefunction for each n
evals, evectors = linalg.eigh(matrixH)

#make graph of of psi(x) with each n
#latex fonts
plt.rc('text',usetex=True)
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

#set xvalues
xvals = linspace(0,L,100)

#calculate psi_n(x) and plot
for n in range(3):
    psi_n = []
    for x in xvals:
        s = 0
        for i in range(100): #because matrixH is 100x100 matrix, eigenvectors are 100 elements long
            s += evectors[n,i] * sin(pi*x*(i+1)/L) #nth row in evectors is the nth eigenvector.
        psi_n.append(s * sqrt(L/2)) #multiply sum by sqrt(L/2) to normalize the wavevector
    plt.plot(xvals, psi_n, label = r'n ='+str(n+1), linewidth=3)
plt.xlabel(r'$x$', fontsize = 16)
plt.ylabel(r'$\psi^n(x)$', fontsize = 24)
plt.legend()
plt.show()