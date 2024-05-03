#Nina Martinez Diers
import numpy as np
from random import random, randrange
from math import exp
import matplotlib.pyplot as plt
from vpython import * #sphere, rate, canvas, display

#Function to calculate energy of array of spin-states with given J.
def energy(J, spin_array):
    #Fixed boundaries
    '''#over/under neighbors
    spin1 = np.delete(spin_array,0,0)
    spin2 = np.delete(spin_array,19,0)
    prod1 = spin1*spin2

    #side-by-side neighbors
    spin3 = np.delete(spin_array,0,1)
    spin4 = np.delete(spin_array,19,1)
    prod2 = spin3*spin4'''

    #periodic boundary conditions:
    #over/under neighbors
    spin1 = np.append(spin_array, np.atleast_2d(spin_array[0,:]), 0)
    spin1 = np.delete(spin1,0,0)
    prod1 = spin_array*spin1
        
    #side-by-side neighbors
    spin2 = np.append(spin_array, np.atleast_2d(spin_array[:,0]).T, 1)
    spin2 = np.delete(spin2,0,1)
    prod2 = spin_array*spin2

    return -J*(np.sum(prod1)+np.sum(prod2))

#Show an animation
#scene = canvas()
#make 2D lattice of spheres. adapted from lattice.py in Mark Newman's Computational Physics
L = 20
v = 10
R = 0.3
s = np.empty((L,L),sphere)
#Make 20x20 lattice randomly populated with equal distribution of spin states +1 or -1:
#Make visual lattice for animation with spin states blue or cyan
spins = np.ones((L,L), dtype=int)
for i in range(L):
    for j in range(L):
        s[i,j] = sphere(pos=vec(i-v,j-v,0),radius=R)
        if random() < 0.5:
            spins[i,j] = -1
            s[i,j].color = color.cyan
        else:
            s[i,j].color = color.blue

N = 1000000 #Number of monte-carlo time-steps
J = 1
T = 1
kB = 1
beta = 1/(kB*T)
E = energy(J, spins)

#Magnetization:
M = [np.sum(spins)/L**2]

for n in range(N):
    accept = False #initialize acceptance condition
    #rate(1000)
    spins_test = spins
    #randomly pick i,j element to switch spin
    i = randrange(20)
    j = randrange(20)
    #flip the sign of (i,j) element in spins_test
    spins_test[i,j] = -1*spins_test[i,j]
    #calculate new energy
    E_test = energy(J,spins_test)
    
    #determine whether to accept or reject:
    if E_test <= E:
        accept = True
    else:
        if random() < exp(-beta*(E_test - E)):
            accept = True
        else:
            #reject, accept remains false and restart calcs
            accept = False
    #if accept, update spins and change sphere colors
    if accept ==True:
        if spins_test[i,j] == 1:
            s[i,j].color = color.blue
        else:
            s[i,j].color = color.cyan
        #s[i,j].color = color.yellow #used to test spin changes better
        E = E_test
        spins = spins_test
    #calculate magnetization, append to list
    M.append(np.sum(spins)/L**2)

#Plot Magentization as a function of time  
#use latex fonts
plt.rc('text',usetex=True)
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

plt.plot(range(N+1),M,'-k')
plt.xlabel(r'Time', fontsize = 24)
plt.ylabel(r'Magnetization', fontsize = 24)
plt.show()