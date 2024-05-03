#Nina Martinez Diers
#Nina Martinez Diers
import numpy as np
from random import random, randrange
from math import exp, cos
import matplotlib.pyplot as plt

#use latex fonts
plt.rc('text',usetex=True)
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

#Function to calculate energy of array of spin-states with given J.
def energy(J, spin_array):
    '''#over/under neighbors
    spin1 = np.delete(spin_array,0,0)
    spin2 = np.delete(spin_array,19,0)
    diff1 = np.subtract(spin1,spin2)

    #side-by-side neighbors
    spin3 = np.delete(spin_array,0,1)
    spin4 = np.delete(spin_array,19,1)
    diff2 = np.subtract(spin3,spin4)'''
    
    #Periodic Boundary Conditions:
    #over/under neighbors
    spin1 = np.append(spin_array, np.atleast_2d(spin_array[0,:]), 0)
    spin1 = np.delete(spin1,0,0)
    diff1 = np.subtract(spin_array,spin1)

    #side-by-side neighbors
    spin2 = np.append(spin_array, np.atleast_2d(spin_array[:,0]).T, 1)
    spin2 = np.delete(spin2,0,1)
    diff2 = np.subtract(spin_array,spin2)
    for i in range(19):
        for j in range(19):
            diff1[i,j]= cos(diff1[i,j])
            diff2[i,j]= cos(diff2[i,j])

    return -J*(np.sum(diff1)+np.sum(diff2))

        
L = 20
N = 10000 #Number of monte-carlo time-steps
J = 1
kB = 1
Temp = np.arange(0.2,1.6,0.1)
sysE = [] #energy of system for diff temp vals

for T in Temp:
    #Make 20x20 lattice randomly populated with ordered distribution of spin states:
    #adapted from lattice.py in Mark Newman's Computational Physics
    spins = np.zeros((L,L), dtype=float)
    accepts = []

    #Energy
    E = energy(J, spins)
    Energy = [E]
    
    beta = 1/(kB*T)
    
    for n in range(N):
        accept = False
        E = energy(J,spins)

        #randomly pick i,j element to switch spin
        i = randrange(20)
        j = randrange(20)
        #generate new spin θi for (i,j) element in spins_test
        spin_save = spins[i,j]
        spins[i,j] = 2*np.pi*random()
        #calculate new energy
        E_test = energy(J,spins)
        dE = E_test - E

        #determine whether to accept or reject:
        if E_test <= E:
            accept = True
        else:
            if random() < exp(-beta*dE):
                accept = True
            else:
                #reject, revert spins to original and restart calcs
                spins[i,j] = spin_save

        #append energy to list
        Energy.append(energy(J,spins))


    #Average last 50 vals of Energy to find the equilibrated energy of the system
    sysE.append(np.average(Energy[-50:])/L**2)

#Plot system energy as a function of temperature  
plt.plot(Temp,sysE,'k')
plt.xlabel(r'Temperature', fontsize = 24)
plt.ylabel(r'System Energy', fontsize = 24)
plt.show()