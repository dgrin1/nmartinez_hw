#Exercise 2.5: Quantum Potential Step

from math import sqrt
#mass of electron:
m = 9.11e-31
#energy of electron in eV
E = 10
#energy of potential step in eV
V = 9
#reduced planck's constant in eV*s
hbar = 6.582e-16

#calculate wavevectors k1 (initial and reflected) and k2 (transmitted)
k1 = sqrt(2*m*E) / hbar
k2 = sqrt(2*m*(E-V)) / hbar

#calculate the probability of transmission and reflection 
T = 4*k1*k2 / (k1+k2)**2
R = ((k1-k2) / (k1+k2))**2

print('The probability of transmission is', T)
print('The probability of reflection is', R)