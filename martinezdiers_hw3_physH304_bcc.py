#Nina Martinez Diers
#Madelung constant for a body-centered cubic lattice
from math import sqrt
import matplotlib.pyplot as plt

r = 1/2 #radius of each atom. This gives us same atom size as for simple cubic lattice
a = (4*r)/sqrt(3) #length of unit cell side. Because body diagonal = diameter of cl ion + diameter cs ion
x = range(11)
Mlist = []

def MadelungBCC(N): # 2N is number of lattice sites per side of crystal
    M = 0
    for i in range(-N,N+1):
        for j in range(-N,N+1):
            for k in range(-N,N+1):
                mCenter = -1 / sqrt((a*i + a/2)**2 + (a*j + a/2)**2 + (a*k + a/2)**2)
                if i == 0 and j == 0 and k == 0: #center atom is reference point
                    mCorner = 0
                else:
                    mCorner = 1 / (a*sqrt(i**2 + j**2 + k**2)) 
                M += mCenter + mCorner
    return M

xlist = []
for i in x:
    xlist.append(2*i)
    Mlist.append(MadelungBCC(i))

#use latex fonts
plt.rc('text',usetex=True)
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

plt.plot(xlist,Mlist,'-k')
plt.xlabel(r'Number of Edge Lattice sites')
plt.ylabel(r'Madelung Constant')
plt.show()