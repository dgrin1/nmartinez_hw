#(a) relaxation method to solve for x = 1-exp(-cx) for c=2 with an accuracy of at least 1e-6
from math import exp
from numpy import std, arange
import matplotlib.pyplot as plt
x = 1
c = 2
epsilon = 1
x_last5 = [] #create list of the last 5 calculations for x
while epsilon > 1e-7: #ensure accuracy to the 6th decimal place
    x = 1 - exp(-c*x)
    x_last5.append(x)
    if len(x_last5) > 5:
        x_last5.pop(0) #remove first element of list to keep last 5 calculations of x
        #standard deviation of last 5 values of x. loop will stop when is accuracy is less than 1e-7
        epsilon = std(x_last5)
print(x_last5, epsilon)


#(b) modify for values of 0≤c≤3 in steps of 0.1. Make a plot of x(c)
def relax_x(c):
    x = 1
    epsilon = 1
    x_last5 = [] #create list of the last 5 calculations for x
    while epsilon > 1e-7: #ensure accuracy to the 6th decimal place
        x = 1 - exp(-c*x)
        x_last5.append(x)
        if len(x_last5) > 5:
            x_last5.pop(0) #remove first element of list to keep last 5 calculations of x
            #standard deviation of last 5 values of x. loop will stop when is accuracy is less than 1e-7
            epsilon = std(x_last5)
    return x

#set fonts
plt.rc('text',usetex=True)
#plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font', **{'family':'serif','serif':['Times New Roman']})

c = arange(0,3.01,0.01) #array of c-values from 0≤c≤3 in steps of 0.1
x_c = []
for i in c:
    x_c.append(relax_x(i)) #list of x as a function of c using relaxation method

#plot x as a function of c
plt.plot(c,x_c,'k',linewidth=3)
plt.xlabel(r'$c$', fontsize=16)
plt.ylabel(r'$x(c)$', fontsize=24)
plt.show()