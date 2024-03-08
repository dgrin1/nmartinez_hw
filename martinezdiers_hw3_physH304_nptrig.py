#compare sine function to numpy sine function
import numpy as np
import matplotlib.pyplot as plt
from math import pi

#sk user for desired fractional error for sine calculation
inputErr = float(input("Desired fractional error:"))

#return sine of given x
def ninaFactorial(n):
    factorial = 1
    while n > 0:
        factorial *= n
        n -= 1
    if n == 0:
        return factorial

def ninaSine(x):
    i = 0
    sinx = 0 #this will be the taylor sum with i terms added.
    userErr = inputErr
    fracErr = userErr + 1 #ensures while loop will start

    while fracErr > userErr: #if fractional error is greater than user error, continue with taylor series expansion
        if i == 0: #this sets initial fractional error to zero after the loop has been initiated instead of arbitrary value
            fracErr = 0
        iTerm = (-1)**i * x**(2*i + 1) / ninaFactorial(2*i + 1) #calculate ith term (starting with i = 0)
        sinx += iTerm #add ith term to sinx
        #check if sinx == 0 using algorithm from Newman on pg129
        epsilon = 1e-12
        if abs(sinx) < epsilon: #if sinx is sufficiently close to zero, break out of while loop (such as at multiples of pi)
            break
        else:
            fracErr = abs(iTerm/sinx)
        i += 1
    return sinx, fracErr

#create list of values to test sine function
xvals = np.linspace(0,5*2*pi,num=5000,endpoint=True) #calculate sinx over 5 periods

ninaSineWave = [] #list of my sine function outputs to plot
numpySineWave = [] #list of numpy sine function outputs to plot
sineDifference = [] #calculate the difference between the two functions
for i in xvals:
    a, b = ninaSine(i) #calculate sine with fractional error of xvals using sine function; 
    ninaSineWave.append(a)
    numpySineWave.append(np.sin(i))
    sineDifference.append(a-np.sin(i))

#set latex font
plt.rc('text',usetex=True)
plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})

#make plot of cosine wave
fig, ((ninawave), (numpywave), (diff)) = plt.subplots(nrows = 3, ncols = 1)
ninawave.plot(xvals, ninaSineWave, 'k', lw=3, label=r'my sine') 
numpywave.plot(xvals,numpySineWave, 'k', lw=3, label=r'numpy sine')
diff.plot(xvals,sineDifference, 'k', lw=3, label=r'difference between sine algorithms')

#labels in latex
plt.xlabel(r'$x$', fontsize=24)
ninawave.set_ylabel(r'my sine', fontsize=24)
numpywave.set_ylabel(r'numpy sine', fontsize=24)
diff.set_ylabel(r'difference', fontsize=24)

#show plot
plt.show()