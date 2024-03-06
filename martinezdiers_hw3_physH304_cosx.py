import numpy as np
import matplotlib.pyplot as plt
from math import pi

#return cosine of given x
inputErr = float(input("What is the desired fractional error?: "))

#calculate the factorial
def ninaFactorial(n):
    factorial = 1
    while n > 0:
        factorial *= n
        n -= 1
    if n == 0:
        return factorial

def ninaCosine(x):
    i = 0
    cosx = 0 #this will be the taylor sum with i terms added.
    userErr = inputErr #change this to user input when code is working
    fracErr = userErr + 1 #ensures while loop will start
    epsilon = 1e-12
    
    while fracErr > userErr: #if fractional error is greater than user error, continue with taylor series expansion
        if i == 0: #this sets initial fractional error to zero after the loop has been initiated instead of arbitrary value
            fracErr = 0
        iTerm = (-1)**i * x**(2*i) / ninaFactorial(2*i) #calculate ith term (starting with i = 0)
        cosx += iTerm #add ith term to cosx
        #check if cosx == 0 using algorithm from Newman on pg129
        if abs(cosx) < epsilon: #if cosx is sufficiently close to zero, break out of while loop (such as at multiples of pi)
            break
        else:
            fracErr = abs(iTerm/cosx)
        i += 1
    return cosx, fracErr
#create list of values to test sine function
xvals = np.linspace(0,5*2*pi,num=5000,endpoint=True) #calculate cosx over 5 periods

cosineWave = [] #list of cosine function outputs to plot
error = [] #list of the fractional error between i and i-1 terms
for i in xvals:
    a, b = ninaCosine(i) #calculate cosine with fractional error of xvals using cosine function; 
    cosineWave.append(a) #add these values to the lists
    error.append(b)


#set latex font
plt.rc('text',usetex=True)
plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})

#make plot of cosine wave
wave = plt.plot(xvals, cosineWave, 'k', lw=3, label=r'$cos(x)$') 
err = plt.plot(xvals,error, 'r', label=r'Fractional Error')

#labels in latex
plt.xlabel(r'$x$', fontsize=16)

#labels in latex
plt.legend(loc=1)

#show plot
plt.show()