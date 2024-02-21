import numpy as np
import matplotlib.pyplot as plt
from math import pi

#return sine of given x
def ninaFactorial(n):
    factorial = 1
    while n > 0:
        factorial *= n
        n -= 1
    if n == 0:
        return factorial
#still need to debug

def ninaSine(x):
    i = 0
    sinx = 0 #this will be the taylor sum with i terms added.
    userErr = 1 #change this to user input when code is working
    fracErr = userErr + 1 #ensures while loop will start
    #sinx is not being expanded past i = 1
def ninaSine(x):
    i = 0
    sinx = 0 #this will be the taylor sum with i terms added.
    userErr = 1 #change this to user input when code is working
    fracErr = userErr + 1 #ensures while loop will start
    
    while fracErr > userErr:
        iTerm = (-1)**i * x**(2*i + 1) / ninaFactorial(2*i + 1) #calculate ith term (starting with i = 0)
        sinx += iTerm #add ith term to sinx
        if sinx <= 1e-7 and sinx >= -1e-7: #set high fractional error to break loop. This allows function not to break at multiples of 2pi
            fracErr = 100
        else:
            fracErr = abs(iTerm/sinx)
        if fracErr >= 100: #set high fractional error to break loop. This allows function not to break at multiples of 2pi
            break
        i += 1
    return sinx, fracErr


#create list of values to test sine function
#xvals = [0,pi/2,pi,3*pi/2,2*pi]
xvals = np.linspace(0,2*pi,num=100,endpoint=False) #True)

sineWave = [] #list of sine function outputs to plot
error = []
for i in xvals:
    a, b = ninaSine(i) #calculate sine with fractional error of xvals using sine function; 
    sineWave.append(a)
    error.append(b)

print(sineWave)
wave = plt.plot(xvals,sineWave,"k")
err = plt.plot(xvals,error, "r")

#plot title
plt.title("Sine Wave")
#axis labels
plt.xlabel("x")
plt.ylabel("sin(x)")
#Make plot
plt.show()