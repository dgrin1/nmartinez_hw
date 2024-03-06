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
    userErr = 0.01 #change this to user input when code is working
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

#incorporate error calculation once debugged with sine
def ninaCosine(x):
    i = 0
    cosx = 0 #this will be the taylor sum with i terms added.
    while i <= 9:
        iTerm = (-1)**i * x**(2*i) / ninaFactorial(2*i) #calculate ith term (starting with i = 0)
        cosx += iTerm #add ith term to sinx
        i += 1
    return cosx

def ninaTangent(x):
    tanx = ninaSine(x) / ninaCosine(x)
    return tanx

#create list of values to test sine function
#xvals = [0,pi/2,pi,3*pi/2,2*pi]
xvals = np.linspace(0,5*2*pi,num=5000,endpoint=True) #calculate sinx over 5 periods

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