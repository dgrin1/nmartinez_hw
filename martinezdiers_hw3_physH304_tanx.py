import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, pi

#ask user for desired fractional error:
inputErr = float(input("What is the desired fractional error?: ")) 

def ninaFactorial(n):
    factorial = 1
    while n > 0:
        factorial *= n
        n -= 1
    if n == 0:
        return factorial

#set error margin like Newman on pg129
epsilon = 1e-12

def ninaSine(x):
    i = 0
    sinx = 0 #this will be the taylor sum with i terms added.
    userErr = inputErr #change this to user input when code is working
    fracErr = userErr + 1 #ensures while loop will start

    while fracErr > userErr: #if fractional error is greater than user error, continue with taylor series expansion
        if i == 0: #this sets initial fractional error to zero after the loop has been initiated instead of arbitrary value
            fracErr = 0
        iTerm = (-1)**i * x**(2*i + 1) / ninaFactorial(2*i + 1) #calculate ith term (starting with i = 0)
        sinx += iTerm #add ith term to sinx
        #check if sinx == 0 using algorithm from Newman on pg129
        if abs(sinx) < epsilon: #if sinx is sufficiently close to zero, break out of while loop (such as at multiples of pi)
            break
        else:
            fracErr = abs(iTerm/sinx)
        i += 1
    return sinx, fracErr

def ninaCosine(x):
    i = 0
    cosx = 0 #this will be the taylor sum with i terms added.
    userErr = inputErr #change this to user input when code is working
    fracErr = userErr + 1 #ensures while loop will start
    
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

def ninaTangent(x):
    sinx,sinErr = ninaSine(x)
    cosx,cosErr = ninaCosine(x)
    tanx = sinx/cosx
    #Using Eq. 4.3 from Newman on pg130 for the variance when dividing variables:
    #sigma^2 = (tanx)^2 * [(sinErr)^2/(sinx)^2 + (cosErr)^2/(cosx)^2]

    if abs(sinx) > epsilon and abs(cosx) > epsilon:
        fracErr = sqrt(tanx**2 * (sinErr**2/sinx**2 + cosErr**2/cosx**2))
    elif abs(sinx) < epsilon and abs(cosx) > epsilon:
        fracErr = sqrt(tanx**2 * cosErr**2/cosx**2)
    elif abs(sinx) > epsilon and abs(cosx) < epsilon:
        fracErr = sqrt(tanx**2 * sinErr**2/sinx**2)
    else:
        fracErr = 0
    return tanx, fracErr

#create list of values to test sine function
#xvals = [0,pi/2,pi,3*pi/2,2*pi] #test points
xvals = np.linspace(0,5*pi,num=1000,endpoint=True)

tangentWave = [] #list of tangent function outputs to plot
error = []
for i in xvals:
    tan, tanErr = ninaTangent(i) #calculate tangent with fractional error of xvals using tangent function; 
    tangentWave.append(tan)
    error.append(tanErr)

#set latex font
plt.rc('text',usetex=True)
plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})

#making figure with two y-axes: https://pythonguides.com/matplotlib-two-y-axes/
#make plot of tangent wave
fig, ax1 = plt.subplots()
wave = ax1.plot(xvals, tangentWave, 'k', lw=3) 

#labels in latex
ax1.set_xlabel(r'$x$',fontsize=16)
ax1.set_ylabel(r'$tan(x)$',fontsize=24)

#plot fractional error with 2nd y axis
ax2 = ax1.twinx()
err = ax2.plot(xvals,error, 'r')

#labels in latex
ax2.set_ylabel(r'Fractional Error', fontsize=24, color = 'r')

#show plot
plt.show()
#plt.savefig("tanxFig.pdf")