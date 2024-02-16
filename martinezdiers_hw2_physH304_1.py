#Exercise 3.1(a)
import numpy as np
from math import pi
import matplotlib.pyplot as plt

'''pseudocode:
load data: data = np.loadtxt(sunspots file, float)
select columns for variables: time = data[:,0] and sunspots =  data[:,1]
make plot: (add title, labels)
plt.plot(time,sunspots)
plt.title()
plt.xlabel() and plt.ylabel()
plt.show()
'''
data = np.loadtxt("sunspots.txt", float) #read data file sunspots
time = data[:,0] #number of months since Jan 1749
sunspots = data[:,1] #number of sunspots observed on the Sun during that month

plt.plot(time,sunspots,"k")
#plot title
plt.title("Number of Sunspots Per Month Since January 1749")
#axis labels
plt.xlabel("Time") #find how to label each month?
plt.ylabel("Number of Sunspots")
#Make plot
plt.savefig("sunspotsFig1.pdf")
plt.clf() #citation: https://www.activestate.com/resources/quick-reads/how-to-clear-a-plot-in-python/

#Exercise 3.1(b)
#Objective: modify the program to only include the first 1000 data points on the graph
#pseudocode:
#when create time and sunspots data sets, add only the first 1000 rows instead of all rows
#dataCol = data[0:1000, column]

import numpy as np
from math import pi
import matplotlib.pyplot as plt

data = np.loadtxt("sunspots.txt", float) #read data file sunspots
time = data[0:1000,0] #number of months since Jan 1749
sunspots = data[0:1000,1] #number of sunspots observed on the Sun during that month

plt.plot(time,sunspots, "k")
#plot title
plt.title("Number of Sunspots Per Month from January 1749 to April 1832")
#axis labels
plt.xlabel("Time") #find how to label each month?
plt.ylabel("Number of Sunspots")
#Make plot
plt.savefig("sunspotsFig2.pdf")
plt.clf()

#Exercise 3.1(c)
'''pseudocode:

def runAvg(list): #calculate running average of all elements in a list
    sum  = 0
    for i in list:
        sum += list[i]
    avg = sum / (2*len(list)+1)
    return avg


Yk = [0,0,0,0] #make list of initial Yk values. Because the running average uses previous 5 values, the running average for the first 4 months is not applicable so we just set the starting value to zero.
for k in range(1,1001):
    append list r with sunspots[k]
    if len(r) > 5:
        r = r[1:5]
    #calculate running average of sunspots
    if len(r) == 5:
        call running average function
        append to list Yk
'''
r = 5 #number of months 
yk = [] #list of sunspots to calculate average
Yk = [] #list of running average
def Avg(list): #calculate the average of all elements in a list
    sum  = 0
    for i in range(0,r):
        sum += list[i] #calculate sum of all sunspots over 5 month period
    avg = sum / (2*len(list)+1) #divide sum by normalization constant to calculate average
    return avg

#Calculate Running Average of sunspots
for k in range(0,1000): #only perform running average for the first 1000 months
    yk.append(sunspots[k]) #append list k with sunspots[n]
    if len(yk) < r: #Because the running average uses previous 5 values, the running average for the first 4 months is not applicable so we just set the starting value to zero.
        Yk.append(0)
    if len(yk) > r:
        yk = yk[1:r+1] #if k contains 5 elements, keep only the last 5 (remove the first element)
    if len(yk) == r: #when updated list is the correct length:
        Yk.append(Avg(yk)) #Add element to the list with the new value for the running average

#Make plot with sunspots per month (for first 1000 months) and running average
plt.plot(time,sunspots,"k", label = "Number of Sunspots Per Month")
plt.plot(time,Yk,"b", label = "5-Month Running Average")
#plot title
plt.title("Sunspots from January 1749 to April 1832")
#axis labels
plt.xlabel("Time") #find how to label each month?
plt.ylabel("Number of Sunspots")
#plot legend
plt.legend(loc = 1)
#Make plot
plt.savefig("sunspotsFig3.pdf")
plt.clf()