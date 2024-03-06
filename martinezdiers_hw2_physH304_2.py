#Exercise 3.7 The Mandelbrot Set
#generate a grid of numbers, check whether they are in the mandelbrot set, and create a density plot showing the mandelbrot set pattern
''' pseudocode:
-Use linspace to generate 100 evenly spaced points between -2 and 2 for x and y
-create 100x100 matrix of zeroes to populate with mandelbrot set numbers
-Use for loops to iterate over all (x,y) points. For each point, calculate c and check whether c is in the mandelbrot set => if z > 2, not mandelbrot; else mandelbrot
-If c is in the mandelbrot set, change (x,y) index in matrix from zero to one. 
    This will indicate that the value of c corresponding to that location in the matrix is in the mandelbrot set.
    -If c is not in the mandelbrot set, keep the zero in the matrix
    -If c is in the mandelbrot set, change the zero to a 1 in the matrix.
-This matrix shows for each xy point if the resulting c is in the mandelbrot set or not
    It is a grid of ones that form the shape of the mandelbrot pattern
-make a density plot of the matrix data using plt.imshow(matrix)
-when that works, increase  linespace to 1000 points to increase resolution'''

import numpy as np
import cmath

#Create 2000 evenly spaced x and y coordinates for the range [-2,2]
x = np.linspace(-2,2,num=1000,endpoint=True)
y = np.linspace(-2,2,num=1000,endpoint=True)

#create matrix of zeroes to initialize the mandelbrot set grid.
mdbSet = np.zeros([1000,1000])

#define function for checking if a number is in the mandelbrot set:
def mdbCheck(c):
    i = 0 #initialize counter
    z = 0 #start the sum from z=0
    while i < 100: #repeat sum 100 times to make sure |z| is always less than 2
        z = z**2 + c #This sum must always be -2<z<2 for c to be a mandelbrot number
        if abs(z) > 2: #when c is not in the mandelbrot set:
            mdbNum = 0 #function returns a zero
            break
        i += 1
    if abs(z) < 2: # when c is in the mandelbrot set:
        mdbNum = 1 #function returns a one
    return mdbNum

#Fill the matrix indicating numbers that are in the mandelbrot set
for i in range(len(x)):
    for n in range(len(y)):
        c = complex(x[i],y[n]) #calculate imaginary number c from x and y coordinates
        mdbSet[i,n] = mdbCheck(c) #Set mdbSet[i,n] to a 0 or 1 based on whether c is in the mandelbrot set.

#Make figure of the mandelbrot set

import matplotlib.pyplot as plt
from matplotlib import cm

#Make density plot from mandelbrot data
plt.imshow(mdbSet, cmap = cm.gray_r, origin="lower")

#Label axis tick marks
ticks = ['-2', '-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2']
plt.xticks(np.linspace(0,1000,num=9,endpoint=True), ticks)
plt.yticks(np.linspace(0,1000,num=9,endpoint=True), ticks)

#Label axes
plt.xlabel("X Values")
plt.ylabel("Y Values")
plt.title("The Mandelbrot Set")

#save file as pdf
plt.savefig("./MandelbrotSet.pdf")
