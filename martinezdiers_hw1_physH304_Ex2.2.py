#Exercise 2.2b
'''
#pseudocode:
-ask for the period of the satellite in minutes using input command.
-convert the period to seconds by multiplying the input number by 60.
-define constants G, M, and R that are used in the satellite altitude calculation.
-plug in constants and period in seconds into the satellite altitude formula to calculate the satellite altitude.
-print the altitude of a satellite with the desired period.
'''

#Import pi for calculation
from math import pi

#T is the desired period in minutes
Tinput = float(input("Enter the period of the satellite in minutes: "))
T = Tinput*60 #convert satellite period from minutes to seconds for altitude calculation

#set constants
G = 6.67e-11 #gravitational constant in m^3*kg^(-1)*s^(-2)
M = 5.97e24 #mass of the Earth in kg
R = 6371e3 #radius of the Earth in m

#calculate altitude of the satellite from the Earth's surface
h = (G*M*T**2/4/pi**2)**(1/3)-R

print("The altitude of a satellite with a period ", Tinput, " minutes is ", h, " meters.")

#Exercise 2.2c
'''
The altitude of a geosynchronous satellite is 35855910 meters.
The altitude of a satellite with a period of 90 minutes is 279321 meters.
The altitude of a satellite with a period of 45 minutes is -2181559 meters.
Because the altitude of a satellite with a period of 45 minutes is a negative number, it is not
possible to have a satellite orbit with this period. This is because the required radius between the
satellite and the center of the Earth is smaller than the radius of the Earth.
'''

#Exercise 2.2d
'''
The difference between the altitude of the satellite with a period of a day versus a sidereal day is 
82148 meters. There is a difference between the sidereal day and the full day because a sidereal day is just the time it takes for the 
earth to rotate while a full day is the time it takes for the earth to do a full rotation relative to the sun. Because of the 
translational motion of the earth relative to the sun, it takes longer for the earth to perform a full rotation relative to the sun rather 
than relative to its axis.
'''