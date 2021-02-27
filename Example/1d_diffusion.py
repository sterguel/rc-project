"""
For more details, refer to https://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/04_Step_3.ipynb

In this script, a Taylor expansion of the solution and forward and backward diffference of the first derivative is used to discretize the second derivative.
"""
import matplotlib.pyplot as plt
import numpy as np
import time, sys



#Initialising simulation grid
nx = 251 #Number of elements along x-axis
x_range = 10
dx = x_range/(nx-1) #Spacings of elements along x-axis
nt = 2500 #Number of time steps
nu = 0.3 #Viscosity
sigma = 0.2
dt = sigma * dx**2 / nu

#Setting initial conditions - Sort of a top-hat function used
a = 4
b = 6
u = np.ones(nx)
u[(int(a/dx)):(int(b/dx + 1))] = 2 #Setting values between A and B to 2

#Visualising initial conditions
plt.plot(np.linspace(0, x_range, nx), u)
plt.show()


un = np.ones(nx) #Usual placeholder array

for t in range(nt): #Time iteration
    un = u.copy()
    for i in range(1,nx-1):
        u[i] = un[i] + nu*(dt/(dx**2))*(un[i+1] - 2*un[i] + un[i-1]) #Using discretization equation; see linked page for more details.


plt.plot(np.linspace(0,x_range,nx),u)
plt.show()