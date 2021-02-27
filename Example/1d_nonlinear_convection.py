import matplotlib.pyplot as plt
import numpy as np
import time, sys



#Initialising simulation grid
nx = 161 #Number of elements along x-axis
x_range = 10
dx = x_range/(nx-1) #Spacings of elements along x-axis
nt = 250 #Number of time steps
dt = 0.01 #Size of time step
c = 1 #Wavespeed


u = np.ones(nx) #Initialising solution
t_range = []
t_soln = [] #Array containing solution at all timesteps

#Setting initial conditions - Sort of a top-hat function used
a = 0.5
b = 0.8
u[(int(a/dx)):(int(b/dx + 1))] = 2 #Setting values between A and B to 2

#Visualising initial conditions
plt.plot(np.linspace(0,x_range,nx),u)
plt.show()

un = np.ones(nx) #Initialing placeholder array that should hold the solution

for i in range(nt): #Iterating through time steps
    un = u.copy()
    for i in range(1,nx): #Iterating through the u array (spatial solution)
        u[i] = un[i] - un[i]*(dt/dx)*(un[i]-un[i-1])




#Plotting solution at end of time iteration
plt.plot(np.linspace(0,x_range,nx),u)
plt.show()
