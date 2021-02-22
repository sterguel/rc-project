import matplotlib.pyplot as plt
import numpy as np
import time, sys



#Initialising simulation grid
nx = 81 #Number of elements along x-axis
x_range = 4
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

un = np.ones(nx) #Initialise temporary array
for t in range(nt): #Iterating over time steps
    un = u.copy() #Creates copy of u at previous timestep
    for i in range(1,nx):
        u[i] = un[i] - c*(dt/dx)*(un[i]-un[i-1]) #Computing new value using value at previous timestep using finite
    t_soln.append(u) #Copying solution to array t_soln
    t_range.append(t) #Copying time after start to array t_range

#Plotting solution after final timestep
plt.plot(np.linspace(0,x_range,nx),u)
plt.show()
