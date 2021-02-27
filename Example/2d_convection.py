"""
Solving 2D linear convection equation; more details on https://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/07_Step_5.ipynb

Bear in mind that this file might take a while to run due to the use of double for-loop to iterate through all coordinate space.
"""



from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt


#Initialising parameters for simulation grid
nx = 201 #Number of elements in x direction
ny = 201 #Number of elements in y direction
nt = 100 #Number of timesteps
cx = 0.5 #x-component of velocity
cy = 0.7 #y-component of velocity
x_range = 4 #Range of x
y_range = 4 #Range of y
dx = x_range/nx #Cell size in x direction
dy = y_range/ny #Cell size in y direction
sigma = 0.2 #Parameter used to scale timestep to order of magnitude of cell size - keep as < 1
dt = sigma * dx

#Initialising axis
x = np.linspace(0,x_range,nx) #Initialising x-axis
y = np.linspace(0,y_range,ny) #Initialising y-axis

u = np.ones((ny,nx)) #Creating 2D array to representing coordinate space
un = np.ones((ny,nx)) #Creating placeholder solution array

#Setting initial conditions
#Location of the initial 2D tophat function
xa = 0.5
xb = 1.0
ya = 0.5
yb = 1.5

u[(int(ya/dx)):(int(yb/dx + 1)),(int(xa/dx)):(int(xb/dx + 1))] = 2 #Setting values between A and B to 2 in both axes

#Visualising initial condition
fig = plt.figure(figsize=(14,8),dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y,u[:])
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.show()

#Looping through timesteps to solve equation
for time in range(nt):
    un = u.copy() #Creating copy of solution in previous timestep
    for j in range(1,ny): #Looping through all y
        for i in range(1,nx): #Looping through all x
            u[j,i] = un[j,i] - cx*(dt/dx)*(un[j,i]-un[j,i-1]) - cy*(dt/dy)*(un[j,i]-un[j-1,i]) #Equation for method of finite differences
            
            #Applying boundary conditions: u=1 for x=0,4 and y=0,4
            u[0,:] = 1 #Boundary condition along first value of y for all x i.e. u=1 at y=0
            u[-1,:] = 1 #Boundary condition along last value of y for all x i.e. u=1 at y=4
            u[:,0] = 1 #Boundary condition along first value of x for all y i.e. u=1 at x=0
            u[:,-1] = 1 #Boundary condition along last value of x for all y i.e. u=1 at x=4

#Visualising solution
#Visualising initial condition
fig = plt.figure(figsize=(14,8),dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y,u[:])
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.show()