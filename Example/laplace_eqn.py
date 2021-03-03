"""
Solving Laplace's equation using finite-difference methods, details outline in: https://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/12_Step_9.ipynb

From this example onwards, array operations instead of nested loops will be used to copmute the solution.
"""
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt


#Initialising parameters for simulation grid
nx = 101 #Number of elements in x direction
ny = 101 #Number of elements in y direction
nt = 5000 #Number of timesteps
x_range = 2 #Range of x
y_range = 2 #Range of y
dx = x_range/nx #Cell size in x direction
dy = y_range/ny #Cell size in y direction
dt = 0.01

#Initialising axis
x = np.linspace(0,x_range,nx) #Initialising x-axis
y = np.linspace(0,y_range,ny) #Initialising y-axis

u = np.zeros((ny,nx)) #Creating 2D array to representing coordinate space
un = np.zeros((ny,nx))

#Setting initial conditions
u[0,:] = 5 #Setting all values of u @ y=0 to 5 for initial condition

#Visualising initial conditions for u
fig = plt.figure(figsize=(14,8),dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y,u[:])
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.show()

for t in range(nt): #Iterating over timesteps
    un = u.copy() #Copying u into temporary array
    u[1:-1,1:-1] = (((dy**2) * (u[1:-1,2:] + u[1:-1,0:-2])) + ((dx**2) * (u[2:,1:-1] + u[0:-2,1:-1]))) / (2 * ((dx**2) + (dy**2))) #Solving using array method

    #Setting boundary conditions - u=0 at boundaries except for at source
    u[0,:] = 5 #u=5 @ ymin; treating as a constant source
    u[-1,:] = 0 #u=0 @ ymax
    u[:,0] = 0 #u=0 @ xmin
    u[:,-1] = 0 #u=0 @xmax

#Plotting u after timesteps
fig = plt.figure(figsize=(14,8),dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y,u[:],cstride=1,cmap=cm.viridis) #Makes the plot look prettier
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.show()