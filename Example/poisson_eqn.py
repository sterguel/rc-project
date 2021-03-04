"""
Solving Poisson's equation, which is just Laplace's equation with a inhomogenoeus term. More details at: https://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/13_Step_10.ipynb

"""
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt


#Initialising parameters for simulation grid
nx = 101 #Number of elements in x direction
ny = 101 #Number of elements in y direction
nt = 1000 #Number of timesteps
x_range = 2 #Range of x
y_range = 2 #Range of y
dx = x_range/nx #Cell size in x direction
dy = y_range/ny #Cell size in y direction
dt = 0.01


#Initialising axis
x = np.linspace(0,x_range,nx) #Initialising x-axis
y = np.linspace(0,y_range,ny) #Initialising y-axis

#Initialising solution array
u = np.zeros((ny,nx))
un = np.zeros((ny,nx))

#PVisualising initial conditions
fig = plt.figure(figsize=(14,8),dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y,u[:],cstride=1,cmap=cm.viridis) #Makes the plot look prettier
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.show()


#Initialising source term
b = np.zeros((ny,nx))
b[int(ny/2),int(nx/2)] = -1
#Source term consists of a constant value at a single point with 0 everywhere else.
#This is chosen as it would make Poisson's equation correspond to the potential of a point charge in an electric or gravitational field.
#Which is a well known example in physics, making it very suitable for testing the PDE solver.

for n in range(nt):
    un = u.copy()
    u[1:-1,1:-1] = (((dy**2) * (un[1:-1,2:] + un[1:-1,:-2])) + ((dx**2) * (un[2:,1:-1] + un[:-2,1:-1])) - (b[1:-1,1:-1] * (dy**2)))/(2 * (dx**2 + dy**2))

    #Setting boundaries to 0
    u[0,:] = 0
    u[-1,:] = 0
    u[:,0] = 0
    u[:,-1] = 0

#Plotting u after timesteps
fig = plt.figure(figsize=(14,8),dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y,u[:],cstride=1,cmap=cm.viridis) #Makes the plot look prettier
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.show()

