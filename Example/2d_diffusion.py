"""
Solving 2D diffusion equation; details on https://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/09_Step_7.ipynb

Bear in mind that this file might take a while to run due to the use of double for-loop to iterate through all coordinate space.
"""
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt


#Initialising parameters for simulation grid
nx = 101 #Number of elements in x direction
ny = 101 #Number of elements in y direction
nt = 250 #Number of timesteps
nu = 2 #Diffusivity i.e. how well matter diffuses
x_range = 2 #Range of x
y_range = 2 #Range of y
dx = x_range/nx #Cell size in x direction
dy = y_range/ny #Cell size in y direction
sigma = 0.2 #Parameter used to scale timestep to order of magnitude of cell size - keep as < 1
dt = sigma * dx * dy

#Initialising axis
x = np.linspace(0,x_range,nx) #Initialising x-axis
y = np.linspace(0,y_range,ny) #Initialising y-axis

u = np.ones((ny,nx)) #Creating 2D array to representing coordinate space
un = np.ones((ny,nx))

#Setting initial conditions
#Initial conditions for u
uxa = 0.7
uxb = 0.9
uya = 0.7
uyb = 0.9

u[(int(uya/dx)):(int(uyb/dx + 1)),(int(uxa/dx)):(int(uxb/dx + 1))] = 8

#Visualising initial conditions for u
fig = plt.figure(figsize=(14,8),dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y,u[:])
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.show()

for t in range(nt):
    un = u.copy()
    for j in range(1,ny-1):
        for i in range(1,nx-1):
            u[j,i] = un[j,i] + nu*(dt/(dx**2))*(u[j,i+1]-(2*u[j,i])+u[j,i-1]) + nu*(dt/(dy**2))*(u[j+1,i]-(2*u[j,i])+u[j-1,i])

            #Setting boundary conditions
            #Setting u,v = 1 at edge as boundary condition
            u[0,:] = 1 #Setting u=1 at y=0
            u[-1,:] = 1 #Setting u=1 at y=2
            u[:,0] = 1 #Setting u=1 at x=0
            u[:,-1] = 1 #Setting u=1 at x=2

#Plotting u
fig = plt.figure(figsize=(14,8),dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y,u[:])
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.show()