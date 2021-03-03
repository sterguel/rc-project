"""
Solving 2D nonlinear convection equation; more details on https://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/08_Step_6.ipynb

This example involves solving 2 coupled differential equations.

Bear in mind that this file might take a while to run due to the use of double for-loop to iterate through all coordinate space.
"""



from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt


#Initialising parameters for simulation grid
nx = 101 #Number of elements in x direction
ny = 101 #Number of elements in y direction
nt = 100 #Number of timesteps
cx = 2 #x-component of velocity
cy = 3 #y-component of velocity
x_range = 2 #Range of x
y_range = 2 #Range of y
dx = x_range/nx #Cell size in x direction
dy = y_range/ny #Cell size in y direction
sigma = 0.2 #Parameter used to scale timestep to order of magnitude of cell size - keep as < 1
dt = sigma * dx

#Initialising axis
x = np.linspace(0,x_range,nx) #Initialising x-axis
y = np.linspace(0,y_range,ny) #Initialising y-axis

u = np.ones((ny,nx)) #Creating 2D array to representing coordinate space
v = np.ones((ny,nx))

un = np.ones((ny,nx)) #Creating placeholder solution array
vn = np.ones((ny,nx))

#Setting initial conditions
#Initial conditions for u
uxa = 0.5
uxb = 1.0
uya = 0.5
uyb = 1.5

u[(int(uya/dx)):(int(uyb/dx + 1)),(int(uxa/dx)):(int(uxb/dx + 1))] = 2 #Setting values between A and B to 2 in both axes

#Initial conditions for v
vxa = 0.5
vxb = 1.0
vya = 0.5
vyb = 1.5

v[(int(vya/dx)):(int(vyb/dx + 1)),(int(vxa/dx)):(int(vxb/dx + 1))] = 1.5 #Setting values between A and B to 2 in both axes

#Visualising initial conditions for u
fig = plt.figure(figsize=(14,8),dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y,u[:])
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.show()

#Visualising initial conditions for v
fig = plt.figure(figsize=(14,8),dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y,v[:])
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.show()


for n in range(nt):
    un = u.copy()
    vn = v.copy()
    for j in range(1,ny): #Looping through y
        for i in range(1,nx): #Looping through x
            #Solving coupled DEs
            u[j,i] = u[j,i] - u[j,i]*(dt/dx)*(u[j,i]-u[j,i-1]) - v[j,i]*(dt/dy)*(u[j,i]-u[j-1,i])
            v[j,i] = v[j,i] = u[j,i]*(dt/dx)*(v[j,i]-v[j,i-1]) - v[j,i]*(dt/dy)*(v[j,i]-v[j-1,i])

            #Setting boundary conditions
            #Setting u,v = 1 at edge as boundary condition
            u[0,:] = 1 #Setting u=1 at y=0
            u[-1,:] = 1 #Setting u=1 at y=2
            u[:,0] = 1 #Setting u=1 at x=0
            u[:,-1] = 1 #Setting u=1 at x=2

            #Doing the same for v
            v[0,:] = 1
            v[-1,:] = 1
            v[:,0] = 1
            v[:,-1] = 1


#Plotting u
fig = plt.figure(figsize=(14,8),dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y,u[:])
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.show()

#Plotting v
fig = plt.figure(figsize=(14,8),dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y,v[:])
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.show()

