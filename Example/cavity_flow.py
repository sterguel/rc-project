"""
Solving the Navier-Stokes equation for cavity flow i.e. container with boundaries on 3 sides and a constant velocity on one.

In this example, the velocity at the top of the container is constant at 1 in the x-direction.

More details on https://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/14_Step_11.ipynb
"""

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt


#Initialising parameters for simulation grid
nx = 41 #Number of elements in x direction
ny = 41 #Number of elements in y direction
nt = 500 #Number of timesteps
x_range = 2 #Range of x
y_range = 2 #Range of y
dx = x_range/nx #Cell size in x direction
dy = y_range/ny #Cell size in y direction
dt = 0.001

#Defining fluid properties
rho = 1 #Density of fluid
nu = 0.1 #Viscosity

#Initialising arrays
u = np.zeros((ny,nx)) #Velocity in x-direction
v = np.zeros((ny,nx)) #Velocity in y-direction
b = np.zeros((ny,nx)) #Placeholder to make life easier with Poisson's equation
p = np.zeros((ny,nx)) #Pressure in fluid

#Setting initial conditions for problem
u[:,-1] = 0 #x-velocity = 0 @ x=0
u[:,0] = 0 #x-velocity = 0 @ x=xmax
u[0,:] = 0 #x-velocity = 0 @ y=0
u[-1,:] = 1 #x-velocity = 1 @ y=ymax. Note the order that this overwrites initial condition at x=0,xmax
#Doing the same for y-velocity
v[:,-1] = 0
v[:,0] = 0
v[0,:] = 0
v[-1,:] = 0

def get_b(b,u,v):
    #Function for computing inhomogeneous term of the Poisson pressure equation
    b[1:-1, 1:-1] = rho * (
        ((1/dt) * ( ((u[1:-1,2:] - u[1:-1,0:-2])/(2*dx)) + ((v[2:,1:-1] - v[0:-2,1:-1])/(2*dy)) ))
        - ( ((u[1:-1,2:] - u[1:-1,0:-2])/(2*dx))**2)
        - 2*( ((u[2:,1:-1] - u[0:-2,1:-1])/(2*dy)) * ((v[1:-1,2:] - v[1:-1,0:-2])/(2*dx)) )
        - ( ((v[2:,1:-1] - v[0:-2,1:-1])/(2*dy))**2 )
        )
    
    return b





def get_pressure(p,u,v,b):
    #Function to solve the Poisson pressure equation
    pn = np.empty_like(p)
    pn = p.copy()
    nit = 50 #Pseudotime to advance equation
    for q in range(nit):
        pn = p.copy()
        p[1:-1,1:-1] = (
            ( ( ((pn[1:-1,2:] + pn[1:-1,0:-2])*(dy**2)) +  ((pn[2:,1:-1] + pn[0:-2,1:-1])*(dx**2)) ) / (2*((dx**2) + (dy**2))) )
            - ( ( ((dx**2)*(dy**2)) / (2*((dx**2) + (dy**2))) ) * b[1:-1,1:-1] ) #The b function as defined above makes my life a lot easier in this line
            )

        #Setting boundary condition for pressures
        p[:,-1] = p[:,-2] #dp/dx=0 at x=xmax
        p[:,0] = p[:,1] #dp/dx=0 at x=0
        p[0,:] = p[1,:] #dp/dy=0 at y=0
        p[-1,:] = 0 #p=0 at top of container
    
    return p



#Solving cavity flow problem
#Initialising placeholder array
un = np.empty_like(u)
vn = np.empty_like(v)

#Actually solving the problem
for n in range(nt): #Iterating through time
    un = u.copy()
    vn = v.copy()
    #Calculating pressure throughout fluid
    b = get_b(b,u,v)
    p = get_pressure(p,u,v,b)


    #FDM to solve for velocity components
    #x-component of velocity
    u[1:-1,1:-1] = (un[1:-1,1:-1] - (un[1:-1,1:-1]*(dt/dx)*(un[1:-1,1:-1] - un[1:-1,0:-2]))
    - (vn[1:-1,1:-1]*(dt/dy)*(un[1:-1,1:-1] - un[0:-2,1:-1]))
    - (((dt)/(rho*2*dx)) * (p[1:-1,2:] - p[1:-1,0:-2]))
    + nu*(((dt/(dx**2))*(un[1:-1,2:] - (2*un[1:-1,1:-1]) + un[1:-1,0:-2])) + ((dt/(dy**2))*(un[2:,1:-1] - (2*un[1:-1,1:-1]) + un[0:-2,1:-1]))))
    #y-component of velocity
    v[1:-1,1:-1] = (vn[1:-1,1:-1] - (vn[1:-1,1:-1]*(dt/dx)*(vn[1:-1,1:-1] - un[1:-1,0:-2]))
    - (vn[1:-1,1:-1]*(dt/dy)*(vn[1:-1,1:-1] - vn[0:-2,1:-1]))
    - (((dt)/(rho*2*dy)) * (p[2:,1:-1] - p[0:-2,1:-1]))
    + nu*(((dt/(dx**2))*(vn[1:-1,2:] - (2*vn[1:-1,1:-1]) + vn[1:-1,0:-2])) + ((dt/(dy**2))*(vn[2:,1:-1] - (2*vn[1:-1,1:-1]) + vn[0:-2,1:-1]))))

    #Setting boundary condition for velocities; set to 0 on all edges of container except for top.
    u[0, :]  = 0
    u[:, 0]  = 0
    u[:, -1] = 0
    u[-1, :] = 1    # set velocity on cavity lid equal to 1
    v[0, :]  = 0
    v[-1, :] = 0
    v[:, 0]  = 0
    v[:, -1] = 0



#Plots velocity vector throughout fluid
x = np.linspace(0,x_range,nx)
y = np.linspace(0,y_range,ny)
X,Y = np.meshgrid(x,y)

fig = plt.figure(figsize=(11,7), dpi=100)
plt.quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2]) 
plt.xlabel('X')
plt.ylabel('Y');
plt.show()