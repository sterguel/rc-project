"""
Solving the Navier-Stokes equation for channel flow driven by a constant pressure.

This example is the cavity flow file copied and modified as much of the code for cavity flow carries over to channel flow with a few modifications.

More details on https://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/15_Step_12.ipynb

To change the shape of the flow channel, there are 2 things that need to be changed:
    - The pressure boundary conditions at the walls of the container (See line 111)
    - The velocity boundary conditions at the walls of the container (See line 178)
"""

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt

#Data save filename and filepath
fname = 'data_unblocked'
output_path = 'Output\\' + fname 

#Initialising parameters for simulation grid
nx = 501 #Number of elements in x direction
ny = 201 #Number of elements in y direction
nt = 500 #Number of timesteps
x_range = 100 #Range of x (Length of modelled section in mm)
y_range = 12 #Range of y (Diameter of vessel in mm)

dx = x_range/nx #Cell size in x direction
dy = y_range/ny

#dy = y_range/ny #Cell size in y direction

dt = 0.000001

#Defining fluid properties
rho = 1.0565e-3 #Density of fluid (g/mm^3)
nu = 2.57 #Viscosity

#Simulation parameters
Fx = 3.629e5 #Force in x-direction (defined in +x direction) (Defined in g*mm/s^2)



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

    #Defining pressure boundary conditions along x-axis
    
    #Periodic boundary conditions for pressure at x=xmax
    b[1:-1,-1] = rho * (
        ((1/dt) * ( ((u[1:-1,0] - u[1:-1,-2])/(2*dx)) + ((v[2:,-1] - v[0:-2,-1])/(2*dy)) ))
        - ( ((u[1:-1,0] - u[1:-1,-2])/(2*dx))**2)
        - 2*( ((u[2:,-1] - u[0:-2,-1])/(2*dy)) * ((v[1:-1,0] - v[1:-1,-2])/(2*dx)) )
        - ( ((v[2:,-1] - v[0:-2,-1])/(2*dy))**2 )
        )

    
    #Periodic boundary conditions for pressure at x=0
    b[1:-1,0] = rho * (
        ((1/dt) * ( ((u[1:-1,1] - u[1:-1,-1])/(2*dx)) + ((v[2:,0] - v[0:-2,0])/(2*dy)) ))
        - ( ((u[1:-1,1] - u[1:-1,-1])/(2*dx))**2)
        - 2*( ((u[2:,0] - u[0:-2,0])/(2*dy)) * ((v[1:-1,1] - v[1:-1,-1])/(2*dx)) )
        - ( ((v[2:,0] - v[0:-2,0])/(2*dy))**2 )
        )
    return b





def get_pressure(p,u,v,b):
    #Function to solve the Poisson pressure equation
    pn = np.empty_like(p)
    pn = p.copy()
    nit = 50 #Pseudotime to advance equation
    for q in range(nit):
        pn = p.copy()
        p[1:-1,1:-1] = ( #Determining pressure
            ( ( ((pn[1:-1,2:] + pn[1:-1,0:-2])*(dy**2)) +  ((pn[2:,1:-1] + pn[0:-2,1:-1])*(dx**2)) ) / (2*((dx**2) + (dy**2))) )
            - ( ( ((dx**2)*(dy**2)) / (2*((dx**2) + (dy**2))) ) * b[1:-1,1:-1] ) #The b function as defined above makes my life a lot easier in this line
            )
        
        #Periodic pressure BC at x=xmax
        p[1:-1,-1] = (
            ( ( ((pn[1:-1,0] + pn[1:-1,-2])*(dy**2)) +  ((pn[2:,-1] + pn[0:-2,-1])*(dx**2)) ) / (2*((dx**2) + (dy**2))) )
            - ( ( ((dx**2)*(dy**2)) / (2*((dx**2) + (dy**2))) ) * b[1:-1,-1] )
        )

        #Periodic pressure BC at x=0
        p[1:-1,0] = (
            ( ( ((pn[1:-1,1] + pn[1:-1,-1])*(dy**2)) +  ((pn[2:,0] + pn[0:-2,0])*(dx**2)) ) / (2*((dx**2) + (dy**2))) )
            - ( ( ((dx**2)*(dy**2)) / (2*((dx**2) + (dy**2))) ) * b[1:-1,0] )
            )

        #Wall boundary condition for pressures dp/dy = 0 at walls y=0,ymax
        #Change this to change the shape of channel
        p[-1,:] = p[-2,:] #dp/dy=0 at y=ymax
        p[0,:] = p[1,:] #dp/dy=0 at y=0
    
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
    + nu*(((dt/(dx**2))*(un[1:-1,2:] - (2*un[1:-1,1:-1]) + un[1:-1,0:-2])) + ((dt/(dy**2))*(un[2:,1:-1] - (2*un[1:-1,1:-1]) + un[0:-2,1:-1])))
    + dt*Fx) #Last line added a force term to represent contribution to the momentum from force being applied.
    #y-component of velocity
    v[1:-1,1:-1] = (vn[1:-1,1:-1] - (vn[1:-1,1:-1]*(dt/dx)*(vn[1:-1,1:-1] - un[1:-1,0:-2]))
    - (vn[1:-1,1:-1]*(dt/dy)*(vn[1:-1,1:-1] - vn[0:-2,1:-1]))
    - (((dt)/(rho*2*dy)) * (p[2:,1:-1] - p[0:-2,1:-1]))
    + nu*(((dt/(dx**2))*(vn[1:-1,2:] - (2*vn[1:-1,1:-1]) + vn[1:-1,0:-2])) + ((dt/(dy**2))*(vn[2:,1:-1] - (2*vn[1:-1,1:-1]) + vn[0:-2,1:-1]))))

    #Setting periodic BCs for velocities
    #Note that the periodic BCs for velocities is commented out.
    #This is because I don't see why velocity should wrap around back to the start for a small, unique section of a blood vessel. Feel free to change if needed.

    #Periodic BC u @ x=xmax
    u[1:-1,-1] = (un[1:-1,-1] - (un[1:-1,-1]*(dt/dx)*(un[1:-1,-1] - un[1:-1,-2]))
    - (vn[1:-1,-1]*(dt/dy)*(un[1:-1,-1] - un[0:-2,-1]))
    - (((dt)/(rho*2*dx)) * (p[1:-1,0] - p[1:-1,-2]))
    + nu*(((dt/(dx**2))*(un[1:-1,0] - (2*un[1:-1,-1]) + un[1:-1,-2])) + ((dt/(dy**2))*(un[2:,-1] - (2*un[1:-1,-1]) + un[0:-2,-1])))
    + dt*Fx)

    #Periodic BC u @x=0
    u[1:-1,0] = (un[1:-1,0] - (un[1:-1,0]*(dt/dx)*(un[1:-1,0] - un[1:-1,-1]))
    - (vn[1:-1,0]*(dt/dy)*(un[1:-1,0] - un[0:-2,0]))
    - (((dt)/(rho*2*dx)) * (p[1:-1,1] - p[1:-1,-1]))
    + nu*(((dt/(dx**2))*(un[1:-1,1] - (2*un[1:-1,0]) + un[1:-1,-1])) + ((dt/(dy**2))*(un[2:,0] - (2*un[1:-1,0]) + un[0:-2,0])))
    + dt*Fx)
    #Periodic BC v @x=xmax
    v[1:-1,-1] = (vn[1:-1,-1] - (vn[1:-1,-1]*(dt/dx)*(vn[1:-1,-1] - un[1:-1,-2]))
    - (vn[1:-1,-1]*(dt/dy)*(vn[1:-1,-1] - vn[0:-2,-1]))
    - (((dt)/(rho*2*dy)) * (p[2:,-1] - p[0:-2,-1]))
    + nu*(((dt/(dx**2))*(vn[1:-1,0] - (2*vn[1:-1,-1]) + vn[1:-1,-2])) + ((dt/(dy**2))*(vn[2:,-1] - (2*vn[1:-1,-1]) + vn[0:-2,-1]))))

    #Periodic BC v @x=0
    v[1:-1,0] = (vn[1:-1,0] - (vn[1:-1,0]*(dt/dx)*(vn[1:-1,0] - un[1:-1,-1]))
    - (vn[1:-1,0]*(dt/dy)*(vn[1:-1,0] - vn[0:-2,0]))
    - (((dt)/(rho*2*dy)) * (p[2:,0] - p[0:-2,0]))
    + nu*(((dt/(dx**2))*(vn[1:-1,1] - (2*vn[1:-1,0]) + vn[1:-1,-1])) + ((dt/(dy**2))*(vn[2:,0] - (2*vn[1:-1,0]) + vn[0:-2,0]))))

    #Wall boundary conditions - change this to change shape of channel
    u[0,:] = 0 #u=0 @ y=0
    u[-1,:] = 0 #u=0 @ y=ymax
    v[0,:] = 0 #v=0 @ y=0
    v[-1,:] = 0 #v=0 @ y=ymax


#Plots velocity vector throughout fluid
x = np.linspace(0,x_range,nx)
y = np.linspace(0,y_range,ny)
X,Y = np.meshgrid(x,y)

fig = plt.figure(figsize=(11,7), dpi=100)
plt.quiver(X[::2, ::10], Y[::2, ::10], u[::2, ::10], v[::2, ::10]) 
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

np.savez(file=output_path,x_vel=u,y_vel=v,pressure=p,x_coord=X,y_coord=Y)