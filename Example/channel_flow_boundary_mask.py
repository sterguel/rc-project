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


#Initialising parameters for simulation grid
nx = 150 #Number of elements in x direction
ny = 101 #Number of elements in y direction
nt = 500 #Number of timesteps
x_range = 8 #Range of x
y_range = 2 #Range of y

dx = x_range/nx #Cell size in x direction

#dy = y_range/ny #Cell size in y direction

dt = 0.0001

#Defining fluid properties
rho = 1 #Density of fluid
nu = 0.1 #Viscosity

#Simulation parameters
Fx = 1 #Force in x-direction (defined in +x direction)



#Initialising arrays
u = np.zeros((ny,nx),dtype=np.float64) #Velocity in x-direction
v = np.zeros((ny,nx),dtype=np.float64) #Velocity in y-direction
b = np.zeros((ny,nx),dtype=np.float64) #Placeholder to make life easier with Poisson's equation
p = np.zeros((ny,nx),dtype=np.float64) #Pressure in fluid

#Setting initial conditions for problem
u[:,-1] = 0 #x-velocity = 0 @ x=0
u[:,0] = 0 #x-velocity = 0 @ x=xmax
u[0,:] = 0 #x-velocity = 0 @ z=0
u[-1,:] = 2 #x-velocity = 1 @ z=ymax. Note the order that this overwrites initial condition at x=0,xmax

#Doing the same for y-velocity
v[:,-1] = 0
v[:,0] = 0
v[0,:] = 0
v[-1,:] = 0

#calculate actual grid
def f(x):
    #issues with rounding errors even though this function is symmetric?
    return round(np.exp(-(x-4)**2),2)
    #return 0
def f_d(x):
    return round(-2 * (x-4)* np.exp(-(x-4)**2),2)
#more y grid points are needed for a non symmetric function

fvs = list(set([f(dx*i) for i in range(nx)]))
if max(fvs) < y_range:
    fvs += list(np.linspace(y_range, max(fvs), ny - len(fvs), False))
fvs.sort()
finvmap = [fvs.index(f(i*dx)) for i in range(nx)]

inside_mask = np.zeros((ny,nx,))
ybottom_mask = np.zeros((ny,nx,))
dxu_mask = np.zeros((ny,nx,))
dxl_mask = np.zeros((ny,nx,))
dyu_mask = np.zeros((ny,nx,))
dyl_mask = np.zeros((ny,nx,))
for j in range(ny):
    for i in range(nx):
        def interior(i,j):
            return (j > finvmap[i]) and (j < ny - 1) and (i > 0) and (i<nx - 1)
        inside_mask[j,i] = interior(i,j)
        ybottom_mask[j,i] = j == finvmap[i]
        dxu_mask[j,i] = interior(i-1,j)
        dxl_mask[j,i] = interior(i+1,j)
        dyu_mask = interior(i,j-1)
        dyl_mask = interior(i,j+1)

dys =  [fvs[i] - fvs[i-1] for i in range(1, len(fvs))]
dys2 = [fvs[i+1] - fvs[i-1] for i in range(1,len(fvs) - 1) ]

dy = np.reshape(nx *dys, ( nx,len(dys),)).T
dy2 = np.reshape(nx *dys2, ( nx,len(dys2),)).T
ocf = 2/(dy[1:,1:-1] * dy[:-1,1:-1] * dy2[:,1:-1])
def get_b(b,u,v):
    #Function for computing inhomogeneous term of the Poisson pressure equation
    b[inside_mask] = rho * (
        ((1/dt) * ( ((u[dxu_mask] - u[dxl_mask])/(2*dx)) + ((v[dyu_mask] - v[dyl_mask])/dy2[inside_mask[1:-1,:]]) ))
        -  ((u[dxu_mask] - u[dxl_mask])/(2*dx))**2
        - 2*( ((u[dyu_mask] - u[dyl_mask])/dy2[inside_mask[1:-1,:]]) * ((v[dxu_mask] - v[dxl_mask])/(2*dx)) )
        - ((v[dyu_mask] - v[dyl_mask])/dy2[inside_mask[1:-1,:]])**2 
        )

    #Defining pressure boundary conditions along x-axis
    
    #Periodic boundary conditions for pressure at x=xmax
    b[1:-1,-1] = rho * (
        ((1/dt) * ( ((u[1:-1,0] - u[1:-1,-2])/(2*dx)) + ((v[:,-1][dyu_mask[:,-1]] - v[0:-2,-1])/dy2[:,-1]) ))
        - ( ((u[1:-1,0] - u[1:-1,-2])/(2*dx))**2)
        - 2*( ((u[2:,-1] - u[0:-2,-1])/dy2[:,-1]) * ((v[1:-1,0] - v[1:-1,-2])/(2*dx)) )
        - ( ((v[2:,-1] - v[0:-2,-1])/dy2[:,-1])**2 )
        )

    
    #Periodic boundary conditions for pressure at x=0
    b[1:-1,0] = rho * (
        ((1/dt) * ( ((u[1:-1,1] - u[1:-1,-1])/(2*dx)) + ((v[2:,0] - v[0:-2,0])/dy2[:,0]) ))
        - ( ((u[1:-1,1] - u[1:-1,-1])/(2*dx))**2)
        - 2*( ((u[2:,0] - u[0:-2,0])/dy2[:,0]) * ((v[1:-1,1] - v[1:-1,-1])/(2*dx)) )
        - ( ((v[2:,0] - v[0:-2,0])/dy2[:,0])**2 )
        )
    return b





def get_pressure(p,u,v,b):
    #Function to solve the Poisson pressure equation
    pn = np.empty_like(p)
    pn = p.copy()
    nit = 50 #Pseudotime to advance equation
    da = 2/(dx**2) + 2/(dy[1:,1:-1] * dy[:-1, 1:-1])
    da = 1/da
    for q in range(nit):
        pn = p.copy()
        #p[1:-1,1:-1] = ( #Determining pressure
        #    ( ( ((pn[1:-1,2:] + pn[1:-1,0:-2])*(dy**2)) +  ((pn[2:,1:-1] + pn[0:-2,1:-1])*(dx**2)) ) / (2*((dx**2) + (dy**2))) )
        #    - ( ( ((dx**2)*(dy**2)) / (2*((dx**2) + (dy**2))) ) * b[1:-1,1:-1] ) #The b function as defined above makes my life a lot easier in this line
        #    )
        
        p[1:-1,1:-1] = da * (-b[1:-1,1:-1] + (pn[1:-1,2:] + pn[1:-1,0:-2])/dx**2 + ocf * (dy[ :-1 , 1:-1] * pn[2:,1:-1] + dy[ 1: , 1:-1] * pn[:-2, 1:-1] ) )


        #Periodic pressure BC at x=xmax
        #p[1:-1,-1] = (
        #    ( ( ((pn[1:-1,0] + pn[1:-1,-2])*(dy**2)) +  ((pn[2:,-1] + pn[0:-2,-1])*(dx**2)) ) / (2*((dx**2) + (dy**2))) )
        #    - ( ( ((dx**2)*(dy**2)) / (2*((dx**2) + (dy**2))) ) * b[1:-1,-1] )
        #)
        p[1:-1,-1] = da[:,0] * (-b[1:-1,-1] + (pn[1:-1,0] + pn[1:-1,-2])/dx**2 + ocf[:,0] * (dy[ :-1 , 0] * pn[2:,-1] + dy[ 1: ,0 ] * pn[:-2, -1] ) )

        #Periodic pressure BC at x=0
        #p[1:-1,0] = (
        #    ( ( ((pn[1:-1,1] + pn[1:-1,-1])*(dy**2)) +  ((pn[2:,0] + pn[0:-2,0])*(dx**2)) ) / (2*((dx**2) + (dy**2))) )
        #    - ( ( ((dx**2)*(dy**2)) / (2*((dx**2) + (dy**2))) ) * b[1:-1,0] )
        #    )
        p[1:-1,0] = da[:,0] * (-b[1:-1,0] + (pn[1:-1,1] + pn[1:-1,-1])/dx**2 + ocf[:,0] * (dy[ :-1 , 0] * pn[2:,0] + dy[ 1: ,0 ] * pn[:-2, 0] ) )

        
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
    #print(n)
    un = u.copy()
    vn = v.copy()
    #Calculating pressure throughout fluid
    b = get_b(b,u,v)
    p = get_pressure(p,u,v,b)


    #FDM to solve for velocity components
    #x-component of velocity
    u[1:-1,1:-1] = (un[1:-1,1:-1] - (un[1:-1,1:-1]*(dt/dx)*(un[1:-1,1:-1] - un[1:-1,0:-2]))
    - (vn[1:-1,1:-1]*(dt/dy[:-1,1:-1])*(un[1:-1,1:-1] - un[0:-2,1:-1]))
    - (dt/(rho*2*dx) * (p[1:-1,2:] - p[1:-1,0:-2]))
    + nu*(((dt/dx**2)*(un[1:-1,2:] - (2*un[1:-1,1:-1]) + un[1:-1,0:-2])) 
    + ( ocf * dt * ( dy[:-1,1:-1]*un[2:,1:-1] - (dy[1:,1:-1] + dy[:-1,1:-1]) *un[1:-1,1:-1] + dy[1:,1:-1] * un[0:-2,1:-1])  ))
    + dt*Fx) #Last line added a force term to represent contribution to the momentum from force being applied.
    #y-component of velocity
    v[1:-1,1:-1] = (vn[1:-1,1:-1] - (vn[1:-1,1:-1]*(dt/dx)*(vn[1:-1,1:-1] - un[1:-1,0:-2]))
    - (vn[1:-1,1:-1]*(dt/dy[:-1,1:-1])*(vn[1:-1,1:-1] - vn[0:-2,1:-1]))
    - (((dt)/(rho*dy2[:,1:-1])) * (p[2:,1:-1] - p[0:-2,1:-1]))
    + nu*(((dt/(dx**2))*(vn[1:-1,2:] - (2*vn[1:-1,1:-1]) + vn[1:-1,0:-2])) 
    +( ocf * dt * ( dy[:-1,1:-1]*vn[2:,1:-1] - (dy[1:,1:-1] + dy[:-1,1:-1]) *vn[1:-1,1:-1] + dy[1:,1:-1] * vn[0:-2,1:-1])  )) )


    #Wall boundary conditions - change this to change shape of channel
    #for i in range(nx):
    #    u[:finvmap[i],i] = 0
    #    v[:finvmap[i],i] = 0
    u[-1,:] = 0 #u=0 @ y=ymax
    v[-1,:] = 0 #v=0 @ y=ymax

    u[0,:] = 0 #u=0 @ y=0
    v[0,:] = 0 #v=0 @ y=0


#Plots velocity vector throughout fluid
x = np.linspace(0,x_range,nx)
y = fvs
X,Y = np.meshgrid(x,y)

fig = plt.figure(figsize=(11,7), dpi=100)
plt.quiver(X[::2, ::10], Y[::2, ::10], u[::2, ::10], v[::2, ::10]) 
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
