"""
For more details, refer to https://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/05_Step_4.ipynb

Burger's equation describes speed of a fluid at positions along a thin, ideal pipe at different times.

Initial conditions is one period of the sawtooth wave, but function used to represent it is not periodic. As such, a periodic boundary condition will be applied.
"""

import numpy as np
import sympy
import matplotlib.pyplot as plt
from sympy import init_printing
from sympy.utilities.lambdify import lambdify
init_printing(use_latex=True)

#Setting up initial conditions - refer to webpage for analytical form
x, nu, t = sympy.symbols('x nu t')

#Defining function used in initial condition
phi = (sympy.exp(-(x - 4 * t)**2 / (4 * nu * (t + 1))) +
       sympy.exp(-(x - 4 * t - 2 * sympy.pi)**2 / (4 * nu * (t + 1))))

#Taking deriative
phiprime = phi.diff(x)

#Turning intital condition into a function
u = (-2*nu*(phiprime/phi)) + 4
ufunc = lambdify((t,x,nu),u) #Turning u into function with parameters t,x,nu

#Initialising simulation grid
nx = 251 #Number of elements along x-axis
x_range = 2*np.pi
dx = x_range/(nx-1) #Spacings of elements along x-axis
nt = 500 #Number of time steps
nu = 0.07 #Viscosity
dt = 0.0002 #Timestep size

#Be careful when choosing the simulation parameters; entire simulation breaks if time or space steps are too big (as one should expect when using FDM to solve DEs)

#Turning initial condition into an array containing its values at discrete points
x = np.linspace(0,x_range,nx)
un = np.empty(nx)
t = 0 #Initial time

u = ufunc(t,x,nu) #Calling function on a np array containing list of x positions to create array representing function (i.e. initial state of system)

#Plotting initial conditions for the sake of visualisation
plt.plot(x,u,label='Initial conditions')
#plt.show()

#Solving by method of finite differences
for time in range(nt):
    un = u.copy()
    for i in range(1,nx-1):
        #Computing next cell using method of finite differences
        u[i] = un[i] - u[i]*(dt/dx)*(un[i]-un[i-1]) + nu*(dt/(dx**2))*(un[i+1]-2*un[i]+un[i-1])
    #Applying boundary condition - Cell at right hand of frame wraps around back to the front (i.e. start of the solution)
    u[0] = un[0] - u[0]*(dt/dx)*(un[0]-un[-2]) + nu*(dt/(dx**2))*(un[1]-2*un[0]+un[-2]) #Index -2 of solution list is the last cell that is computed by the loop through space.
    u[-1] = u[0]


#Plotting solution
plt.plot(x,u,label='Time-dependent solution')
plt.legend()
plt.show()

