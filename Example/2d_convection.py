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
yb = 1.0

u[(int(xa/dx)):(int(xb/dx + 1)),(int(ya/dx)):(int(yb/dx + 1))] = 2 #Setting values between A and B to 2 in both axes

#Visualising initial condition
fig = plt.figure(figsize=(14,8),dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x,y)
surf = ax.plot_surface(X,Y,u[:])
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.show()