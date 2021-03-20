"""
Reads saved data and uses an interpolation method to create an evenly spaced grid of data points.
This is done in a way such that all the old plotting code should work with the rescaled arrays.
Make a copy of this file and edit that rather than directly editing this file.

Array definitions

Original data
-------------
u - x-velocity
v - y-velocity
p - pressure
X - x-coordinates of the data points
Y - y-coordinates of the data points

Rescaled data
-------------
u_rescale - x-velocity
v_rescale - y-velocity
p_rescale - pressure
X_rescale - x-coordinates of the data points
Y_rescale - y-coordinates of the data points
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

#I/O configuration for unblocked channel
fname = 'data_unblocked.npz'
fpath = 'Output\\' + fname

#Loading the data from the file into numpy arrays
data = np.load(fpath)
u  = data['x_vel']
v = data['y_vel']
p = data['pressure']
X = data['x_coord']
Y = data['y_coord']

#streamplot with coloured velocity for unblocked data
fig, ax = plt.subplots()

stream = ax.streamplot(X,Y,u,v, density = 1, color=u, cmap=plt.cm.Reds, )
ax.set_title('Fluid Velocities in Unblocked Channel')
plt.xlabel("Length along vessel")
plt.ylabel("Length across vessel")
colours = plt.colorbar(stream.lines)
colours.set_label("Velocity of fluid")
plt.show()

#streamplot with coloured pressures for unblocked data
fig, ax = plt.subplots()

stream = ax.streamplot(X,Y,u,v, density = 1, color=p, cmap=plt.cm.Reds, )
ax.set_title('Fluid Velocities in Unblocked Channel')
plt.xlabel("Length along vessel")
plt.ylabel("Length across vessel")
colours = plt.colorbar(stream.lines)
colours.set_label("Pressure of fluid particles")
plt.show()


#I/O configuration for blocked channels
fname = 'data_blocked.npz'
fpath = 'Output\\' + fname

#Loading the data from the file into numpy arrays
data = np.load(fpath)
u  = data['x_vel']
v = data['y_vel']
p = data['pressure']
X = data['x_coord']
Y = data['y_coord']

#Grid parameters
x_min = 0
x_max = 100
nx = 500
y_min = 0
y_max = 12
ny = 150

#Initialising points for evenly spaced data points

rescaled_x_axis = np.linspace(x_min,x_max,nx)
rescaled_y_axis = np.linspace(y_min,y_max,ny)
X_rescale,Y_rescale = np.meshgrid(rescaled_x_axis,rescaled_y_axis) #Creating rescaled coordinate grid

original_coords = (X.flatten(),Y.flatten())
rescaled_coords = (X_rescale,Y_rescale)


#Using cubic interpolation to obtain data points in new grid
u_rescale = griddata(original_coords, u.flatten(), rescaled_coords, method='cubic')
v_rescale = griddata(original_coords, v.flatten(), rescaled_coords, method='cubic')
p_rescale = griddata(original_coords, p.flatten(), rescaled_coords, method='cubic')

#Filtering computational artefacts from interpolation
#Filter thresholds
u_threshold = 0.005

#Setting velocities under threshold to 0 (As extremely small values are likely a consequence of floating point precision)
#Made assumption that computational artefacts are characterised by x-velocity
u_rescale[abs(u_rescale)<u_threshold] = 0
v_rescale[abs(u_rescale)<u_threshold] = 0
#Setting corresponding pressures to 0
p_rescale[abs(u_rescale)<u_threshold] = 0

#Plotting code goes below this line.
#Calculating magnitude of velocity (speed)
v_mag = np.sqrt(u_rescale**2 + v_rescale**2)

#Streamplot with coloured velocities
fig, ax = plt.subplots()
stream = ax.streamplot(X_rescale, Y_rescale, u_rescale, v_rescale, density = 1, color=p_rescale, cmap=plt.cm.Reds, )
ax.set_title('Fluid Velocities in blocked Channel')
plt.xlabel("Length along vessel")
plt.ylabel("Length across vessel")
colours = plt.colorbar(stream.lines)
colours.set_label("Pressure of fluid particles")
plt.show()

#Streamplot with coloured pressures
fig, ax = plt.subplots()
stream = ax.streamplot(X_rescale, Y_rescale, u_rescale, v_rescale, density = 1, color=v_mag, cmap=plt.cm.Reds, )
ax.set_title('Fluid Velocities in blocked Channel')
plt.xlabel("Length along vessel")
plt.ylabel("Length across vessel")
colours = plt.colorbar(stream.lines)
colours.set_label("Velocities of fluid particles")
plt.show()

