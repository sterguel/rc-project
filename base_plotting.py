"""
File containing the basic code to read/write out of the npz files that the simulated data are saved to.
Make a copy of this file and edit that rather than directly editing this file.

Array definitions

u - x-velocity
v - y-velocity
p - pressure
x_coord - x-coordinates of the data points
y_coord - y-coordinates of the data points
"""

import numpy as np
import matplotlib.pyplot as plt

#I/O configuration
fname = 'data_1.npz'
fpath = 'Output\\' + fname

#Loading the data from the file into numpy arrays
data = np.load(fpath)
u  = data['x_vel']
v = data['y_vel']
p = data['pressure']
X = data['x_coord']
Y = data['y_coord']

#Plotting code goes below this line.
