# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 12:31:07 2018

@author: jingroup

Script plots location of mature neurons using Mollweide projection of sphere
"""
import matplotlib.pyplot as plt
import math
import numpy as np
import utils
import reading
import os

R = 1.0 # radius of sphere

num_iter = 10000
tolerance = 1e-7

dirname = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/events1/"
trial_number = 1700

fileMature = os.path.join(dirname, "mature_" + str(trial_number) + ".bin")
fileCoordRA = os.path.join(dirname, "RA_xy_" + str(trial_number) + ".bin")
fileTraining = os.path.join(dirname, "training_neurons.bin")

coord = reading.read_coordinates(fileCoordRA)
(_, _, mature_indicators) = reading.read_mature_indicators(fileMature)
training_neurons = reading.read_training_neurons(fileTraining)

mature_neurons = np.where(mature_indicators == 1)[0]

coord_mature = coord[mature_neurons]

### calculate latitude and longitude of mature HVC-RA neurons ####
longitude, latitude = utils.calculate_longAndLat(coord_mature)

print min(longitude) / np.pi
print max(longitude) / np.pi
print min(latitude) / (np.pi/2.)
print max(latitude) / (np.pi/2.)

Mollweide_RA = np.empty(shape=(len(mature_neurons),2), dtype=np.float32) # array with Mollweide coordinates of HVC(I) neurons
Mollweide_RA.fill(np.nan)


for i in range(len(mature_neurons)):
    Mollweide_RA[i,0], Mollweide_RA[i,1] = utils.Mollweide_projection(longitude[i], latitude[i], tolerance, num_iter)


colors = []
s = []

for i in mature_neurons:
    if i in training_neurons:
        colors.append('m')
        s.append(30)
    else:
        colors.append('b')
        s.append(10)

fig = plt.figure()
#ax = Axes3D(fig)
ax = fig.add_subplot(111)
ax.scatter(Mollweide_RA[:,0], Mollweide_RA[:,1], c=colors, s=s)

# plot meridians
num_meridians = 10
longitude = np.linspace(-np.pi, np.pi, num_meridians)

for i in range(num_meridians):    
    meridian = utils.get_meridian(longitude[i], 1000, tolerance, num_iter)
    ax.plot(meridian[0], meridian[1], c='k')

# plot parallels
num_parallels = 10
latitude = np.linspace(-np.pi/2.0, np.pi/2.0, num_parallels)

for i in range(num_parallels):    
    parallel = utils.get_parallel(latitude[i], 1000, tolerance, num_iter)
    ax.plot(parallel[0], parallel[1], c='k')

#ax.scatter(coord_I[:,0], coord_I[:,1], coord_I[:,2], c='r')
#ax.scatter(coord_RA[:,0], coord_RA[:,1], coord_RA[:,2], c='b')

plt.show()

