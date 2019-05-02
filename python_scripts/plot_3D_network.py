# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 16:26:01 2018

@author: jingroup
"""
import reading
import os
import matplotlib.pyplot as plt
import numpy as np

def set_axes_radius(ax, origin, radius):
    ax.set_xlim3d([origin[0] - radius, origin[0] + radius])
    ax.set_ylim3d([origin[1] - radius, origin[1] + radius])
    ax.set_zlim3d([origin[2] - radius, origin[2] + radius])

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])

    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    set_axes_radius(ax, origin, radius)

from mpl_toolkits.mplot3d import Axes3D

dirname = "/mnt/hodgkin/eugene/Output/networks/chainGrowth/network100RA27I"
fileTraining = "/mnt/hodgkin/eugene/Output/networks/chainGrowth/network100RA27I/training_neurons_random.bin"

#dirname = "/mnt/hodgkin/eugene/results/immature/clusters/matTrans62/"


training_neurons = set(reading.read_training_neurons(fileTraining))

print training_neurons

#coord_HVCRA = reading.read_coordinates(os.path.join(dirname, "RA_xy_4200.bin"))
coord_HVCRA = reading.read_coordinates(os.path.join(dirname, "RA_xy.bin"))

coord_HVCI = reading.read_coordinates(os.path.join(dirname, "I_xy.bin"))    

color_HVCRA = ['r' for i in range(coord_HVCRA.shape[0])]

for i in training_neurons:
    color_HVCRA[i] = 'g'
    
f = plt.figure()
ax = Axes3D(f)
ax.scatter(coord_HVCRA[:,0], coord_HVCRA[:,1], coord_HVCRA[:,2], c=color_HVCRA)
ax.scatter(coord_HVCI[:,0], coord_HVCI[:,1], coord_HVCI[:,2], c='b')

ax.set_aspect('equal')         # important!

ax.set_xlim3d([-1.0, 1.0])
ax.set_ylim3d([-1.0, 1.0])
ax.set_zlim3d([-1.0, 1.0])

#set_axes_equal(ax)

plt.show()
