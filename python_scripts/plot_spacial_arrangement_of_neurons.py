# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 20:25:10 2018

@author: jingroup

Script analyzes spatial arrangment of neurons on sphere
"""
import matplotlib.pyplot as plt
import reading
import os
import space
import numpy as np

directory = "/home/eugene/Output/networks/chainGrowth/network2000/"

file_RA_coord = os.path.join(directory, "RA_xy.bin")
file_I_coord = os.path.join(directory, "I_xy.bin")

coord_RA = reading.read_coordinates(file_RA_coord)
coord_I = reading.read_coordinates(file_I_coord)

# find distribution of smallest distances between interneurons
N_I = coord_I.shape[0]

smallest_distances = np.empty(N_I, np.float32)

for i in range(N_I):
    min_distance = 1e6

    for j in range(N_I):
        if i != j:
            distance = space.distance_on_sphere(coord_I[i], coord_I[j], 1.0)
            if distance < min_distance:
                min_distance = distance

    smallest_distances[i] = min_distance

##########################
### spatial distribution
##########################
plt.figure()
plt.hist(smallest_distances)
plt.xlabel('Distance to nearest HVC-I neighbour')
plt.ylabel('# of occurences')
plt.xlim([0.0,0.30])


########################
###   HVC-I neurons
########################

f = plt.figure()

ax1 = f.add_subplot(221)
ax1.scatter(coord_I[:,0], coord_I[:,1])
ax1.set_xlim([-1.0, 1.0])
ax1.set_ylim([-1.0, 1.0])
ax1.set_title('HVC-I xy-plane')
ax1.axis('equal')

ax2 = f.add_subplot(222)
ax2.scatter(coord_I[:,1], coord_I[:,2])
ax2.set_xlim([-1.0, 1.0])
ax2.set_ylim([-1.0, 1.0])
ax2.set_title('HVC-I yz-plane')
ax2.axis('equal')

ax3 = f.add_subplot(223)
ax3.scatter(coord_I[:,0], coord_I[:,2])
ax3.set_xlim([-1.0, 1.0])
ax3.set_ylim([-1.0, 1.0])
ax3.set_title('HVC-I xz-plane')
ax3.axis('equal')

########################
###   HVC-RA neurons
########################

f = plt.figure()

ax1 = f.add_subplot(221)
ax1.scatter(coord_RA[:,0], coord_RA[:,1])
ax1.set_xlim([-1.0, 1.0])
ax1.set_ylim([-1.0, 1.0])
ax1.set_title('HVC-RA xy-plane')
ax1.axis('equal')

ax2 = f.add_subplot(222)
ax2.scatter(coord_RA[:,1], coord_RA[:,2])
ax2.set_xlim([-1.0, 1.0])
ax2.set_ylim([-1.0, 1.0])
ax2.set_title('HVC-RA yz-plane')
ax2.axis('equal')

ax3 = f.add_subplot(223)
ax3.scatter(coord_RA[:,0], coord_RA[:,2])
ax3.set_xlim([-1.0, 1.0])
ax3.set_ylim([-1.0, 1.0])
ax3.set_title('HVC-RA xz-plane')
ax3.axis('equal')


plt.show()