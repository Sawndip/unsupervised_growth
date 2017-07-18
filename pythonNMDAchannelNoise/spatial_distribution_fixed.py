# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 09:39:36 2016

@author: jingroup
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May 11 12:41:50 2016

@author: eugene

Script plots spatial distributions of fixed synapses
"""
import analyze_spatial_dist as spatial
import reading
import math
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import os

dirname = "/home/eugene/Output/networks/sphere_170717_hodgkin/"
dim = 3 # network dimensionality

latest_checkpoint = spatial.find_latest_checkpoint(dirname) 

print "latest checkpoint = ",latest_checkpoint

#dirname = "/home/eugene/lionX/clustered/network/"
RA_xy = os.path.join(dirname, "RA_xy_" + str(latest_checkpoint) + "_.bin")
I_xy = os.path.join(dirname, "I_xy.bin")
RA2I = os.path.join(dirname, "RA_I_connections_" + str(latest_checkpoint) + "_.bin")
I2RA = os.path.join(dirname, "I_RA_connections_" + str(latest_checkpoint) + "_.bin")


SQUARE_SIDE = 100
SIDE = SQUARE_SIDE * math.sqrt(2)
R = 1.0 # radius of sphere

def distance_on_sphere(v1, v2, R):
    """
    Computes distance between two points on sphere of radius R
    """
    return np.arccos(np.dot(v1, v2)) / R

coord_I= reading.read_coordinates(dim, I_xy)
coord_RA = reading.read_coordinates(dim, RA_xy)

(N_RA, RA_targets, RA_targets_G) = reading.read_connections(RA2I)
(N_I, I_targets, I_targets_G) = reading.read_connections(I2RA)

num_targetsRAI = [len(t) for t in RA_targets]
num_targetsIRA = [len(t) for t in I_targets]
targetsRAI = [t for sublist in RA_targets for t in sublist]
targetsIRA = [t for sublist in I_targets for t in sublist]

# Here we calculate separation between I neurons
distances_between_all_I = np.empty(shape=(N_I, N_I), dtype=np.float32)
distances_between_all_I.fill(1e6)

for i in xrange(N_I):
    for j in xrange(i+1, N_I):
        if dim == 2:
            distances_between_all_I[i][j] = np.linalg.norm(coord_I[i]-coord_I[j]) / SIDE
        if dim == 3:
            distances_between_all_I[i][j] = distance_on_sphere(coord_I[i], coord_I[j], R)

distances_between_all_I[np.tril_indices(N_I, -1)] = distances_between_all_I.T[np.tril_indices(N_I, -1)]

min_distances = np.empty(N_I, np.float32)
for i in range(N_I):
    min_distances[i] = np.min(distances_between_all_I[i])

distance_between_interneurons = np.mean(min_distances)

print "Average distance between interneurons: ",distance_between_interneurons



#mean_distance_RAandI = np.mean(distances_between_all_RA_and_I)

# Here we calculate mean distance between RA and their I targets

distances_between_RA_and_their_I_targets = []

for i in xrange(N_RA):
    if len(RA_targets[i]) > 1:
        for target_ind in RA_targets[i]:
            if dim == 2:     
                distances_between_RA_and_their_I_targets.append(np.linalg.norm(coord_RA[i]-coord_I[target_ind]) / (SIDE*distance_between_interneurons))
            if dim == 3:
                distances_between_RA_and_their_I_targets.append(distance_on_sphere(coord_RA[i], coord_I[target_ind], R) / distance_between_interneurons)

# Here we calculate mean distance between I and their RA targets

distances_between_I_and_their_RA_targets = []

for i in xrange(N_I):
    if len(I_targets[i]) > 1:
        for target_ind in I_targets[i]:
            if dim == 2:     
                distances_between_I_and_their_RA_targets.append(np.linalg.norm(coord_I[i]-coord_RA[target_ind])/(SIDE*distance_between_interneurons))
            if dim == 3:
                distances_between_I_and_their_RA_targets.append(distance_on_sphere(coord_I[i], coord_RA[target_ind], R) / distance_between_interneurons)
       

# let's count number of input I connections

numInputI = Counter(targetsIRA).values()
numWithIinputs = len(numInputI)
numInputI.extend([0]*(N_RA - numWithIinputs))

# let's count number of input RA connections
numInputRA = Counter(targetsRAI).values()
numWithRAinputs = len(numInputRA)
numInputRA.extend([0]*(N_I - numWithRAinputs))


print "Mean distance for RA to I connections normalized by distance between interneurons: ", np.mean(distances_between_RA_and_their_I_targets)
print "Standard deviation for distances from RA to I: ", np.std(distances_between_RA_and_their_I_targets)
print "Mean distance for I to RA connections normalized by distance between interneurons: ", np.mean(distances_between_I_and_their_RA_targets)
print "Standard deviation for distances from I to RA: ", np.std(distances_between_I_and_their_RA_targets)

numBins = 30

f1 = plt.figure()
ax = f1.add_subplot(111)
ax.hist(distances_between_RA_and_their_I_targets, numBins)
ax.set_xlabel("Normalized distance")
ax.set_ylabel("Number of connections")
ax.set_title("Spatial distribution of connections from RA onto I neurons")
ax.set_xlim([0,5])

numBins = 25

f2 = plt.figure()
ax = f2.add_subplot(111)
ax.hist(distances_between_I_and_their_RA_targets, numBins)
ax.set_xlabel("Normalized distance")
ax.set_ylabel("Number of connections")
ax.set_title("Spatial distribution of connections from I onto RA neurons")
ax.set_xlim([0,5])

numBins = 25


f3 = plt.figure()
ax = f3.add_subplot(111)
ax.hist(num_targetsRAI, numBins)
ax.set_xlabel("# of targets")
ax.set_title("Distribution of number of targets from RA to I neurons")

f4 = plt.figure()
ax = f4.add_subplot(111)
ax.hist(num_targetsIRA, numBins)
ax.set_xlabel("# of targets")
ax.set_title("Distribution of number of targets from I to RA neurons")

f5 = plt.figure()
ax = f5.add_subplot(111)
ax.hist(numInputI, numBins)
ax.set_xlabel("# of inputs from interneurons")
ax.set_title("Distribution of number of inputs from I to RA neurons")

f6 = plt.figure()
ax = f6.add_subplot(111)
ax.hist(numInputRA, numBins)
ax.set_xlabel("# of inputs from RA neurons")
ax.set_title("Distribution of number of inputs from RA to I neurons")

plt.show()
