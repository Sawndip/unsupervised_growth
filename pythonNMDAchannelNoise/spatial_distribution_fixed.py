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

filename = "/home/eugene/Output/super.net"

RA_xy = "/home/eugene/Output/RA_xy.bin"
I_xy = "/home/eugene/Output/I_xy.bin"
RA2I = "/home/eugene/Output/RA_I_connections.bin"
I2RA = "/home/eugene/Output/I_RA_connections.bin"


import reading
import math
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

SQUARE_SIDE = 100
SIDE = SQUARE_SIDE * math.sqrt(2)



def d(x1, y1, x2, y2):
    """ Computes distance between two points in space"""
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

def std(x):
    """Computes standard deviation"""
    mean = np.mean(x)
    
    sum = 0    
    for e in x:
        sum += (e - mean)**2
        
    return math.sqrt(sum/len(x))
        



(x_I, y_I) = reading.read_coordinates(I_xy)
(x_RA, y_RA) = reading.read_coordinates(RA_xy)



print "(x_I, y_I):", zip(x_I, y_I)
print "(x_RA, y_RA):", zip(x_RA, y_RA)


(N_RA, RA_targets, RA_targets_G) = reading.read_connections(RA2I)
(N_I, I_targets, I_targets_G) = reading.read_connections(I2RA)

num_targetsRAI = [len(t) for t in RA_targets]
num_targetsIRA = [len(t) for t in I_targets]
targetsRAI = [t for sublist in RA_targets for t in sublist]
targetsIRA = [t for sublist in I_targets for t in sublist]

# Here we calculate mean distance between all RA and I neurons

all_dist_RAandI = []

for i in xrange(N_RA):
    for j in xrange(N_I):
        dist = d(x_RA[i], y_RA[i], x_I[j], y_I[j])
        all_dist_RAandI.append(dist)

mean_distance_RAandI = np.mean(all_dist_RAandI)



# Here we calculate mean distance between all RA and RA neurons
        
all_dist_RAandRA = []


for i in xrange(N_RA):
    for j in xrange(i, N_RA):
        if j!=i:
            dist = d(x_RA[i], y_RA[i], y_RA[j], y_RA[j])
        all_dist_RAandRA.append(dist)

mean_distance_RAandRA = np.mean(all_dist_RAandRA)



# Here we calculate mean distance between RA and their I targets

dist_RA2I = []

for i in xrange(N_RA):
    if len(RA_targets[i]) > 1:
        for target_ind in RA_targets[i]:
            dist_RA2I.append(d(x_RA[i], y_RA[i], x_I[target_ind], y_I[target_ind])/SIDE)


# Here we calculate mean distance between I and their RA targets

dist_I2RA = []

for i in xrange(N_I):
    if len(I_targets[i]) > 1:
        for target_ind in I_targets[i]:
            dist_I2RA.append(d(x_I[i], y_I[i], x_RA[target_ind], y_RA[target_ind])/SIDE)
            

# let's count number of input I connections

numInputI = Counter(targetsIRA).values()
numWithIinputs = len(numInputI)
numInputI.extend([0]*(N_RA - numWithIinputs))

# let's count number of input RA connections
numInputRA = Counter(targetsRAI).values()
numWithRAinputs = len(numInputRA)
numInputRA.extend([0]*(N_I - numWithRAinputs))

#print dist_RA2I
#print "test:", d(17,31,33,45)
print "Mean distance between all RA and I neurons normalized by side of square area: ", mean_distance_RAandI/SIDE
print "Standard deviation for all distances between RA and I: ", std(all_dist_RAandI)/SIDE
print "Mean distance between all RA and RA neurons normalized by side of square area: ", mean_distance_RAandRA/SIDE
print "Standard deviation for all distances between RA and RA: ", std(all_dist_RAandRA)/SIDE

print "Mean distance for RA to I connections normalized by side of square area: ", np.mean(dist_RA2I)
print "Standard deviation for distances from RA to I: ", std(dist_RA2I)
print "Mean distance for I to RA connections normalized by side of square area: ", np.mean(dist_I2RA)
print "Standard deviation for distances from I to RA: ", std(dist_I2RA)

numBins = 30

f1 = plt.figure()
ax = f1.add_subplot(111)
ax.hist(dist_RA2I, numBins)
ax.set_xlabel("Normalized distance")
ax.set_ylabel("Number of connections")
ax.set_title("Spatial distribution of connections from RA onto I neurons")
ax.set_xlim([0,1])

numBins = 25

f2 = plt.figure()
ax = f2.add_subplot(111)
ax.hist(dist_I2RA, numBins)
ax.set_xlabel("Normalized distance")
ax.set_ylabel("Number of connections")
ax.set_title("Spatial distribution of connections from I onto RA neurons")
ax.set_xlim([0,1])

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