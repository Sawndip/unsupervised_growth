# -*- coding: utf-8 -*-
"""
Created on Wed May 11 12:41:50 2016

@author: eugene
"""
import os
import math

dirname = "/home/eugene/Output/networks/gabaMaturation020517/"
#dirname = "/home/eugene/lionX/clustered/network/"
RA_xy = os.path.join(dirname, "RA_xy.bin")
I_xy = os.path.join(dirname, "I_xy.bin")
RA2I = os.path.join(dirname, "RA_I_connections.bin")
I2RA = os.path.join(dirname, "I_RA_connections.bin")
#RARA = os.path.join(dirname, "RA_RA_super_connections.bin")
RARA = "/home/eugene/Output/networks/gabaMaturation020517/RA_RA_super_connections.bin"
#RARA = "/home/eugene/lionX/clustered/r0.05/gabaMaturation180417_0/RA_RA_super_connections.bin"

import reading
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import indirect_connections as connect

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
(N_RA, RA_super_targets, RA_super_targets_G) = reading.read_connections(RARA)



print "RA_super_targets: ", RA_super_targets


num_targetsRAI = [len(t) for t in RA_targets]
num_targetsIRA = [len(t) for t in I_targets]
targetsRAI = [t for sublist in RA_targets for t in sublist]
targetsIRA = [t for sublist in I_targets for t in sublist]

# code for super.net file
#==============================================================================
# startingString = N_RA + 2
# 
# with open(filename, "r") as f:
#     lines = f.readlines()
# f.close()
# 
# data = lines[startingString:]
# #print data
# 
# source = []
# target = []
# 
# for line in data:
#     line = line.strip()
#     columns = line.split()
#     source.append(int(columns[0]) - 1)
#     target.append(int(columns[1]) - 1)
# 
#==============================================================================
#print "source: ", source
#print "target: ", target

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


# Here we calculate mean distance between RA and their sypersynapses (version for super.net file)

#==============================================================================
# dist_RA2RA = []
# 
# for i in xrange(len(source)):
#     target_ind = target[i]
#     source_ind = source[i]
#     dist_RA2RA.append(d(x_RA[source_ind], y_RA[source_ind], x_RA[target_ind], y_RA[target_ind])/SIDE)
# #print "RA_targets: ", RA_targets
#==============================================================================
#print "I_targets_G: ", I_targets_G

# Here we calculate mean distance between RA and their sypersynapses (version for super.net file)

dist_RA2RA = []

for i in xrange(N_RA):
    if len(RA_super_targets[i]) > 0:
        for target_ind in RA_super_targets[i]:
            dist_RA2RA.append(d(x_RA[i], y_RA[i], x_RA[target_ind], y_RA[target_ind])/SIDE)

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
            

# Here we plot spatial distribution of output connections for each HVC(RA) neurons
training_neurons = [0, 1, 2, 3]

chain_layers = connect.find_chain_layers(RA_super_targets, training_neurons)

chain_neurons = [i for layer in chain_layers for i in layer]

print "Chain neurons: ",chain_neurons

source_neurons = [] # list of source neuron ids
synapse_distances = [] # list of output synapse distances for each source neuron

mean_synapse_distances = [] # average output synapse distance for each neuron
std_synapse_distances = [] # standard deviation of output synapse distance for each neuron

target_RA_x = [] # x-coordinates of target HVC(RA) neurons
target_RA_y = [] # y-coordinates of target HVC(RA) neurons

target_I_x = [] # x-coordinates of target HVC(I) neurons
target_I_y = [] # y-coordinates of target HVC(I) neurons


for neuron in chain_neurons:
    if len(RA_super_targets[neuron]) > 0:
        source_neuron = []
        synapse_distance = []
        coord_RA_x = []
        coord_RA_y = []
        
        # get coordinates of HVC(RA) targets
        for target_ind in RA_super_targets[neuron]:
            source_neuron.append(neuron)
            synapse_distance.append(d(x_RA[neuron], y_RA[neuron], x_RA[target_ind], y_RA[target_ind])/SIDE)
            coord_RA_x.append(x_RA[target_ind])
            coord_RA_y.append(y_RA[target_ind])
            

        source_neurons.append(source_neuron)
        synapse_distances.append(synapse_distance)
        
        # calculate mean and variance of synapse distances
        mean_synapse_distances.append(np.mean(synapse_distance))
        std_synapse_distances.append(np.std(synapse_distance))
        
        target_RA_x.append(coord_RA_x)
        target_RA_y.append(coord_RA_y)
        
        coord_I_x = []
        coord_I_y = []
        
         # get coordinates of HVC(I) targets
        for target_ind in RA_targets[neuron]:
            coord_I_x.append(x_I[target_ind])
            coord_I_y.append(y_I[target_ind])
            
        target_I_x.append(coord_I_x)
        target_I_y.append(coord_I_y)
        
            
        

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
print "Mean distance for RA to RA connections normalized by side of square area: ", np.mean(dist_RA2RA)
if len(dist_RA2RA) > 1:
    print "Standard deviation for distances from RA to RA: ", std(dist_RA2RA)

numBins = 30

f1 = plt.figure()
plt.suptitle("Spatial distributions")
ax = f1.add_subplot(311)
ax.hist(dist_RA2I, numBins)
#ax.set_xlabel("Normalized distance")
ax.set_ylabel("Number of connections")
ax.set_title("HVC(RA) to HVC(I)")
ax.set_xlim([0,1])

numBins = 25

ax = f1.add_subplot(312)
ax.hist(dist_I2RA, numBins)
#ax.set_xlabel("Normalized distance")
ax.set_ylabel("Number of connections")
ax.set_title("HVC(I) to HVC(RA)")
ax.set_xlim([0,1])

numBins = 35

ax = f1.add_subplot(313)
ax.hist(dist_RA2RA, numBins)
ax.set_xlabel("Normalized distance")
ax.set_ylabel("Number of connections")
ax.set_title("HVC(RA) to HVC(RA)")
ax.set_xlim([0,1])

#==============================================================================
# f4 = plt.figure()
# ax = f4.add_subplot(111)
# ax.hist(num_targetsRAI, numBins)
# ax.set_xlabel("# of targets")
# ax.set_title("Distribution of number of targets from RA to I neurons")
# 
# f5 = plt.figure()
# ax = f5.add_subplot(111)
# ax.hist(num_targetsIRA, numBins)
# ax.set_xlabel("# of targets")
# ax.set_title("Distribution of number of targets from I to RA neurons")
# 
# f5 = plt.figure()
# ax = f5.add_subplot(111)
# ax.hist(num_targetsIRA, numBins)
# ax.set_xlabel("# of targets")
# ax.set_title("Distribution of number of targets from I to RA neurons")
# 
#==============================================================================


#==============================================================================
# f6 = plt.figure()
# ax = f6.add_subplot(111)
# 
# for i in range(len(source_neurons)):
#     ax.plot(synapse_distances[i], source_neurons[i], marker='x', markersize=12, linewidth=2)
#         
# ax.set_xlabel("synapse distance")
# ax.set_ylabel("neuron id")
# ax.set_title("Distribution of output synaptic distances for each neuron")
# 
# numBins = 35
# 
#==============================================================================
#==============================================================================
# f7 = plt.figure()
# ax = f7.add_subplot(111)
# ax.hist(std_synapse_distances, numBins)
# ax.set_xlabel("$\sigma_{distance}$")
# ax.set_title("Distribution of standard deviations of synaptic distances")
# 
#==============================================================================
#==============================================================================
# f8 = plt.figure()
# ax = f8.add_subplot(111)
# ax.scatter(mean_synapse_distances, std_synapse_distances)
# ax.set_xlabel("mean synapse distance")
# ax.set_ylabel("$\sigma_{distance}$")
# 
# 
# ax.set_title("Standard deviation of synapse distance vs mean synapse distance")
# 
#==============================================================================

#f6 = plt.figure()
#ax = f6.add_subplot(111)
print "number of neurons with supersynaptic targets: ",len(source_neurons)
print source_neurons

num_neurons_in_plot = 16

num_cols = 4
num_rows = 4

num_plots = len(source_neurons) / num_neurons_in_plot # number of fully covered plots

num_neurons_last_plot = len(source_neurons) - num_plots * num_neurons_in_plot
num_rows_last_plot = num_neurons_last_plot / num_rows + 1
num_cols_last_plot = num_neurons_last_plot % 4


print "num full plots: ",num_plots
print "num rows in last plot: ",num_rows_last_plot
print "num cols in last plot: ",num_cols_last_plot



# plot full plots first
for i in range(num_plots):    
    f, axarr = plt.subplots(num_rows, num_cols)
    start_id = i*num_neurons_in_plot
    
    
    for j in range(num_rows):
        for k in range(num_cols):
            source_id = source_neurons[start_id+j*num_cols+k][0]
            
            for m in range(len(target_RA_x[start_id+j*num_cols+k])):
                
                axarr[j,k].plot([x_RA[source_id], target_RA_x[start_id+j*num_cols+k][m]], [y_RA[source_id], target_RA_y[start_id+j*num_cols+k][m]], color='red')            
            
            axarr[j,k].scatter(target_I_x[start_id+j*num_cols+k], target_I_y[start_id+j*num_cols+k], color='red', marker='x', s=50)            
            
            #axarr[j,k].scatter(target_x[start_id+j*num_cols+k], target_y[start_id+j*num_cols+k], marker='x', s=50)
            axarr[j,k].set_xlim([0, SQUARE_SIDE])
            axarr[j,k].set_ylim([0, SQUARE_SIDE])

# plotting the last figure
f, axarr = plt.subplots(num_rows, num_cols)
start_id = num_plots*num_neurons_in_plot

for j in range(num_rows_last_plot):
    for k in range(num_cols):
        
        if j == num_rows_last_plot - 1 and k >= num_cols_last_plot:
            break
        
        source_id = source_neurons[start_id+j*num_cols+k][0]
        
        for m in range(len(target_RA_x[start_id+j*num_cols+k])):
            
            axarr[j,k].plot([x_RA[source_id], target_RA_x[start_id+j*num_cols+k][m]], [y_RA[source_id], target_RA_y[start_id+j*num_cols+k][m]], color='red')            
        
        axarr[j,k].scatter(target_I_x[start_id+j*num_cols+k], target_I_y[start_id+j*num_cols+k], color='red', marker='x', s=50)            
        
        #axarr[j,k].scatter(target_x[start_id+j*num_cols+k], target_y[start_id+j*num_cols+k], marker='x', s=50)
        axarr[j,k].set_xlim([0, SQUARE_SIDE])
        axarr[j,k].set_ylim([0, SQUARE_SIDE])
             
#ax.set_xlabel("synapse distance")
#ax.set_ylabel("neuron id")
#ax.set_title("Distribution of output synaptic distances for each neuron")


plt.show()