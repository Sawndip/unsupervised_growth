# -*- coding: utf-8 -*-
"""
Created on Wed May 11 12:41:50 2016

@author: eugene
"""
import analyze_spatial_dist as spatial
import reading
import math
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import os
import space

#dirname = "/home/eugene/Output/networks/sphere_180717_huxley/"
dirname = "/home/eugene/results/noDelays/replacement/sphere/180717_lionx_2/"
arrangement = "sphere" # spatial arrangement of neurons

distance = space.distance_function[arrangement] # define distance function for current arrangement
dim = space.num_coordinates[arrangement] # number of coordinates to define neuronal location in space

latest_checkpoint = spatial.find_latest_checkpoint(dirname) 

print "latest checkpoint = ",latest_checkpoint

#dirname = "/home/eugene/lionX/clustered/network/"
RA_xy = os.path.join(dirname, "RA_xy_" + str(latest_checkpoint) + "_.bin")
I_xy = os.path.join(dirname, "I_xy.bin")
RA2I = os.path.join(dirname, "RA_I_connections_" + str(latest_checkpoint) + "_.bin")
I2RA = os.path.join(dirname, "I_RA_connections_" + str(latest_checkpoint) + "_.bin")
RARA = os.path.join(dirname, "RA_RA_super_connections_" + str(latest_checkpoint) + "_.bin")


coord_I= reading.read_coordinates(dim, I_xy)
coord_RA = reading.read_coordinates(dim, RA_xy)

(N_RA, RA_targets, RA_targets_G) = reading.read_connections(RA2I)
(N_I, I_targets, I_targets_G) = reading.read_connections(I2RA)
(_, RA_super_targets, RA_super_targets_G) = reading.read_connections(RARA)

num_targetsRAI = [len(t) for t in RA_targets]
num_targetsIRA = [len(t) for t in I_targets]
targetsRAI = [t for sublist in RA_targets for t in sublist]
targetsIRA = [t for sublist in I_targets for t in sublist]

# Here we calculate separation between I neurons
distances_between_all_I = np.empty(shape=(N_I, N_I), dtype=np.float32)
distances_between_all_I.fill(1e6)

for i in xrange(N_I):
    for j in xrange(i+1, N_I):
        distances_between_all_I[i][j] = distance(coord_I[i], coord_I[j])

distances_between_all_I[np.tril_indices(N_I, -1)] = distances_between_all_I.T[np.tril_indices(N_I, -1)]

min_distances = np.empty(N_I, np.float32)
for i in range(N_I):
    min_distances[i] = np.min(distances_between_all_I[i])

distance_between_interneurons = np.mean(min_distances)

print "Average distance between interneurons: ",distance_between_interneurons

# Here we calculate distances between RA and their I targets

distances_between_RA_and_their_I_targets = []

for i in xrange(N_RA):
    if len(RA_targets[i]) > 1:
        for target_ind in RA_targets[i]:
            distances_between_RA_and_their_I_targets.append(distance(coord_RA[i], coord_I[target_ind]) / distance_between_interneurons)
         
# Here we calculate distances between I and their RA targets

distances_between_I_and_their_RA_targets = []

for i in xrange(N_I):
    if len(I_targets[i]) > 1:
        for target_ind in I_targets[i]:
            distances_between_I_and_their_RA_targets.append(distance(coord_I[i], coord_RA[target_ind]) / distance_between_interneurons)
            
# Here we calculate distances between HVC(RA) neurons and their super targets
distances_between_RA_and_their_RA_targets = []

for i in xrange(N_RA):
    if len(RA_super_targets[i]) > 1:
        for target_ind in RA_super_targets[i]:
            distances_between_RA_and_their_RA_targets.append(distance(coord_RA[i],coord_RA[target_ind]) / distance_between_interneurons)
            
    
print "Mean distance for RA to I connections normalized by distance between interneurons: ", np.mean(distances_between_RA_and_their_I_targets)
print "Standard deviation for distances from RA to I: ", np.std(distances_between_RA_and_their_I_targets)
print "Mean distance for I to RA connections normalized by distance between interneurons: ", np.mean(distances_between_I_and_their_RA_targets)
print "Standard deviation for distances from I to RA: ", np.std(distances_between_I_and_their_RA_targets)
print "Mean distance for RA to RA connections normalized by distance between interneurons: ", np.mean(distances_between_RA_and_their_RA_targets)
print "Standard deviation for distances from RA to RA: ", np.std(distances_between_RA_and_their_RA_targets)



            

#==============================================================================
# # Here we plot spatial distribution of output connections for each HVC(RA) neurons
# training_neurons = [0, 1, 2, 3]
# 
# chain_layers = connect.find_chain_layers(RA_super_targets, training_neurons)
# 
# chain_neurons = [i for layer in chain_layers for i in layer]
# 
# print "Chain neurons: ",chain_neurons
# 
# source_neurons = [] # list of source neuron ids
# synapse_distances = [] # list of output synapse distances for each source neuron
# 
# mean_synapse_distances = [] # average output synapse distance for each neuron
# std_synapse_distances = [] # standard deviation of output synapse distance for each neuron
# 
# target_RA_x = [] # x-coordinates of target HVC(RA) neurons
# target_RA_y = [] # y-coordinates of target HVC(RA) neurons
# 
# target_I_x = [] # x-coordinates of target HVC(I) neurons
# target_I_y = [] # y-coordinates of target HVC(I) neurons
# 
# 
# for neuron in chain_neurons:
#     if len(RA_super_targets[neuron]) > 0:
#         source_neuron = []
#         synapse_distance = []
#         coord_RA_x = []
#         coord_RA_y = []
#         
#         # get coordinates of HVC(RA) targets
#         for target_ind in RA_super_targets[neuron]:
#             source_neuron.append(neuron)
#             synapse_distance.append(d(x_RA[neuron], y_RA[neuron], x_RA[target_ind], y_RA[target_ind])/SIDE)
#             coord_RA_x.append(x_RA[target_ind])
#             coord_RA_y.append(y_RA[target_ind])
#             
# 
#         source_neurons.append(source_neuron)
#         synapse_distances.append(synapse_distance)
#         
#         # calculate mean and variance of synapse distances
#         mean_synapse_distances.append(np.mean(synapse_distance))
#         std_synapse_distances.append(np.std(synapse_distance))
#         
#         target_RA_x.append(coord_RA_x)
#         target_RA_y.append(coord_RA_y)
#         
#         coord_I_x = []
#         coord_I_y = []
#         
#          # get coordinates of HVC(I) targets
#         for target_ind in RA_targets[neuron]:
#             coord_I_x.append(x_I[target_ind])
#             coord_I_y.append(y_I[target_ind])
#             
#         target_I_x.append(coord_I_x)
#         target_I_y.append(coord_I_y)
#         
#==============================================================================
            
numBins = 30

f1 = plt.figure()
ax = f1.add_subplot(311)
ax.hist(distances_between_RA_and_their_I_targets, numBins)
#ax.set_xlabel("Normalized distance")
ax.set_ylabel("Number of connections")
ax.set_title("HVC(RA) to HVC(I)")
ax.set_xlim([0,10])

numBins = 25


ax = f1.add_subplot(312)
ax.hist(distances_between_I_and_their_RA_targets, numBins)
#ax.set_xlabel("Normalized distance")
ax.set_ylabel("Number of connections")
ax.set_title("HVC(I) to HVC(RA)")
ax.set_xlim([0,10])


numBins = 35

ax = f1.add_subplot(313)
ax.hist(distances_between_RA_and_their_RA_targets, numBins)
ax.set_xlabel("Normalized distance")
ax.set_ylabel("Number of connections")
ax.set_title("HVC(RA) to HVC(RA)")
ax.set_xlim([0,10])

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
#==============================================================================
# print "number of neurons with supersynaptic targets: ",len(source_neurons)
# print source_neurons
# 
# num_neurons_in_plot = 16
# 
# num_cols = 4
# num_rows = 4
# 
# num_plots = len(source_neurons) / num_neurons_in_plot # number of fully covered plots
# 
# num_neurons_last_plot = len(source_neurons) - num_plots * num_neurons_in_plot
# num_rows_last_plot = num_neurons_last_plot / num_rows + 1
# num_cols_last_plot = num_neurons_last_plot % 4
# 
# 
# print "num full plots: ",num_plots
# print "num rows in last plot: ",num_rows_last_plot
# print "num cols in last plot: ",num_cols_last_plot
# 
# 
# 
# # plot full plots first
# for i in range(num_plots):    
#     f, axarr = plt.subplots(num_rows, num_cols)
#     start_id = i*num_neurons_in_plot
#     
#     
#     for j in range(num_rows):
#         for k in range(num_cols):
#             source_id = source_neurons[start_id+j*num_cols+k][0]
#             
#             for m in range(len(target_RA_x[start_id+j*num_cols+k])):
#                 
#                 axarr[j,k].plot([x_RA[source_id], target_RA_x[start_id+j*num_cols+k][m]], [y_RA[source_id], target_RA_y[start_id+j*num_cols+k][m]], color='red')            
#             
#             axarr[j,k].scatter(target_I_x[start_id+j*num_cols+k], target_I_y[start_id+j*num_cols+k], color='red', marker='x', s=50)            
#             
#             #axarr[j,k].scatter(target_x[start_id+j*num_cols+k], target_y[start_id+j*num_cols+k], marker='x', s=50)
#             axarr[j,k].set_xlim([0, SQUARE_SIDE])
#             axarr[j,k].set_ylim([0, SQUARE_SIDE])
# 
# # plotting the last figure
# f, axarr = plt.subplots(num_rows, num_cols)
# start_id = num_plots*num_neurons_in_plot
# 
# for j in range(num_rows_last_plot):
#     for k in range(num_cols):
#         
#         if j == num_rows_last_plot - 1 and k >= num_cols_last_plot:
#             break
#         
#         source_id = source_neurons[start_id+j*num_cols+k][0]
#         
#         for m in range(len(target_RA_x[start_id+j*num_cols+k])):
#             
#             axarr[j,k].plot([x_RA[source_id], target_RA_x[start_id+j*num_cols+k][m]], [y_RA[source_id], target_RA_y[start_id+j*num_cols+k][m]], color='red')            
#         
#         axarr[j,k].scatter(target_I_x[start_id+j*num_cols+k], target_I_y[start_id+j*num_cols+k], color='red', marker='x', s=50)            
#         
#         #axarr[j,k].scatter(target_x[start_id+j*num_cols+k], target_y[start_id+j*num_cols+k], marker='x', s=50)
#         axarr[j,k].set_xlim([0, SQUARE_SIDE])
#         axarr[j,k].set_ylim([0, SQUARE_SIDE])
#==============================================================================
             
#ax.set_xlabel("synapse distance")
#ax.set_ylabel("neuron id")
#ax.set_title("Distribution of output synaptic distances for each neuron")


plt.show()