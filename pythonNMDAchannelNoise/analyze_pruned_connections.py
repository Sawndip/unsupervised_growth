#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 18:18:54 2018

@author: jingroup

Script analyzes which connections were pruned during simulation
"""
import reading
import os
import matplotlib.pyplot as plt
import numpy as np

#dirname = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/maturationTransition2/"
dirname = "/home/eugene/results/immature/clusters/matTrans19/"
end_trial = 18000
trialStep = 100



# go through all previous trials and record all supersynaptic connections
num_timepoints = end_trial / trialStep + 1
    
all_supersynapses_ever = set()

current_trial = 0
timepoint = 0

while current_trial <= end_trial:   
    fileSuper = os.path.join(dirname, "RA_RA_super_connections_" + str(current_trial) + ".bin")
    
    (N, _, super_synapses) = reading.read_synapses(fileSuper)
    
    for source_id in range(N):
        for target_id in super_synapses[source_id]:
            all_supersynapses_ever.add((source_id, target_id))
            
    timepoint += 1
    current_trial += trialStep
    

# check which supersynapses are not supersynapses at the end_trial

fileAxonalDelaysRA2RA = os.path.join(dirname, "axonal_delays_RA2RA_" + str(end_trial) + ".bin")
fileSuper = os.path.join(dirname, "RA_RA_super_connections_" + str(end_trial) + ".bin")

(N, _, super_synapses) = reading.read_synapses(fileSuper)
(_, _, axonal_delays_RA2RA) = reading.read_axonal_delays(fileAxonalDelaysRA2RA)

supersynapses_final = set()
delays_final = []

for i in range(N):
    for target in super_synapses[i]:
        supersynapses_final.add((i, target))
        delays_final.append(axonal_delays_RA2RA[i][target])
        
pruned_super = set()
pruned_delays = []
        
for (source_id, target_id) in all_supersynapses_ever:
    if (source_id, target_id) not in supersynapses_final:
        pruned_super.add((source_id, target_id))
        pruned_delays.append(axonal_delays_RA2RA[source_id][target_id])
        
#print pruned_super

all_axonal_delays_RA2RA = [delay for delays in axonal_delays_RA2RA for delay in delays]


plt.figure()
hist, bin_edges = np.histogram(all_axonal_delays_RA2RA, bins=50)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
plt.step(bin_centers, hist / float(np.sum(hist)), label="HVC(RA) -> HVC(RA) random", where="pre")

hist, bin_edges = np.histogram(pruned_delays, bins=50)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
plt.step(bin_centers, hist / float(np.sum(hist)), label="pruned super", where="pre")

hist, bin_edges = np.histogram(delays_final, bins=50)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
plt.step(bin_centers, hist / float(np.sum(hist)), label="super", where="pre")


plt.xlabel('Delay (ms)')
plt.ylabel("Counts")
plt.legend()


##############################################################################
#### Check if any pruning of supersynapses happened before target matured ####
##############################################################################

#==============================================================================
# # first find trial when neurons became mature
# num_timepoints = end_trial / trialStep + 1
#     
# current_trial = 0
# timepoint = 0
# 
# maturation_time = np.zeros(N)
# maturation_time[:] = np.nan
# #print maturation_time
# 
# while current_trial <= end_trial:   
#    fileMature = os.path.join(dirname, "mature_" + str(current_trial) + ".bin")
#    
#    (_, _, mature_indicators) = reading.read_mature_indicators(fileMature)
#     
#    for i in range(N):
#         if np.isnan(maturation_time[i]) and mature_indicators[i] == 1: 
#             maturation_time[i] = current_trial
#             
#    timepoint += 1
#    current_trial += trialStep
# 
# 
# fileTraining = os.path.join(dirname, "training_neurons.bin")
# 
# #print maturation_time
# 
# # record supersynapses that emerged and disappeared before target maturation
# all_supersynapses_emerged_before_target_maturation = set()
# 
# 
# current_trial = 0
# timepoint = 0
# while current_trial <= end_trial:   
#     fileSuper = os.path.join(dirname, "RA_RA_super_connections_" + str(current_trial) + ".bin")
#     
#     (N, _, super_synapses) = reading.read_synapses(fileSuper)
#     
#     for source_id in range(N):
#         for target_id in super_synapses[source_id]:
#             if (not np.isnan(maturation_time[target_id])) and (current_trial >= maturation_time[target_id]):
#                 all_supersynapses_emerged_before_target_maturation.add((source_id, target_id))
#             
#     timepoint += 1
#     current_trial += trialStep
# 
# supersynapses_before_target_maturation_final = set()
# delays_before_target_maturation_final = []
# 
# for target_id in range(N):
#     if not np.isnan(maturation_time[target_id]):
#         fileSuper = os.path.join(dirname, "RA_RA_super_connections_" + str(int(maturation_time[target_id])) + ".bin")
#     else:
#         fileSuper = os.path.join(dirname, "RA_RA_super_connections_" + str(end_trial) + ".bin")
#     
#     (N, _, super_synapses) = reading.read_synapses(fileSuper)
#     
#     for source_id in range(N):
#         if target_id in super_synapses[source_id]:
#             supersynapses_before_target_maturation_final.add((source_id, target_id))
#             delays_before_target_maturation_final.append(axonal_delays_RA2RA[source_id][target_id])
#             
# pruned_before_target_maturation_super = set()
# pruned_before_target_maturation_delays = []
# 
# 
#         
# for (source_id, target_id) in all_supersynapses_emerged_before_target_maturation:
#     if (source_id, target_id) not in supersynapses_before_target_maturation_final:
#         pruned_before_target_maturation_super.add((source_id, target_id))
#         pruned_before_target_maturation_delays.append(axonal_delays_RA2RA[source_id][target_id])
# 
# 
#==============================================================================

# record supersynapses that emerged and disappeared before target maturation

was_super_indicators = np.zeros((N, N), dtype=bool)

num_timepoints = end_trial / trialStep + 1
current_trial = 0
timepoint = 0

pruned_before_target_maturation_super = set()
pruned_before_target_maturation_delays = []

while current_trial <= end_trial:   
    fileSuper = os.path.join(dirname, "RA_RA_super_connections_" + str(current_trial) + ".bin")
    fileMature = os.path.join(dirname, "mature_" + str(current_trial) + ".bin")
   
    (_, _, mature_indicators) = reading.read_mature_indicators(fileMature)    
    (N, _, super_synapses) = reading.read_synapses(fileSuper)
    
    super_indicators = np.zeros((N, N), dtype=bool)
    
    super_indicators[:, np.where(mature_indicators == 1)] = 1
    
    
    for source_id in range(N):
        for target_id in super_synapses[source_id]:
            if mature_indicators[target_id] == 0:
                was_super_indicators[source_id][target_id] = 1
                super_indicators[source_id][target_id] = 1
    
    pruned_ind = np.where((was_super_indicators == 1) & (super_indicators == 0))
    
    #print pruned_ind
    
    num_pruned = len(pruned_ind[0])
    
    for i in range(num_pruned):
        row = pruned_ind[0][i]
        col = pruned_ind[1][i]
        pruned_before_target_maturation_super.add((row, col))
        pruned_before_target_maturation_delays.append(axonal_delays_RA2RA[row][col])
   
    
  #  for source_id in range(N):
   #     for target_id in range(N):
    #        if was_super_indicators[source_id][target_id] == 1 and super_indicators[source_id][target_id] == 0:
     #          pruned_before_target_maturation_super.add((source_id, target_id))
      #         pruned_before_target_maturation_delays.append(axonal_delays_RA2RA[source_id][target_id])
   
    timepoint += 1
    current_trial += trialStep


print pruned_before_target_maturation_delays
  
#print pruned_super
plt.figure()
plt.title('Before target maturation')
hist, bin_edges = np.histogram(all_axonal_delays_RA2RA, bins=50)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
plt.step(bin_centers, hist / float(np.sum(hist)), label="HVC(RA) -> HVC(RA) random", where="pre")

hist, bin_edges = np.histogram(pruned_before_target_maturation_delays, bins=50)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
plt.step(bin_centers, hist / float(np.sum(hist)), label="pruned super", where="pre")

#hist, bin_edges = np.histogram(delays_before_target_maturation_final, bins=50)
#bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
#plt.step(bin_centers, hist / float(np.sum(hist)), label="super", where="pre")


plt.xlabel('Delay (ms)')
plt.ylabel("Counts")
plt.legend()

plt.show()