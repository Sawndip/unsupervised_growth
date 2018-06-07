#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 21 11:51:44 2018

@author: jingroup

Script analyzes topology of grown network using both
connectivity and spike pattern
"""
import reading
import os
import numpy as np
import utils
import matplotlib.pyplot as plt
import networkx as nx

datadir = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/chain_0.5_delay_3.0"
#datadir = "/home/eugene/results/immature/clusters/matTrans31/"

#trial_number = 20600
trial_number = 0

fileConnections = os.path.join(datadir, "RA_RA_super_connections_" + str(trial_number) + ".bin")
fileDelaysRA2RA = os.path.join(datadir, "axonal_delays_RA2RA_" + str(trial_number) + ".bin")
fileMature = os.path.join(datadir, "mature_" + str(trial_number) + ".bin")
fileTraining = os.path.join(datadir, "training_neurons.bin")

#fileSpikeTimes = "/home/eugene/results/immature/clusters/test/matTrans40/test_spike_times_soma_10.bin"
#fileJitter = "/home/eugene/results/immature/clusters/test/matTrans31/jitter.bin"
fileJitter = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/chain_0.5_delay_3.0/jitter.bin"



(N_RA, _, super_synapses) = reading.read_synapses(fileConnections)
(_, _, mature_indicators) = reading.read_mature_indicators(fileMature)
training_neurons = reading.read_training_neurons(fileTraining)

# reading spike times of particular trial
#(_, _, spike_times_raw, neuron_spiked) = reading.read_time_info(fileSpikeTimes)

# reading mean spike times
_, num_test_trials, \
    probability_soma_spike, average_num_soma_spikes_in_trial, first_spike_times, std_first_soma_spike_time,\
    probability_dend_spike, average_num_dend_spikes_in_trial, mean_first_dend_spike_time, std_first_dend_spike_time = reading.read_jitter(fileJitter) 


print "Number of HVC(RA) neurons = ",N_RA

#==============================================================================
# first_spike_times = np.empty(N_RA, np.float32)
# first_spike_times.fill(-1.0)
# 
# for n, time in zip(neuron_spiked, spike_times_raw):
#     first_spike_times[n[0]] = time[0]
#==============================================================================

# sort all connections to achieve fast look-up
for i in range(N_RA):
    super_synapses[i].sort()

windows = [0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]

mean_fraction_shared = []
mean_fraction = []

for window in windows:
    fraction_shared_window = []
    fraction_window = []
    for i in range(N_RA):
        if (mature_indicators[i]) and (first_spike_times[i] > 0) and (i not in training_neurons):
            center_time = first_spike_times[i]
            synch_neurons = np.where((first_spike_times > center_time - window/2.) & (first_spike_times < center_time + window/2.))[0]
            synch_neurons = synch_neurons[mature_indicators[synch_neurons] == 1]
    
            #print synch_neurons
            #print "Center time: ",center_time
            #print first_spike_times[synch_neurons]
            
            all_inputs = set() # all inputs received by synchronous neurons 
            
            inputs_to_synch = []
            
            for n in synch_neurons:
                inputs = set()
                
                for k in range(N_RA):
                    try:
                        ind = utils.index(super_synapses[k], n)
                        inputs.add(k)
                        all_inputs.add(k)
                    except ValueError:
                        continue
                    
                inputs_to_synch.append(inputs)
            
            
            shared_inputs = set.intersection(*inputs_to_synch)
            
            fraction_of_all_inputs = [float(len(inp)) / float(len(all_inputs)) for inp in inputs_to_synch]
               
            fraction_shared_window.append(float(len(shared_inputs)) / float(len(all_inputs)))
            fraction_window.append(np.mean(fraction_of_all_inputs))
            #print inputs_to_synch
    print "Window: {0}; mean fraction shared: {1}; mean fraction: {2}".format(window, np.mean(fraction_shared_window), np.mean(fraction_window))
    #print fraction_shared_window
    #print fraction_window
    
    
    mean_fraction_shared.append(np.mean(fraction_shared_window))
    mean_fraction.append(np.mean(fraction_window))
   
plt.figure()
plt.plot(windows, mean_fraction_shared, '-ro', label='shared')
plt.plot(windows, mean_fraction, '-bo', label='individual')

plt.xlabel('Window size (ms)')
plt.ylabel('Mean fraction')

# make a directed graph and estimate distance to training neurons
DG = nx.DiGraph()

DG.add_nodes_from(range(N_RA))

for i in range(N_RA):
    edges = [(i,t) for t in super_synapses[i]]
    
    DG.add_edges_from(edges)

distances = np.zeros(N_RA, np.int32)
    
#distances = np.empty(N_RA, np.int32)
#distances.fill(np.nan)

for i in range(N_RA):
    if i not in training_neurons:
        smallest_distance = 10000
        for j in training_neurons:
            try:
                d = nx.shortest_path_length(DG, source=j, target=i)
                if d < smallest_distance:
                    smallest_distance = d
            except nx.NetworkXNoPath:
                continue

        if smallest_distance == 10000:
            print "No path to training neurons was found for neuron",i
        else:
            distances[i] = smallest_distance

indsorted = np.argsort(first_spike_times[(mature_indicators > 0) & (first_spike_times > 0)])

plt.figure()
plt.plot(first_spike_times[(mature_indicators > 0) & (first_spike_times > 0)][indsorted], distances[(mature_indicators > 0) & (first_spike_times > 0)][indsorted], '-o')
plt.ylabel('Distance to training neurons')
plt.xlabel('Burst time (ms)')

plt.show()