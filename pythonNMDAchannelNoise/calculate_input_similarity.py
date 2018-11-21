#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 15:17:27 2018

@author: jingroup

Script find Jaccard index for neurons that burst synchronously
"""
import reading
import matplotlib.pyplot as plt
import numpy as np
import os
import utils

CURRENT_INJECTION = 50.0

def find_inputs(connections):
    """
    Find input connections for all neurons
    """
    N = len(connections)
    
    inputs = [set() for i in range(N)]
    
    for i in range(N):
        for target in connections[i]:
            inputs[target].add(i)
            
    return inputs


def get_sorted_burstOnsets_mature(filenameSpikes, filenameMature):
    """
    Get sorted burst onset times of mature neurons
    """
    (_, _, mature_indicators) = reading.read_mature_indicators(filenameMature)
    spike_times_s,  neuron_id_s,  ordered_spikes_s, neuron_ordered_id_s = utils.getSpikes(filenameSpikes)
    
    
    mature_neurons = set(np.where(mature_indicators == 1)[0])    
    
    first_spike_times_mature = []
    id_mature_spiked = []
    
    for spikes, neuron_id in zip(spike_times_s, neuron_id_s):
        if neuron_id[0] in mature_neurons:
            first_spike_times_mature.append(spikes[0])
            id_mature_spiked.append(neuron_id[0])
    
    
    first_spike_times_mature = np.array(first_spike_times_mature)
    id_mature_spiked = np.array(id_mature_spiked)
    
    min_first_spike_time_mature = np.min(first_spike_times_mature)
    
    first_spike_times_mature -= min_first_spike_time_mature
    
    
    ind_sorted = np.argsort(first_spike_times_mature)
        
    return first_spike_times_mature[ind_sorted], id_mature_spiked[ind_sorted]
    

def find_average_Jaccard_index(fileConnections, burst_times_sorted, ind_sorted_bursted, window):
    """
    Find average Jaccard index for network
    """
    (_, _, active_synapses) = reading.read_synapses(fileConnections)


    inputs = find_inputs(active_synapses)
    
    
    
    mean_Jaccard_window = []

    for i, ind in enumerate(ind_sorted_bursted):
        burst_time_neuron = burst_times_sorted[i]
        
        ind_neib = np.where((burst_times_sorted >= burst_time_neuron - window) & (burst_times_sorted <= burst_time_neuron + window))[0]
        
        Jaccard = []
        num_inputs = len(inputs[ind])
        
        if num_inputs != 0:
            for jj in ind_neib:
                inter = inputs[ind].intersection(inputs[ind_sorted_bursted[jj]])
                union = inputs[ind].union(inputs[ind_sorted_bursted[jj]])
                
                Jaccard.append(float(len(inter)) / float(len(union)))
        
            mean_Jaccard_window.append(np.mean(Jaccard))

    return np.mean(mean_Jaccard_window)


data_dir = "/mnt/hodgkin/eugene/results/immature/clusters/"
test_dir = "/mnt/hodgkin/eugene/results/immature/clusters/test/"


simnames = ["matTrans68", "matTrans62", "matTrans69", "matTrans66", "matTrans64"]
trials = [37400, 58400, 27800, 71400, 57200]
velocities = [0.5, 1.0, 1.33, 2.0, 10.0]

############################################
## Calculate Jaccard indices for networks ##
############################################
windows = np.array([0.5*i for i in range(40)]) # was 25

mean_Jaccard = []

for simname, trial in zip(simnames, trials):
    print "simname = {0}; trial = {1}".format(simname, trial)
    
    fileConnections = os.path.join(data_dir, simname + "/RA_RA_active_connections_"+str(trial)+".bin")
    fileMature = os.path.join(data_dir, simname + "/mature_"+str(trial)+".bin")
    
    fileSomaSpikes = os.path.join(test_dir, simname + "/trial" + str(trial) + "/test_spike_times_soma_5.bin")
    
    burst_times_sorted, ind_sorted_bursted = get_sorted_burstOnsets_mature(fileSomaSpikes, fileMature)
    
    mean_Jaccard_index = np.empty(len(windows), np.float32)
    
    for i, window in enumerate(windows):
        print window
        mean_Jaccard_index[i] = find_average_Jaccard_index(fileConnections, burst_times_sorted, ind_sorted_bursted, window)

    mean_Jaccard.append(mean_Jaccard_index)

f = plt.figure()
ax = f.add_subplot(111)
for i, v in enumerate(velocities):
    ax.plot(windows, mean_Jaccard[i], label="velocity {0}x".format(v))
ax.set_xlabel('Window (ms)')
ax.set_ylabel('Similarity between input connections')
plt.legend()
plt.show()

#outfile = os.path.join(prefix_dir+network_dir, "Jaccard.npz")
#np.savez(outfile, windows=windows, mean_Jaccard_index=mean_Jaccard_index)
