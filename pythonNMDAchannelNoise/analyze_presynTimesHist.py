#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 12:00:05 2018

@author: jingroup

Script plots histogram of postsynaptic - presynaptic times
"""
import reading
import matplotlib.pyplot as plt
import numpy as np

fileSpikes = "/home/eugene/results/immature/clusters/test/matTrans19/test_spike_times_soma_1.bin"
fileActive = "/home/eugene/results/immature/clusters/matTrans19/RA_RA_active_connections_22000.bin"
fileSuper = "/home/eugene/results/immature/clusters/matTrans19/RA_RA_super_connections_22000.bin"


(trial_number, simulation_time, spike_times_raw, neuron_fired) = reading.read_time_info(fileSpikes)
(N_RA, _, active_synapses) = reading.read_synapses(fileActive)
(N_RA, _, super_synapses) = reading.read_synapses(fileSuper)

print "Number of HVC(RA) neurons = ",N_RA

#print neuron_fired

threshold = 1.0

all_first_spike_times = np.empty(N_RA, np.float32)
all_first_spike_times.fill(-100.0)

num_bursts = 0

for n, time in zip(neuron_fired, spike_times_raw):
    num_bursts += 1
    all_first_spike_times[n[0]] = time[0]

print list(all_first_spike_times[:200])

id_last_fired = np.argmax(all_first_spike_times)

print "Total number of bursts = ",num_bursts

print "Number of silent neurons: ",np.shape(np.where(all_first_spike_times[:(id_last_fired+1)] < 0)[0])[0]
print "id of last fired HVC(RA) = ",id_last_fired
print "Max burst time relative to current injection: ",np.max(all_first_spike_times)

#print burst_times

dt = []

for i in range(N_RA):
    if all_first_spike_times[i] > 0:
        for target in super_synapses[i]:
            if all_first_spike_times[target] > 0:
                dt.append(all_first_spike_times[target] - all_first_spike_times[i])
                    
print "Mean dt = ",np.mean(dt)
print "Std dt = ",np.std(dt)
                   
nbins = 50

plt.figure()
plt.hist(dt, bins=nbins)
plt.xlabel('Postsyn - presyn first spike time (ms)')
plt.ylabel('Counts')

plt.show()