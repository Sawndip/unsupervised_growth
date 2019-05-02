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
import os

#dirname = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/maturationTransition4/"
dirname = "/home/eugene/results/immature/clusters/matTrans44/"

trial_number = 20400

fileActive = os.path.join(dirname, "RA_RA_active_connections_" + str(trial_number) + ".bin")
fileSuper = os.path.join(dirname, "RA_RA_super_connections_" + str(trial_number) + ".bin")
fileMature = os.path.join(dirname, "mature_" + str(trial_number) + ".bin")

fileSpikes = "/home/eugene/results/immature/clusters/test/matTrans44/test_spike_times_soma_10.bin"
#fileSpikes = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/test/maturationTransition4/test_spike_times_soma_7.bin"


(_,_,mature_indicators) = reading.read_mature_indicators(fileMature)

mature_neurons = np.where(mature_indicators == 1)[0]


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
    if mature_indicators[i] == 1 and all_first_spike_times[i] > 0:
        for target in super_synapses[i]:
            if mature_indicators[target] == 1 and all_first_spike_times[target] > 0:
                time_difference = all_first_spike_times[target] - all_first_spike_times[i]
                if time_difference < -50.0 or time_difference > 50.0:
                    print "time difference {0} for {1} -> {2}".format(time_difference, i, target)
                else:
                    dt.append(time_difference)
                    
print "Mean dt = ",np.mean(dt)
print "Std dt = ",np.std(dt)
                   
nbins = 50

plt.figure()
plt.hist(dt, bins=nbins)
plt.xlabel('Postsyn - presyn first spike time (ms)')
plt.ylabel('Counts')

plt.show()