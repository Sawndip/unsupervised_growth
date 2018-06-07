# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 11:28:55 2017

@author: jingroup

Script analyzes deliveries relative to the birst onset: onset of first spike
"""

import reading
import numpy as np
import os
import matplotlib.pyplot as plt

#dirname = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/maturationTransition4/"
dirname = "/home/eugene/results/immature/clusters/matTrans44/"

trial_number = 20400

fileConnections = os.path.join(dirname, "RA_RA_active_connections_" + str(trial_number) + ".bin")
#fileSpikeTimes = os.path.join(dirname, "spike_times_soma_" + str(trial_number) + ".bin")
fileDelaysRA2RA = os.path.join(dirname, "axonal_delays_RA2RA_" + str(trial_number) + ".bin")
fileMature = os.path.join(dirname, "mature_" + str(trial_number) + ".bin")
fileWeights = os.path.join(dirname, "weights_" + str(trial_number) + ".bin")
fileTraining = os.path.join(dirname, "training_neurons.bin")


#fileSpikeTimes = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/test/maturationTransition4/test_spike_times_soma_7.bin"
fileSpikeTimes = "/home/eugene/results/immature/clusters/test/matTrans44/test_spike_times_soma_10.bin"

MARGIN_LATE = 0.0 # margin for the burst coming late

(N_RA, _, active_synapses) = reading.read_synapses(fileConnections)
(_, _, spike_times_raw, neuron_spiked) = reading.read_time_info(fileSpikeTimes)
(_, _, axonal_delays_RA2RA) = reading.read_axonal_delays(fileDelaysRA2RA)
(_, _, weights) = reading.read_weights(fileWeights)
(_, _, mature_indicators) = reading.read_mature_indicators(fileMature)
training_neurons = reading.read_training_neurons(fileTraining)

print "Number of HVC(RA) neurons = ",N_RA

first_spike_times = np.empty(N_RA, np.float32)
first_spike_times.fill(-1.0)

for n, time in zip(neuron_spiked, spike_times_raw):
    first_spike_times[n[0]] = time[0]


delivered_times_and_id = []
arrivals = []

#for i in range(N_RA):
#    mature_indicators[i] = 1

#==============================================================================
# num_inputs = np.zeros(N_RA, np.int32)
# G_input = np.zeros(N_RA, np.float32)
# 
# for i in range(N_RA):
#     for target, G in zip(RA2RA_targets_ID[i], RA2RA_syn_G[i]):
#         num_inputs[target] += 1
#         G_input[target] += G
# 
#==============================================================================
#==============================================================================
# print "Source: "        
# print first_spike_times[5452]
# print num_inputs[5452]
# print G_input[5452]    
# print 5650 in RA2RA_targets_ID[5452]
# print axonal_delays_RA2RA[5452][RA2RA_targets_ID[5452].index(5650)]
# 
# print "Target: "
# print first_spike_times[5650]
# print num_inputs[5650]
# print G_input[5650]    
#==============================================================================

example_neurons = [548, 293, 210, 856, 793, 186, 27, 29, 927]
num_example = len(example_neurons)


example_arrivals = [[] for i in range(len(example_neurons))]

for i in range(N_RA):
    if first_spike_times[i] > 0:
        for j, target in enumerate(active_synapses[i]):
            if target not in training_neurons and first_spike_times[target] > 0 and mature_indicators[target] == 1 and mature_indicators[i] == 1:
                time_difference = first_spike_times[i] + axonal_delays_RA2RA[i][target] - first_spike_times[target]
                
                if time_difference > 20:
                    print "time difference = {0} source: {1} target: {2} delay = {3} weight = {4}".format(time_difference, i, target, axonal_delays_RA2RA[i][target], weights[i][target])
                else:
                    arrivals.append(time_difference)
                
                if target in example_neurons:
                    ind = example_neurons.index(target)
                    example_arrivals[ind].append(time_difference)


num_arrivals = len(arrivals)
num_late_arrivals = sum(arrival > MARGIN_LATE for arrival in arrivals)

print "Total number of deliveries: ",num_arrivals
print "Number of late deliveries: ",num_late_arrivals
print "Fraction of late deliveries: ",float(num_late_arrivals) / float(num_arrivals)


plt.figure()

plt.hist(arrivals, bins = 200)
#plt.title('Statistics of deliveries')
plt.xlabel('spike delivery time - target burst onset time (ms)')
plt.ylabel('# of deliveries')
plt.xlim([-15,15])


f, axarr = plt.subplots(num_example)

for i in range(num_example):
    axarr[i].hist(example_arrivals[i], bins = 10, label="{0}".format(example_neurons[i]))
    axarr[i].legend()
#plt.title('Statistics of deliveries')
plt.xlabel('spike delivery time - target burst onset time (ms)')
plt.ylabel('# of deliveries')
plt.xlim([-15,15])
plt.show()


#(_, simulation_time, burst_times_raw, neuron_fired) = reading.read_time_info(fileSpikePattern)
#burst_labels = reading.read_burst_labels(fileBurstLabels)

#print delivered_times_and_id

#print burst_times
#print neuron_fired
#N_RA = 6000



