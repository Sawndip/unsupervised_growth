# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 09:36:48 2018

@author: jingroup

Script analyzes excitatory input to neurons in the chain
"""
import matplotlib.pyplot as plt
import reading
import os
import numpy as np
import utils


def get_hist_for_discrete_integers(data):
    """
    Outputs histogram for discrete integer values
    """
    left_of_first_bin = data.min() - 0.5
    right_of_last_bin = data.max() + 0.5

    return 1.0, left_of_first_bin, right_of_last_bin


CONVERTION_CONSTANT = 10.0

trial_number = 13600
#dirname = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/noImmatureOut8/"
dirname = "/home/eugene/results/immature/clusters/11/"

fileSoma = os.path.join(dirname, "spike_times_soma_" + str(trial_number) + ".bin")
fileWeights = os.path.join(dirname, "weights_" + str(trial_number) + ".bin")
fileActiveSynapses = os.path.join(dirname, "RA_RA_active_connections_" + str(trial_number) + ".bin")
fileSuperSynapses = os.path.join(dirname, "RA_RA_super_connections_" + str(trial_number) + ".bin")
fileRemodeled = os.path.join(dirname, "remodeled_indicators_" + str(trial_number) + ".bin")
fileMature = os.path.join(dirname, "mature_" + str(trial_number) + ".bin")
fileTraining = os.path.join(dirname, "training_neurons.bin")

####################################
### Find connected neurons in chain
####################################
(N_RA, _, remodeled_indicators) = reading.read_remodeled_indicators(fileRemodeled)
(_, _, mature_indicators) = reading.read_mature_indicators(fileMature)
(_, _, active_synapses) = reading.read_synapses(fileActiveSynapses)
(_, _, super_synapses) = reading.read_synapses(fileSuperSynapses)


####### Calculate number of active and super inputs and outputs ############
num_active_inputs = np.zeros(N_RA, np.int32)
num_super_inputs = np.zeros(N_RA, np.int32)
num_active_outputs = np.zeros(N_RA, np.int32)
num_super_outputs = np.zeros(N_RA, np.int32)

for i in range(N_RA):
    num_active_outputs[i] = len(active_synapses[i])
    num_super_outputs[i] = len(super_synapses[i])

    for target in active_synapses[i]:
        num_active_inputs[target] += 1
    
    for target in super_synapses[i]:
        num_super_inputs[target] += 1
        
mature_neurons = np.where(mature_indicators == 1)[0]

print "Mature neurons: ", ','.join(map(str, mature_neurons)) 

# select mature neurons
neurons_in_chain = np.where(mature_indicators == 1)[0]

#==============================================================================
# training_neurons = reading.read_training_neurons(fileTraining)
# 
#==============================================================================
#==============================================================================
# neurons_in_chain = [i for i in training_neurons]
# 
# 
# 
# # start with training neurons
# # follow supersynapses
# 
# i = 0
# while i < len(neurons_in_chain):
#     if len(super_synapses[neurons_in_chain[i]]) > 0:
#         for j in super_synapses[neurons_in_chain[i]]:
#             if j not in neurons_in_chain:
#                 neurons_in_chain.append(j)
#         
#     i += 1
#     
#==============================================================================
print "Neurons in chain: ",neurons_in_chain



(trial_number, simulation_time, spike_times_soma, neuron_fired_soma) = reading.read_time_info(fileSoma)

set_neurons_in_chain = set(neurons_in_chain)

#print "Dedritic spikes: ", neuron_fired_dend
#print "Dendritic spike times: ", spike_times_dend

#print "Somatic spikes: ", neuron_fired_soma
#print "Somatic spike times: ", spike_times_soma

####################################
### Order neurons by their first spike time
####################################
ordered_soma_spikes_raw, ordered_soma_raw = zip(*sorted(zip(spike_times_soma, neuron_fired_soma)))

#print ordered_soma_spikes_raw
#print ordered_soma_raw

ordered_first_spike_times = [spikes[0] for spikes in ordered_soma_spikes_raw]
ordered_by_first_spikes = [ordered[0] for ordered in ordered_soma_raw]

#print ordered_by_first_spikes


####################################
### Calculate input weights to neurons
####################################
(_, _, weights) = reading.read_weights(fileWeights) 

training_neurons = reading.read_training_neurons(fileTraining)

super_targets = []
### Show super synapses for training neurons ###
for i in training_neurons:
    for synapse in super_synapses[i]:
        print "super {0} -> {1}; w = {2}".format(i, synapse, weights[i][synapse])
        super_targets.append(synapse)

from collections import Counter
counter = Counter(super_targets)

convergence = [i for i in counter.values()]

nbins = 10
f = plt.figure()

ax1 = f.add_subplot(111)
hist, bin_edges = np.histogram(convergence, bins=nbins)
#print bin_edges
#print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
ax1.step(bin_centers, hist, where="pre")
ax1.set_xlabel('# of convergent super inputs')
ax1.set_ylabel('# of neurons')



first_spike_time_of_chain_neurons = []
input_weights_of_chain_neurons = []

num_active_inputs_of_chain_neurons = []
num_super_inputs_of_chain_neurons = []
num_active_outputs_of_chain_neurons = []
num_super_outputs_of_chain_neurons = []

for neuron, spike_time in zip(ordered_by_first_spikes, ordered_first_spike_times):
    if neuron in set_neurons_in_chain:
        first_spike_time_of_chain_neurons.append(spike_time)
        input_weights_of_chain_neurons.append(CONVERTION_CONSTANT*np.sum(weights[:,neuron]/1000.0, axis=0))

        num_active_inputs_of_chain_neurons.append(num_active_inputs[neuron])
        num_active_outputs_of_chain_neurons.append(num_active_outputs[neuron])
        num_super_inputs_of_chain_neurons.append(num_super_inputs[neuron])
        num_super_outputs_of_chain_neurons.append(num_super_outputs[neuron])


#print first_spike_time_of_chain_neurons
#print input_weights_of_chain_neurons

f = plt.figure()

ax1 = f.add_subplot(511)
ax1.plot(first_spike_time_of_chain_neurons, input_weights_of_chain_neurons)
ax1.set_ylabel('In exc w (nS)')
ax1.set_xticklabels([])

ax2 = f.add_subplot(512)
ax2.plot(first_spike_time_of_chain_neurons, num_active_inputs_of_chain_neurons)
ax2.set_ylabel('# of act in')
ax2.set_xticklabels([])

ax3 = f.add_subplot(513)
ax3.plot(first_spike_time_of_chain_neurons, num_super_inputs_of_chain_neurons)
ax3.set_ylabel('# of sup in')
ax3.set_xticklabels([])

ax4 = f.add_subplot(514)
ax4.plot(first_spike_time_of_chain_neurons, num_active_outputs_of_chain_neurons)
ax4.set_ylabel('# of act out')
ax4.set_xticklabels([])

ax5 = f.add_subplot(515)
ax5.plot(first_spike_time_of_chain_neurons, num_super_outputs_of_chain_neurons)
ax5.set_ylabel('# of sup out')

ax5.set_xlabel('Time (ms)')


###### Plot input weight histogram ###############
nbins = 20

f = plt.figure()

ax1 = f.add_subplot(111)
hist, bin_edges = np.histogram(input_weights_of_chain_neurons, bins=nbins)
#print bin_edges
#print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
ax1.step(bin_centers, hist, where="pre")
ax1.set_xlabel('Input excitatory weight (nS)')
ax1.set_ylabel('# of neurons')

###### Plot in and out degree histograms ###########
nbins = 20

f = plt.figure()

ax1 = f.add_subplot(221)
#hist, bin_edges = np.histogram(num_active_inputs[neurons_in_chain], bins=nbins)
#print bin_edges
#print hist
#bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.

#ax1.step(bin_centers, hist, where="pre")
d, left_of_first_bin, right_of_last_bin = get_hist_for_discrete_integers(num_active_inputs[neurons_in_chain])
ax1.hist(num_active_inputs[neurons_in_chain], np.arange(left_of_first_bin, right_of_last_bin + d, d))
ax1.set_title("# active inputs")

ax2 = f.add_subplot(222)
#hist, bin_edges = np.histogram(num_super_inputs[neurons_in_chain], bins=nbins)
#print bin_edges
#print hist
#bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
#ax2.step(bin_centers, hist, where="pre")

d, left_of_first_bin, right_of_last_bin = get_hist_for_discrete_integers(num_super_inputs[neurons_in_chain])
ax2.hist(num_super_inputs[neurons_in_chain], np.arange(left_of_first_bin, right_of_last_bin + d, d))
ax2.set_title("# super inputs")

ax3 = f.add_subplot(223)
#hist, bin_edges = np.histogram(num_active_outputs[neurons_in_chain], bins=nbins)
#print bin_edges
#print hist
#bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
#ax3.step(bin_centers, hist, where="pre")
d, left_of_first_bin, right_of_last_bin = get_hist_for_discrete_integers(num_active_outputs[neurons_in_chain])
ax3.hist(num_active_outputs[neurons_in_chain], np.arange(left_of_first_bin, right_of_last_bin + d, d))
ax3.set_title("# active outputs")

ax4 = f.add_subplot(224)
#hist, bin_edges = np.histogram(num_super_outputs[neurons_in_chain], bins=nbins)
#print bin_edges
#print hist
#bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
#ax4.step(bin_centers, hist, where="pre")

d, left_of_first_bin, right_of_last_bin = get_hist_for_discrete_integers(num_super_outputs[neurons_in_chain])
ax4.hist(num_super_outputs[neurons_in_chain], np.arange(left_of_first_bin, right_of_last_bin + d, d))
ax4.set_title("# super outputs")


###### Plot total number of active and super synapses #####
trialStep = 50

total_num_active, total_num_super = utils.get_num_active_and_super_synapses(dirname, trial_number, trialStep)

f = plt.figure()

ax1 = f.add_subplot(211)
ax1.plot([float(i)*trialStep for i in range(trial_number / trialStep + 1)], total_num_active)
ax1.set_ylabel("# active synapses")

ax2 = f.add_subplot(212)
ax2.plot([float(i)*trialStep for i in range(trial_number / trialStep + 1)], total_num_super)
ax2.set_ylabel("# super synapses")
ax2.set_xlabel('time in trials')

layer = set(training_neurons)

print "Training_neurons:", layer
print "In degree active = ",num_active_inputs[list(layer)]
print "Out degree active = ",num_active_outputs[list(layer)]
print "In degree super = ",num_super_inputs[list(layer)]
print "Out degree super = ",num_super_outputs[list(layer)]

prev_layer = set()   

while len(layer) <= 10 and len(layer) > 0:
    num_active_out = 0
    num_super_out = 0
    
    new_layer = set()
    
    print "Current layer: ",layer
    print "Prev layer: ",prev_layer
        
    
    for i in layer:
        #print len(super_synapses[i])
        #print len(active_synapses[i])
        num_active_out += len(active_synapses[i])
        num_super_out += len(super_synapses[i])
        
        for synapse in super_synapses[i]:
            if synapse not in prev_layer:
                new_layer.add(synapse)
    
    print "Next layer: ",new_layer
    
        
    prev_layer = layer
    layer = new_layer
        
    
    print "Num active out = ",num_active_out
    print "Num super out = ",num_super_out
    print "In degree active = ",num_active_inputs[list(layer)]
    print "Out degree active = ",num_active_outputs[list(layer)]
    print "In degree super = ",num_super_inputs[list(layer)]
    print "Out degree super = ",num_super_outputs[list(layer)]
   
    
#for i in range(N_RA):
#    if 553 in active_synapses[i]:
#        print "active: {0} -> 553; w = {1}".format(i, weights[i][553])

    #if 553 in super_synapses[i]:
     #   print "super: {0} -> 553; w = {1}".format(i, weights[i][553])


source_neurons = range(N_RA)
#target_neurons = [448]
target_neurons = [897, 803, 465, 203, 913, 210, 186, 785, 125, 927]

for i in source_neurons:
    if i not in training_neurons:
        for j in target_neurons:
            if j in active_synapses[i]:
                print "active: {0} -> {1}; w = {2}".format(i, j, weights[i][j])

            if j in super_synapses[i]:
                print "super: {0} -> {1}; w = {2}".format(i, j, weights[i][j])


for i in range(N_RA):
    if 2 in active_synapses[i]:
        print "{0} -> 2".format(i)
    
    if 93 in active_synapses[i]:
        print "{0} -> 93".format(i)



plt.show()

#spike_times_soma = [t for sublist in list(spike_times_soma) for t in sublist]
#neuron_fired_soma = [ind for sublist in list(neuron_fired_soma) for ind in sublist]


#ordered_soma_spikes, ordered_neurons = zip(*sorted(zip(spike_times_soma, neuron_fired_soma)))

#ordered_soma_spikes = [t for sublist in list(ordered_soma_spikes_raw) for t in sublist]
#ordered_soma_neurons = [t for sublist in list(ordered_soma_raw) for t in sublist]