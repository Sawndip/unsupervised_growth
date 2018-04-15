# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 09:37:24 2018

@author: jingroup

Script plots input weights to certain neurons
"""
import reading
import matplotlib.pyplot as plt
import utils
import os

dirname = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/noImmatureOut7/"
trial_number = 2000
trialStep = 50

fileAxonalDelaysRA2RA = os.path.join(dirname, "axonal_delays_RA2RA_" + str(trial_number) + ".bin")
(_, _, axonal_delays_RA2RA) = reading.read_axonal_delays(fileAxonalDelaysRA2RA)
# 

#source_neurons = [289, 803, 548, 293, 801, 73, 402, 210, 820, 856, 793, 186, 27, 29]
#target_neurons = [927]

#==============================================================================
# input_weights = utils.get_source_input_weight_time_sequence(dirname, trial_number, trialStep, source_neurons, target_neurons)
# 
# num_targets = len(target_neurons)
# 
# time = [i*trialStep for i in range(trial_number / trialStep +1)]
# 
# f, axarr = plt.subplots(num_targets, sharex=True)
# #axarr[i].set_ylim([-0.1,1.1])
# 
# for i in range(num_targets):
#     axarr[i].plot(time, input_weights[i], label='target {0}'.format(target_neurons[i]))
#     axarr[i].legend(loc=2)
#     
# axarr[num_targets - 1].set_xlabel('time in # of trials')
#     
#==============================================================================

#==============================================================================
# 
# ###############################################################
# ### Input weights from source_neurons to a single target ######
# ###############################################################
# source_neurons = [195, 228, 860, 461, 525, 146, 726, 873, 252, 893, 800, 643, 612, 293, 806, 460, 108, 941, 784, 433, 771, 507]
# target_neuron = 293
# 
# first_group = set([195, 228, 860, 461, 525, 146, 726, 873, 252, 893])
# second_group = set([800, 643, 612, 293, 806, 460, 108, 941, 784, 433, 565, 771, 507])
# 
# 
# synaptic_weights = utils.get_weight_time_sequence(dirname, trial_number, trialStep, source_neurons, target_neuron)
# maturation_indicators = utils.get_maturation_time_sequence(dirname, trial_number, trialStep, [target_neuron])
# 
# num_source = len(source_neurons)
# 
# time = [i*trialStep for i in range(trial_number / trialStep +1)]
# 
# f, axarr = plt.subplots(num_source + 1, sharex=True)
# #axarr[i].set_ylim([-0.1,1.1])
# 
# for i in range(num_source):
#     if source_neurons[i] in first_group:
#         label = '1st g source {0}; delay = {1:.1f}'.format(source_neurons[i], axonal_delays_RA2RA[source_neurons[i]][target_neuron])
#     elif source_neurons[i] in second_group:
#         label = '2nd g source {0}; delay = {1:.1f}'.format(source_neurons[i], axonal_delays_RA2RA[source_neurons[i]][target_neuron])
#     else:
#         label = 'source {0}; delay = {1:.1f}'.format(source_neurons[i], axonal_delays_RA2RA[source_neurons[i]][target_neuron])
#     
#     axarr[i].plot(time, synaptic_weights[i], label=label)
#     axarr[i].set_ylabel('w')
#     axarr[i].legend(loc=2)
# 
# axarr[num_source].plot(time, maturation_indicators[0])    
# axarr[num_source].set_xlabel('time in # of trials')
# axarr[num_source].set_ylabel('mature')
# axarr[num_source].set_ylim([-0.1, 1.1])
#==============================================================================


 ##########################################################
 ## Connection from single source to target_neurons
 ##########################################################
N_RA = 1000
 
source_neuron = 228
 
 ##### find all supersynaptic targets ever existing from source neuron
all_supersynapse_targets = set()
 
current_trial = 0
timepoint = 0
 
while current_trial <= trial_number:   
    fileSuper = os.path.join(dirname, "RA_RA_super_connections_" + str(current_trial) + ".bin")
     
    (N, _, super_synapses) = reading.read_synapses(fileSuper)
     
    for target_id in super_synapses[source_neuron]:
        all_supersynapse_targets.add(target_id)
             
    timepoint += 1
    current_trial += trialStep
 
target_neurons = list(all_supersynapse_targets)
print "Target neurons: ",target_neurons
 

for target in target_neurons:
    print "{0} -> {1} delay = {2}".format(source_neuron, target, axonal_delays_RA2RA[source_neuron][target])
 
synaptic_weights = utils.get_weight_time_sequence_from_source(dirname, trial_number, trialStep, source_neuron, target_neurons)
maturation_indicators = utils.get_maturation_time_sequence(dirname, trial_number, trialStep, target_neurons)
#spike_times_sequence = utils.get_first_spike_time_sequence_relative_to_training(dirname, trial_number, trialStep, N_RA, target_neurons, [source_neuron])
 
stdp_timeDif_sequence = utils.get_stdp_timeDif_sequence(dirname, trial_number, trialStep, N_RA, target_neurons, source_neuron)
 
 
 
num_target = len(target_neurons)
 
time = [i*trialStep for i in range(trial_number / trialStep +1)]
 
f, axarr = plt.subplots(num_target, sharex=True)
#axarr[i].set_ylim([-0.1,1.1])
 
for i in range(num_target):
    label = 'id {0} delay {1:.1f} ms'.format(target_neurons[i], axonal_delays_RA2RA[source_neuron][target_neurons[i]])
     
    axarr[i].plot(time, synaptic_weights[i], label=label)
    axarr[i].set_ylabel('w')
    axarr[i].tick_params('y', color='b')
    axarr[i].set_ylabel('w', color='b')
    axarr[i].legend(loc=2)
     
       
    ax = axarr[i].twinx()
    ax.plot(time, stdp_timeDif_sequence[i], color='r')    
    ax.tick_params('y', color='r')
    ax.set_ylabel('$\Delta t$', color='r')
    ax.set_ylim([5, 25])
 
   
    
plt.show()