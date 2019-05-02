#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  8 23:12:40 2018

@author: jingroup

Script finds supersynaptic connections from immature neurons outside chain
to neurons in the chain
"""

import matplotlib.pyplot as plt
import reading
import os
import numpy as np
import utils



#CONVERTION_CONSTANT = 10.0

trial_number = 5500
#dirname = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/noImmatureOut8/"
dirname = "/home/eugene/results/immature/clusters/11/"

fileActiveSynapses = os.path.join(dirname, "RA_RA_active_connections_" + str(trial_number) + ".bin")
fileSuperSynapses = os.path.join(dirname, "RA_RA_super_connections_" + str(trial_number) + ".bin")
fileMature = os.path.join(dirname, "mature_" + str(trial_number) + ".bin")
fileWeights = os.path.join(dirname, "weights_" + str(trial_number) + ".bin")
fileTraining = os.path.join(dirname, "training_neurons.bin")

(_, _, mature_indicators) = reading.read_mature_indicators(fileMature)
(_, _, active_synapses) = reading.read_synapses(fileActiveSynapses)
(_, _, super_synapses) = reading.read_synapses(fileSuperSynapses)
(_, _, weights) = reading.read_weights(fileWeights) 
training_neurons = reading.read_training_neurons(fileTraining)


N_RA = len(active_synapses)
###############################################################
### Calculate number of active and super in and out degrees ###
###############################################################
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

print "Active in degree for neuron {0} : {1}".format(5, num_active_inputs[5])
print "Super in degree for neuron {0} : {1}".format(5, num_super_inputs[5])

################################
### Find neurons in the chain ##
################################

chain_neurons = set()

new_added = set(training_neurons)

while len(new_added) > 0:
    source_id = next(iter(new_added))
    for target in super_synapses[source_id]:
        if target not in chain_neurons:
            new_added.add(target)
        
    new_added.remove(source_id)
    chain_neurons.add(source_id)
    
print chain_neurons

####################################################
### Find supersynapses from neurons not in chain ###
####################################################
for i in range(len(super_synapses)):
    if i not in chain_neurons and len(super_synapses[i]) > 0:
        print "Supersynapses from neuron {0} : {1}".format(i, super_synapses[i])
        
        
        
        
targets = [5]

for i in range(len(active_synapses)):
    if i not in chain_neurons:
        for target in targets:
            if target in active_synapses[i]:
                print "Active synapse from neuron {0} -> {1} : {2}".format(i, target, weights[i][target])


print mature_indicators[5]
print mature_indicators[993]

print weights[993][5]

for i in range(len(super_synapses)):
    if i not in training_neurons:
        for target in targets:
            if target in active_synapses[i]:
                print "Active synapse {0} -> {1}; weight = {2}".format(i, target, weights[i][target])
    
        
    
