# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 19:07:46 2017

@author: jingroup

Script analyzes differences between grown chains
"""
import indirect_connections as ic
import reading
import numpy as np

SUPERSYNAPTIC_THRESHOLD = 0.025
ACTIVATION_THRESHOLD = 0.0005

new_id = 2800
old_id = 2500

dataDir = "/home/eugene/Output/networks/gabaMaturation100217/"

file_weights_new = dataDir + "weights_" + str(new_id) + "_.bin"
file_weights_old = dataDir + "weights_" + str(old_id) + "_.bin"

file_super_new = dataDir + "RA_RA_super_connections_" + str(new_id) + "_.bin"
file_super_old = dataDir + "RA_RA_super_connections_" + str(old_id) + "_.bin"

file_active_new = dataDir + "RA_RA_active_connections_" + str(new_id) + "_.bin"
file_active_old = dataDir + "RA_RA_active_connections_" + str(old_id) + "_.bin"

file_mature_new = dataDir + "mature_" + str(new_id) + "_.bin"
file_mature_old = dataDir + "mature_" + str(old_id) + "_.bin"


def get_synapses(output_connections_id, output_connections_weights):
    """
    Extracts synapses and connected neurons from array with
    output connections
    """
    synapses = []
    synaptic_weights = []
    connected = set()

    for i, targets in enumerate(output_connections_id):
        if len(targets) > 0:
            connected.add(i)
        for j, target in enumerate(targets):
            synapses.append((i, target))
            synaptic_weights.append(output_connections_weights[i][j])
            connected.add(target)
    
    return (synapses, synaptic_weights, connected)


trial_number, gaba_potential_new, firing_rate_new, remodeled_new, mature_new = reading.read_maturation_info(file_mature_new)
trial_number, gaba_potential_old, firing_rate_old, remodeled_old, mature_old = reading.read_maturation_info(file_mature_old)

(N_RA, trial_number, weights_new) = reading.read_weights(file_weights_new) 
(N_RA, trial_number, weights_old) = reading.read_weights(file_weights_old) 

activily_connected = set()
       
# get supersynapses              
(N_RA_old, output_supersynapses_id_old, output_supersynapses_weights_old) = reading.read_connections(file_super_old)
(N_RA_new, output_supersynapses_id_new, output_supersynapses_weights_new) = reading.read_connections(file_super_new)

supersynapses_id_old, _, strongly_connected_old = get_synapses(output_supersynapses_id_old, output_supersynapses_weights_old)
supersynapses_id_new, _, strongly_connected_new = get_synapses(output_supersynapses_id_new, output_supersynapses_weights_new)

# get active synapses
(N_RA_old, output_active_synapses_id_old, output_active_synapses_weights_old) = reading.read_connections(file_active_old)
(N_RA_new, output_active_synapses_id_new, output_active_synapses_weights_new) = reading.read_connections(file_active_new)

active_synapses_id_old, active_synapses_weights_old, actively_connected_old = get_synapses(output_active_synapses_id_old, output_active_synapses_weights_old)
active_synapses_id_new, active_synapses_weights_new, actively_connected_new = get_synapses(output_active_synapses_id_new, output_active_synapses_weights_new)

# print all gained and lost supersynapses
new_supersynapses = set(supersynapses_id_new) - set(supersynapses_id_old)
lost_supersynapses = set(supersynapses_id_old) - set(supersynapses_id_new)

print "new supersynapses: ",new_supersynapses
print "lost supersynapses: ",lost_supersynapses

# print all gained and lost active synapses
new_active_synapses = set(active_synapses_id_new) - set(active_synapses_id_old)
lost_active_synapses = set(active_synapses_id_old) - set(active_synapses_id_new)

print "# new active synapses: ",len(new_active_synapses)
print "# lost active synapses: ",len(lost_active_synapses)


strongly_connected_difference = set(strongly_connected_new).symmetric_difference(strongly_connected_old)

print "difference between strongly connected neurons: ",strongly_connected_difference

neuron_id = 391

print gaba_potential_old[neuron_id]
print gaba_potential_new[neuron_id]
print firing_rate_old[neuron_id]
print firing_rate_new[neuron_id]

# get all input active synapses to neuron with id neuron_id

input_id = []
input_weight = []

for i, (source, target) in enumerate(active_synapses_id_new):
    #print (source, target)
    if target == neuron_id:
        input_id.append(source)
        input_weight.append(active_synapses_weights_new[i])
        
print "Active new input synapses to neuron {0}:".format(neuron_id), input_id
print "Active new input weights to neuron {0}:".format(neuron_id), input_weight

# check indirect connections from active input neurons to neuron with id neuron_id
RA2I = dataDir + "RA_I_connections.bin"
I2RA = dataDir + "I_RA_connections.bin"
   
(N_RA, RA2I_targets, RA2I_targets_G) = reading.read_connections(RA2I)
(N_I, I2RA_targets, I2RA_targets_G) = reading.read_connections(I2RA)

print "indirect connections from active new inputs: ", ic.get_indirect_connections_to_target(input_id, [neuron_id], RA2I_targets, I2RA_targets)

# check indirect connections neurons that formed supersynapses to neuron with id neuron_id

supersynapse_source = []

for (source, target) in new_supersynapses:
    if target == neuron_id:
        supersynapse_source.append(source)

print "indirect connections from neurons that formed supersynapses: ", ic.get_indirect_connections_to_target(supersynapse_source, [neuron_id], RA2I_targets, I2RA_targets)