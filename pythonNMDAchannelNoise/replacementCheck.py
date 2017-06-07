# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 11:07:31 2017

@author: jingroup

Script checks the outcome of replacement
"""

import reading
import matplotlib.pyplot as plt
import numpy as np

num_file = 9 # file number
neuron_id = 7 # id of some replaced neuron



filename_RAI_before = "/home/eugene/Output/networks/networkTest/RA_I_connections.bin"
filename_RAI_after = "/home/eugene/Output/networks/networkTest/RA_I_connections_after_replacement_trial_" + str(num_file) + "_.bin"

filename_IRA_before = "/home/eugene/Output/networks/networkTest/I_RA_connections.bin"
filename_IRA_after = "/home/eugene/Output/networks/networkTest/I_RA_connections_after_replacement_trial_" + str(num_file) + "_.bin"

filename_active_before = "/home/eugene/Output/networks/networkTest/RA_RA_active_connections_before_replacement_trial_" + str(num_file) + "_.bin"
filename_active_after = "/home/eugene/Output/networks/networkTest/RA_RA_active_connections_after_replacement_trial_" + str(num_file) + "_.bin"

filename_maturation_before = "/home/eugene/Output/networks/networkTest/mature_before_replacement_trial_" + str(num_file) + "_.bin"
filename_maturation_after = "/home/eugene/Output/networks/networkTest/mature_after_replacement_trial_" + str(num_file) + "_.bin"

filename_weights_before = "/home/eugene/Output/networks/networkTest/weights_before_replacement_trial_" + str(num_file) + "_.bin"
filename_weights_after = "/home/eugene/Output/networks/networkTest/weights_after_replacement_trial_" + str(num_file) + "_.bin"


filename_replaced = "/home/eugene/Output/networks/networkTest/replaced_neurons.bin"

replacement_time, replaced_neurons = reading.read_replaced_neurons(filename_replaced)

print replacement_time
print "Replaced neurons: ",replaced_neurons[num_file-1]
print "Replacement time: ", replacement_time[num_file-1]

def change_in_fixed_connections(targets_before, targets_after, replaced):
    neurons_with_changed_targets = np.where(np.array(targets_before) != np.array(targets_after))
    
    not_replaced_with_changed_targets = np.setdiff1d(neurons_with_changed_targets, replaced, assume_unique=False)
    
    if len(not_replaced_with_changed_targets) > 0:
        print not_replaced_with_changed_targets
    else:
        print "Only replaced neurons changed their interneuron targets"
    #print np.array(targets_before)[np.where(true_array)]
    #print np.array(targets_after)[np.where(true_array)]

# check RA -> I and I -> RA connections
(N_RA, targets_RAI_ID_before, targets_RAI_G_before) = reading.read_connections(filename_RAI_before)
(_, targets_RAI_ID_after, targets_RAI_G_after) = reading.read_connections(filename_RAI_after)

(_, targets_IRA_ID_before, targets_IRA_G_before) = reading.read_connections(filename_IRA_before)
(_, targets_IRA_ID_after, targets_IRA_G_after) = reading.read_connections(filename_IRA_after)

change_in_fixed_connections(targets_RAI_ID_before, targets_RAI_ID_after, replaced_neurons)

print targets_RAI_ID_before[neuron_id]
print targets_RAI_ID_after[neuron_id]


(_, targets_ID_before, targets_G_before) = reading.read_connections(filename_active_before)
(_, targets_ID_after, targets_G_after) = reading.read_connections(filename_active_after)



print targets_ID_before[neuron_id]
print targets_ID_after[neuron_id]

_, _, weights_before = reading.read_weights(filename_weights_before)
_, _, weights_after = reading.read_weights(filename_weights_after)

#print weights_before[neuron_id]
#print weights_after[neuron_id]

source = 0

print [i for i, w in enumerate(weights_before[source]) if w > 0]

print weights_before[source][neuron_id]
print weights_after[source][neuron_id]


print weights_before[neuron_id] == weights_after[neuron_id]




_, gaba_potential_before, firing_rate_short_before, firing_rate_long_before, remodeled_before = reading.read_maturation_info(filename_maturation_before)
_, gaba_potential_after, firing_rate_short_after, firing_rate_long_after, remodeled_after = reading.read_maturation_info(filename_maturation_after)



print firing_rate_short_before[neuron_id]
print firing_rate_short_after[neuron_id]

print firing_rate_long_before[neuron_id]
print firing_rate_long_after[neuron_id]
#==============================================================================
# filenameStates = "/home/eugene/Output/networks/gabaMaturation010317/maturation_time_sequence.bin"
# 
# (target, t, remodeled, gaba_potential, firing_rate) = reading.read_maturation_time_sequence(filenameStates)
# 
# 
# # plot firing rate and gaba potential
# xmax = 7000 # max trial to show
# 
# f = plt.figure()
# 
# ax1 = f.add_subplot(211)
# 
# ax1.set_title("Firing rate of neuron {0}".format(neuron_id))
# ax1.plot(t, firing_rate[neuron_id])
# ax1.set_ylabel("r")
# ax1.set_xlim([0, xmax])
# 
# ax2 = f.add_subplot(212)
# 
# ax2.set_title("GABA reverse potential of {0}".format(neuron_id))
# ax2.plot(t, gaba_potential[neuron_id])
# ax2.set_xlabel("t (# trial)")
# ax2.set_ylabel("$E_{GABA}$ mV")
# ax2.set_xlim([0, xmax])
# 
# 
# plt.show()   
#==============================================================================
