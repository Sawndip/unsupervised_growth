# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 11:07:31 2017

@author: jingroup

Script checks the outcome of replacement
"""

import reading
import matplotlib.pyplot as plt
import numpy as np

replacement_trial = 500 # trial number when replacement occured
reference_trial = 500 # reference trial

neuron_id = 53 # id of some replaced neuron

dirname = "/mnt/hodgkin_home/eugene/Output/networks/networkTest/"

filename_RAxy_before = dirname + "RA_xy_" + str(reference_trial) + "_.bin"
filename_RAxy_after = dirname + "RA_xy_after_replacement_trial_" + str(replacement_trial) + "_.bin"

filename_RAI_before = dirname + "RA_I_connections_" + str(reference_trial) + "_.bin"
filename_RAI_after = dirname + "RA_I_connections_after_replacement_trial_" + str(replacement_trial) + "_.bin"

filename_IRA_before = dirname + "I_RA_connections_" + str(reference_trial) + "_.bin"
filename_IRA_after = dirname + "I_RA_connections_after_replacement_trial_" + str(replacement_trial) + "_.bin"

filename_active_before = dirname + "RA_RA_active_connections_before_replacement_trial_" + str(replacement_trial) + "_.bin"
filename_active_after = dirname + "RA_RA_active_connections_after_replacement_trial_" + str(replacement_trial) + "_.bin"

filename_super_before = dirname + "RA_RA_super_connections_before_replacement_trial_" + str(replacement_trial) + "_.bin"
filename_super_after = dirname + "RA_RA_super_connections_after_replacement_trial_" + str(replacement_trial) + "_.bin"

filename_maturation_before = dirname + "mature_before_replacement_trial_" + str(replacement_trial) + "_.bin"
filename_maturation_after = dirname + "mature_after_replacement_trial_" + str(replacement_trial) + "_.bin"

filename_weights_before = dirname + "weights_before_replacement_trial_" + str(replacement_trial) + "_.bin"
filename_weights_after = dirname + "weights_after_replacement_trial_" + str(replacement_trial) + "_.bin"

filename_last_dend_spike_times_before = dirname + "last_dendritic_spike_times_before_replacement_trial_" + str(replacement_trial) + "_.bin"
filename_last_dend_spike_times_after = dirname + "last_dendritic_spike_times_after_replacement_trial_" + str(replacement_trial) + "_.bin"


filename_replaced = dirname + "replaced_neurons.bin"

replacement_time, replaced_neurons = reading.read_replaced_neurons(filename_replaced)


print "Replaced neurons: ",replaced_neurons
print "Replacement time: ", replacement_time

def change_in_fixed_connections(targets_before, targets_after, replaced):
    neurons_with_changed_targets = np.where(np.array(targets_before) != np.array(targets_after))
    
    not_replaced_with_changed_targets = np.setdiff1d(neurons_with_changed_targets, replaced, assume_unique=False)
    
    if len(not_replaced_with_changed_targets) > 0:
        print not_replaced_with_changed_targets
    else:
        print "Only replaced neurons changed their interneuron targets"
    #print np.array(targets_before)[np.where(true_array)]
    #print np.array(targets_after)[np.where(true_array)]

# check RA coordinates
(xx_before, yy_before) = reading.read_coordinates(filename_RAxy_before)
(xx_after, yy_after) = reading.read_coordinates(filename_RAxy_after)


print "Coordinates before: ", (xx_before[neuron_id], yy_before[neuron_id])
print "Coordinates after: ", (xx_after[neuron_id], yy_after[neuron_id])


# check RA -> I and I -> RA connections
(N_RA, targets_RAI_ID_before, targets_RAI_G_before) = reading.read_connections(filename_RAI_before)
(_, targets_RAI_ID_after, targets_RAI_G_after) = reading.read_connections(filename_RAI_after)

(_, targets_IRA_ID_before, targets_IRA_G_before) = reading.read_connections(filename_IRA_before)
(_, targets_IRA_ID_after, targets_IRA_G_after) = reading.read_connections(filename_IRA_after)

change_in_fixed_connections(targets_RAI_ID_before, targets_RAI_ID_after, replaced_neurons[0])

print "inhibitory targets before:", targets_RAI_ID_before[neuron_id]
print "inhibitory targets after:", targets_RAI_ID_after[neuron_id]

# check active connections
(_, targets_ID_before, targets_G_before) = reading.read_connections(filename_active_before)
(_, targets_ID_after, targets_G_after) = reading.read_connections(filename_active_after)



print "active synapses before: ",targets_ID_before[neuron_id]
print "active synapses after: ",targets_ID_after[neuron_id]

# check super connections
(_, targets_ID_before, targets_G_before) = reading.read_connections(filename_super_before)
(_, targets_ID_after, targets_G_after) = reading.read_connections(filename_super_after)


print "super synapses before: ",targets_ID_before[neuron_id]
print "super synapses after: ",targets_ID_after[neuron_id]


# check weights

_, _, weights_before = reading.read_weights(filename_weights_before)
_, _, weights_after = reading.read_weights(filename_weights_after)

print "weights before: ",weights_before[neuron_id]
print "weights after: ",weights_after[neuron_id]

source = 0

print [i for i, w in enumerate(weights_before[source]) if w > 0]

print weights_before[source][neuron_id]
print weights_after[source][neuron_id]


print weights_before[neuron_id] == weights_after[neuron_id]


# check maturation parameters

_, gaba_potential_before, firing_rate_short_before, firing_rate_long_before, remodeled_before = reading.read_maturation_info(filename_maturation_before)
_, gaba_potential_after, firing_rate_short_after, firing_rate_long_after, remodeled_after = reading.read_maturation_info(filename_maturation_after)

print "gaba potential before: ",gaba_potential_before[neuron_id]
print "gaba potential after: ",gaba_potential_after[neuron_id]


print "remodeled before: ",remodeled_before[neuron_id]
print "remodeled after: ",remodeled_after[neuron_id]


print "short firing rate before: ",firing_rate_short_before[neuron_id]
print "short firing rate after: ",firing_rate_short_after[neuron_id]

print "long firing rate before: ",firing_rate_long_before[neuron_id]
print "long firing rate after: ",firing_rate_long_after[neuron_id]


# check last dendritic spikes

last_dend_spike_times_before = reading.read_last_dend_spike_times(filename_last_dend_spike_times_before)
last_dend_spike_times_after = reading.read_last_dend_spike_times(filename_last_dend_spike_times_after)

print last_dend_spike_times_before == last_dend_spike_times_after

print "last burst times before :", last_dend_spike_times_before
print "last burst times after :", last_dend_spike_times_after

#print "last burst times new :", last_dend_spike_times_new
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
