# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 19:36:16 2017

@author: jingroup

Script reads data from files
"""

import reading

filename_before_RA = "/home/eugene/Output/networks/gabaMaturation100217/RA_I_connections.bin"

(N_RA_before, targets_ID_before_RA, targets_G_before_RA) = reading.read_connections(filename_before_RA)

print N_RA_before
#print targets_ID_before


filename_after_RA = "/home/eugene/Output/networks/gabaMaturation100217/RA_I_connections_after_addition.bin"

(N_RA_after, targets_ID_after_RA, targets_G_after_RA) = reading.read_connections(filename_after_RA)

print N_RA_after

print targets_ID_before_RA == targets_ID_after_RA[:N_RA_before]
#print targets_ID_after[N_RA_before:]

#print targets_ID_before

filename_before_I = "/home/eugene/Output/networks/gabaMaturation100217/I_RA_connections.bin"

(N_I_before, targets_ID_before_I, targets_G_before_I) = reading.read_connections(filename_before_I)

print N_I_before

filename_after_I = "/home/eugene/Output/networks/gabaMaturation100217/I_RA_connections_after_addition.bin"

(N_I_after, targets_ID_after_I, targets_G_after_I) = reading.read_connections(filename_after_I)

print N_I_after

#print targets_ID_before
I_neuron_id = 30

print targets_ID_before_I[I_neuron_id]
print targets_ID_after_I[I_neuron_id]

filename_weights_before = "/home/eugene/Output/networks/gabaMaturation100217/weights_before_neuron_addition.bin"
filename_weights_after = "/home/eugene/Output/networks/gabaMaturation100217/weights_after_neuron_addition.bin"


(N_RA_before, trial_number, weights_before) = reading.read_weights(filename_weights_before)
(N_RA_after, trial_number, weights_after) = reading.read_weights(filename_weights_after)

print N_RA_before
print N_RA_after

neuron_id = 122
target = 22

equality = np.zeros((N_RA_before, 1), dtype=bool)

for i in range(N_RA_before):
    equality[i] = np.logical_and.reduce(weights_before[i] == weights_after[i][:N_RA_before])

weights_after_smaller = weights_after[0:N_RA_before, 0:N_RA_before]
print np.logical_and.reduce(weights_before == weights_after[0:N_RA_before, 0:N_RA_before])
print weights_before == weights_after[0:N_RA_before, 0:N_RA_before]

print np.array_equal(weights_before, weights_after[0:N_RA_before, 0:N_RA_before])
#print equality

#print weights_before[18]
#print weights_after[18][:N_RA_before]

print weights_before[neuron_id][target]
print weights_after[neuron_id][target]


filename_id_array_before = "/home/eugene/Output/networks/gabaMaturation100217/global_index_array.bin"
filename_id_array_after = "/home/eugene/Output/networks/gabaMaturation100217/global_index_array_after_neuron_addition.bin"


(N_RA_before, Id_RA_global_before) = reading.read_global_index_array(filename_id_array_before)
(N_RA_after, Id_RA_global_after) = reading.read_global_index_array(filename_id_array_after)

print N_RA_before
print N_RA_after
print Id_RA_global_before
print Id_RA_global_after


filename_mature_before = "/home/eugene/Output/networks/gabaMaturation100217/mature_before_neuron_addition.bin"
filename_mature_after = "/home/eugene/Output/networks/gabaMaturation100217/mature_after_neuron_addition.bin"

trial_number, gaba_potential_before, firing_rate_before, remodeled_before, mature_before = reading.read_maturation_info(filename_mature_before)
trial_number, gaba_potential_after, firing_rate_after, remodeled_after, mature_after = reading.read_maturation_info(filename_mature_after)

print gaba_potential_before == gaba_potential_after[:N_RA_before]

print firing_rate_before == firing_rate_after[:N_RA_before]
print remodeled_before == remodeled_after[:N_RA_before]
print mature_before == mature_after[:N_RA_before]

print gaba_potential_after
print firing_rate_after