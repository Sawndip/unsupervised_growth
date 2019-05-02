# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 11:07:31 2017

@author: jingroup

Script checks the outcome of replacement
"""

import reading
import matplotlib.pyplot as plt
import numpy as np
import space

replacement_trial = 2499 # trial number when replacement occured
reference_trial = 2400 # reference trial

neuron_id = 53 # id of some replaced neuron

GABA_IMMATURE = -56.0 # mature value of GABA reversal potential in mV
BURST_TIME_REPLACED = -200.0 # last burst time of replaced neuron

arrangement = "sphere"
dim = space.num_coordinates[arrangement] # dimensionality
dirname = "/home/eugene/Output/networks/sphere_170717_hodgkin/"

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

filename_replacement_history_before = dirname + "replacement_history_before_replacement_trial_" + str(replacement_trial) + "_.bin"
filename_replacement_history_after = dirname + "replacement_history_after_replacement_trial_" + str(replacement_trial) + "_.bin"

filename_replaced = dirname + "replaced_neurons.bin"

replacement_time, replaced_neurons = reading.read_replaced_neurons(filename_replaced)

ind = replacement_time.index(replacement_trial)

print "Replaced neurons: ",replaced_neurons[ind]
print "Replacement time: ", replacement_time[ind]



def change_in_fixed_connections(targets_before, targets_after, replaced):
    neurons_with_changed_targets = np.where(np.array(targets_before) != np.array(targets_after))
    
    not_replaced_with_changed_targets = np.setdiff1d(neurons_with_changed_targets, replaced, assume_unique=False)
    
    if len(not_replaced_with_changed_targets) > 0:
        print not_replaced_with_changed_targets
    else:
        print "Only replaced neurons changed their interneuron targets"
    #print np.array(targets_before)[np.where(true_array)]
    #print np.array(targets_after)[np.where(true_array)]

#########################################
# Check HVC(RA) coordinates replacement
#########################################

coord_before = reading.read_coordinates(dim, filename_RAxy_before)
coord_after = reading.read_coordinates(dim, filename_RAxy_after)

coord_equal = coord_before == coord_after

N_RA = coord_before.shape[0]

not_replaced = [] # neurons not replaced

for i in range(N_RA):
    if i not in replaced_neurons[6]:
        not_replaced.append(i)

print "Neurons not replaced: ",not_replaced

replacement_correct = True

for i in range(N_RA):
    if i in not_replaced:
        if np.logical_and.reduce(coord_equal[i]) == False:
            replacement_correct = False
            print "Coordinates changed for neuron {0} that was not replaced".format(i)
            break
    else:
        if np.logical_and.reduce(coord_equal[i]) == True:
            replacement_correct = False
            print "Coordinates didn't change for neuron {0} that was replaced".format(i)
            break
        
print "Coordinate replacement went correctly: ",replacement_correct

#########################################
# Check HVC(RA) -> HVC(I) connections replacement
#########################################
(N_RA, targets_RAI_ID_before, targets_RAI_G_before) = reading.read_connections(filename_RAI_before)
(_, targets_RAI_ID_after, targets_RAI_G_after) = reading.read_connections(filename_RAI_after)

replacement_correct = True

for i in range(N_RA):
    if i in not_replaced:
        if targets_RAI_ID_before[i] != targets_RAI_ID_after[i]:
            replacement_correct = False
            print "Connections to HVC(I) neurons changed for neuron {0} that was not replaced".format(i)
            break
    else:
         if targets_RAI_ID_before[i] == targets_RAI_ID_after[i]:
            replacement_correct = False
            print "Connections to HVC(I) neurons didn't change for neuron {0} that was replaced".format(i)
            break

print "HVC(RA) -> HVC(I) connections replacement went correctly: ",replacement_correct

#########################################
# Check HVC(I) -> HVC(RA) connections replacement
#########################################
(N_I, targets_IRA_ID_before, targets_IRA_G_before) = reading.read_connections(filename_IRA_before)
(_, targets_IRA_ID_after, targets_IRA_G_after) = reading.read_connections(filename_IRA_after)

replacement_correct = True

for i in range(N_I):
    # get difference in targets for each interneuron and see if this difference
    # is only in the connections to replaced neurons
    target_difference = set(targets_IRA_ID_before[i]).symmetric_difference(set(targets_IRA_ID_after[i]))
    
    for j in target_difference:
        if j in not_replaced:
            replacement_correct = False
            print "Connections to not replaced HVC(RA) neurons changed for HVC(I) neuron {0}".format(i)
            break
            

print "HVC(I) -> HVC(RA) connections replacement went correctly: ",replacement_correct

#########################################
# Check active connections replacement
#########################################
(_, targets_ID_before, targets_G_before) = reading.read_connections(filename_active_before)
(_, targets_ID_after, targets_G_after) = reading.read_connections(filename_active_after)

replacement_correct = True

for i in range(N_RA):
    if i in not_replaced:
        # get difference in targets for each HVC(RA) neuron that was not replaced
        # and see if this difference is only in the connections to replaced neurons
        target_difference = set(targets_ID_before[i]).symmetric_difference(set(targets_ID_after[i]))
        
        for j in target_difference:
            if j in not_replaced:
                replacement_correct = False
                print "Active connections from not replaced HVC(RA) neurons to not replaced HVC(RA) neuron changed for neuron {0}".format(i)
                break
    else:
        # for replaced neurons just check if target list is empty
        if len(targets_ID_after[i]) > 0:
            replacement_correct = False
            print "Active connections from replaced HVC(RA) neurons exist for neuron {0}".format(i)
            break

print "HVC(RA) -> HVC(RA) active connections replacement went correctly: ",replacement_correct

#########################################
# Check super connections replacement
#########################################
(_, targets_ID_before, targets_G_before) = reading.read_connections(filename_super_before)
(_, targets_ID_after, targets_G_after) = reading.read_connections(filename_super_after)

replacement_correct = True

for i in range(N_RA):
    if i in not_replaced:
        # get difference in targets for each HVC(RA) neuron that was not replaced
        # and see if this difference is only in the connections to replaced neurons
        target_difference = set(targets_ID_before[i]).symmetric_difference(set(targets_ID_after[i]))
        
        for j in target_difference:
            if j in not_replaced:
                replacement_correct = False
                print "Super connections from not replaced HVC(RA) neurons to not replaced HVC(RA) neuron changed for neuron {0}".format(i)
                break
    else:
        # for replaced neurons just check if target list is empty
        if len(targets_ID_after[i]) > 0:
            replacement_correct = False
            print "Super connections from replaced HVC(RA) neurons exist for neuron {0}".format(i)
            break

print "HVC(RA) -> HVC(RA) super connections replacement went correctly: ",replacement_correct

#########################################
# Check weight replacement
#########################################
_, _, weights_before = reading.read_weights(filename_weights_before)
_, _, weights_after = reading.read_weights(filename_weights_after)


ind_matrix_with_changes = np.where(np.array(weights_before) != np.array(weights_after))

replacement_correct = True
        
for i in range(ind_matrix_with_changes[0].shape[0]):
    if ind_matrix_with_changes[0][i] in not_replaced:
        # check if weight changed from not replaced to not replaced neuron
        if ind_matrix_with_changes[1][i] in not_replaced:
            replacement_correct = False
            print "Weight from not replaced HVC(RA) neuron to not HVC(RA) neuron changed for neuron {0}".format(ind_matrix_with_changes[0][i])
            break
        
        
    
for i in range(N_RA):
    # check if all weights from replaced neurons are set to zero
    if i not in not_replaced:
        if np.count_nonzero(weights_after[i]) > 0:
            replacement_correct = False
            print "Weight from replaced HVC(RA) neuron are not zero for neuron {0}".format(i)
            break
    # check if weight didn't change from not replaced to replaced neuron
    else:
        for j in range(N_RA):
            if j not in not_replaced:
                if not np.isclose(weights_before[i][j], 0.0) and np.isclose(weights_before[i][j], weights_after[i][j]):
                    replacement_correct = False
                    print "Weight from not replaced HVC(RA) neuron to replaced HVC(RA) neuron didn't change for neuron {0}".format(i)
                    break
        
print "Weights replacement went correctly: ",replacement_correct

#########################################
# Check maturation replacement
#########################################

# check maturation parameters

_, gaba_potential_before, firing_rate_short_before, firing_rate_long_before, remodeled_before = reading.read_maturation_info(filename_maturation_before)
_, gaba_potential_after, firing_rate_short_after, firing_rate_long_after, remodeled_after = reading.read_maturation_info(filename_maturation_after)

replacement_correct = True

for i in range(N_RA):
    if i in not_replaced:
        # check if GABA potential remained the same for no replaced neurons
        if not np.isclose(gaba_potential_before[i], gaba_potential_after[i]):
            replacement_correct = False
            print "GABA potential for not replaced HVC(RA) neuron changed for neuron {0}".format(i)
            break
        
        # check if firing rate short remained the same for no replaced neurons
        if not np.isclose(firing_rate_short_before[i], firing_rate_short_after[i]):
            replacement_correct = False
            print "Firing rate short for not replaced HVC(RA) neuron changed for neuron {0}".format(i)
            break
        
        # check if firing rate long remained the same for no replaced neurons
        if not np.isclose(firing_rate_long_before[i], firing_rate_long_after[i]):
            replacement_correct = False
            print "Firing rate long for not replaced HVC(RA) neuron changed for neuron {0}".format(i)
            break
        
        # check if remodeled indicators remained the same for no replaced neurons
        if remodeled_after[i] != remodeled_before[i]:
            replacement_correct = False
            print "Remodeled for not replaced HVC(RA) neuron changed for neuron {0}".format(i)
            break
    else:
        # check if gaba potential was set to gaba immature for replaced neurons
        if not np.isclose(gaba_potential_after[i], GABA_IMMATURE):
            replacement_correct = False
            print "GABA potential of replaced HVC(RA) neuron didn't change to GABA_IMMATURE {0} for neuron {1}".format(GABA_IMMATURE, i)
            break
    
        # check if firing_rate_short was set to zero for replaced neurons
        if not np.isclose(firing_rate_short_after[i], 0.0):
            replacement_correct = False
            print "Firing rate short of replaced HVC(RA) neuron didn't change to 0 for neuron {0}".format(i)
            break
        
         # check if firing_rate_long was set to zero for replaced neurons
        if not np.isclose(firing_rate_long_after[i], 0.0):
            replacement_correct = False
            print "Firing rate long of replaced HVC(RA) neuron didn't change to 0 for neuron {0}".format(i)
            break
        
         # check if remodeled indicator was set to zero for replaced neurons
        if remodeled_after[i] != 0:
            replacement_correct = False
            print "Remodeled state of replaced HVC(RA) neuron didn't change to 0 for neuron {0}".format(i)
            break
        
print "Maturation replacement went correctly: ",replacement_correct


#########################################
# Check last dendritic burst times replacement
#########################################

last_dend_spike_times_before = reading.read_last_dend_spike_times(filename_last_dend_spike_times_before)
last_dend_spike_times_after = reading.read_last_dend_spike_times(filename_last_dend_spike_times_after)

replacement_correct = True

for i in range(N_RA):
    if i in not_replaced:
        # check if last dendritic burst time changed for not replaced neurons
        if not np.isclose(last_dend_spike_times_before[i], last_dend_spike_times_after[i]):
            replacement_correct = False
            print "Last dendritic burst time for not replaced HVC(RA) neuron changed for neuron {0}".format(i)
            break
       
    else:
        # check if last dendritic burst time was set to BURST_TIME_REPLACED for replaced neurons
        if not np.isclose(last_dend_spike_times_after[i], BURST_TIME_REPLACED):
            replacement_correct = False
            print "Last dendritic burst time of replaced HVC(RA) neuron didn't change to BURST_TIME_REPLACED {0} for neuron {1}".format(BURST_TIME_REPLACED, i)
            break

print "Last dendritic burst time replacement went correctly: ",replacement_correct

#########################################
# Check replacement history replacement
#########################################
time_from_previous_replacement_before = reading.read_replacement_history(filename_replacement_history_before)
time_from_previous_replacement_after = reading.read_replacement_history(filename_replacement_history_after)

replacement_correct = True

for i in range(N_RA):
    if i in not_replaced:
        # check if replacement history remained the same for not replaced neurons
        if time_from_previous_replacement_before[i] != time_from_previous_replacement_after[i]:
            replacement_correct = False
            print "Time from previous replacement for not replaced HVC(RA) neuron changed for neuron {0}".format(i)
            break
       
    else:
        # check if replacement history was set to zero for replaced neurons
        if time_from_previous_replacement_after[i] != 0:
            replacement_correct = False
            print "Time from previous replacement of replaced HVC(RA) neuron didn't change to zero for neuron {0}".format(i)
            break

print "Times from previous replacement replacement went correctly: ",replacement_correct