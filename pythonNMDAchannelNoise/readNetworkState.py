# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 15:27:05 2017

@author: jingroup

Script reads the network state
"""

import reading
import indirect_connections as connect

trial = 7000 # starting point

dirname = "/mnt/hodgkin_home/eugene/Output/networks/gabaMaturation090617_huxley/"

filename_RAxy = dirname + "RA_xy_" + str(trial) + "_.bin"
filename_RAI = dirname + "RA_I_connections_" + str(trial) + "_.bin"
filename_IRA = dirname + "I_RA_connections_" + str(trial) + "_.bin"
filename_active = dirname + "RA_RA_active_connections_" + str(trial) + "_.bin"
filename_super = dirname + "RA_RA_super_connections_" + str(trial) + "_.bin"
filename_maturation = dirname + "mature_" + str(trial) + "_.bin"
filename_weights = dirname + "weights_" + str(trial) + "_.bin"
filename_last_dend_spike_times = dirname + "last_dendritic_spike_times_before_replacement_trial_" + str(trial) + "_.bin"


######################
# HVC(RA) <-> HVC(I) connections #
######################

(N_RA, targets_ID_RAI, targets_G_RAI) = reading.read_connections(filename_RAI)
(N_I, targets_ID_IRA, targets_G_IRA) = reading.read_connections(filename_IRA)

######################
# HVC(RA) <-> HVC(RA) connections #
######################

(N_RA, trial_number, weights) = reading.read_weights(filename_weights)
######################
# Read supersynapses #
######################

(N_RA, targets_ID_super, targets_G_super) = reading.read_connections(filename_super)

#print N_RA
#print targets_ID_super

print "\nSupersynapses:\n"

strongly_connected = set()

for i in range(len(targets_ID_super)):
    if len(targets_ID_super[i]) > 0:    
        strongly_connected.add(i)
    for j in range(len(targets_ID_super[i])):
        print "Supersynapse {0} -> {1} Weight = {2}".format(i, targets_ID_super[i][j], targets_G_super[i][j])
        strongly_connected.add(targets_ID_super[i][j])
        

print "Strongly connected neurons: ",strongly_connected

##
# print some weights

print weights[0][117]
print weights[9][117]
print weights[121][117]
print weights[277][117]

print weights[162][117]


######################
# Read maturation info #
######################

trial_number, gaba_potential, firing_rate_short, firing_rate_long, remodeled = reading.read_maturation_info(filename_maturation)

print "\nMaturation state of strongly connected neurons:\n"

for i in strongly_connected:
    print "Neuron {0}; gaba potential = {1}; firing_rate_short = {2}; firing_rate_long = {3}; remodeled = {4}".format(i, gaba_potential[i], firing_rate_short[i], firing_rate_long[i], remodeled[i])
    

source = [0, 9, 121, 277]
targets = [222]

connections_from_source, num_connections_from_source, \
num_convergent_inputs_to_interneurons_connected_to_targets = connect.find_indirect_connections_to_target(source, targets, targets_ID_RAI, targets_ID_IRA)

print connections_from_source
print num_connections_from_source
print num_convergent_inputs_to_interneurons_connected_to_targets