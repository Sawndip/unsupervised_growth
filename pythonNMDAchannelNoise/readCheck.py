# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 16:44:37 2017

@author: jingroup

Script checks that data was read correctly after chain growth was restarted
"""
import reading
import numpy as np
import space

directory = "/home/eugene/Output/networks/sphere_170717_hodgkin/"
start_number = 22400 # starting point
arrangement = "sphere"
dim = space.num_coordinates[arrangement] # dimensionality 

###################################
# check coordinates of neurons
###################################

filename_old = directory + "RA_xy_" + str(start_number) + "_.bin"
filename_new = directory + "RA_xy_" + str(start_number) + "_NEW.bin"

coord_old = reading.read_coordinates(dim, filename_old)
coord_new = reading.read_coordinates(dim, filename_new)

print "Old and new coordinates match: ",np.logical_and.reduce(np.logical_and.reduce(coord_old == coord_new))



###################################
# check RA -> I connections
###################################

filename_old = directory + "RA_I_connections_" + str(start_number) + "_.bin"
filename_new = directory + "RA_I_connections_" + str(start_number) + "_NEW.bin"

(N_RA, targets_ID_old, targets_G_old) = reading.read_connections(filename_old)
(N_RA, targets_ID_new, targets_G_new) = reading.read_connections(filename_new)

print "Old and new RA -> I connections match: ",targets_ID_old == targets_ID_new
print "Old and new RA -> I connection strengths match: ",targets_G_old == targets_G_new

###################################
# check I -> RA connections
###################################

filename_old = directory + "I_RA_connections_" + str(start_number) + "_.bin"
filename_new = directory + "I_RA_connections_" + str(start_number) + "_NEW.bin"

(N_RA, targets_ID_old, targets_G_old) = reading.read_connections(filename_old)
(N_RA, targets_ID_new, targets_G_new) = reading.read_connections(filename_new)

print "Old and new I -> RA connections match: ",targets_ID_old == targets_ID_new
print "Old and new I -> RA connection strengths match: ",targets_G_old == targets_G_new

###################################
# check active connections
###################################

filename_old = directory + "RA_RA_active_connections_" + str(start_number) + "_.bin"
filename_new = directory + "RA_RA_active_connections_" + str(start_number) + "_NEW.bin"

(N_RA, targets_ID_old, targets_G_old) = reading.read_connections(filename_old)

#==============================================================================
# print N_RA
# print targets_ID_old
# 
# for i in range(len(targets_ID_old)):
#     for j in range(len(targets_ID_old[i])):
#         print "Active supersynapse {0} -> {1}".format(i, targets_ID_old[i][j])
# 
#==============================================================================

(N_RA, targets_ID_new, targets_G_new) = reading.read_connections(filename_new)

#print N_RA

#print targets_ID_new[0]
#print targets_G_new

print "Old and new active connections match: ",targets_ID_old == targets_ID_new
print "Old and new active connection strengths match: ",targets_G_old == targets_G_new

#print targets_G_old[0]
#print targets_G_new[0]


###################################
# check super connections
###################################

filename_old = directory + "RA_RA_super_connections_" + str(start_number) + "_.bin"
filename_new = directory + "RA_RA_super_connections_" + str(start_number) + "_NEW.bin"

(N_RA, targets_ID_old, targets_G_old) = reading.read_connections(filename_old)

#==============================================================================
# print N_RA
# print targets_ID_old
# 
# for i in range(len(targets_ID_old)):
#     for j in range(len(targets_ID_old[i])):
#         print "Active supersynapse {0} -> {1}".format(i, targets_ID_old[i][j])
# 
#==============================================================================

(N_RA, targets_ID_new, targets_G_new) = reading.read_connections(filename_new)

#print N_RA
#print targets_ID_new
#print targets_G_new

print "Old and new super connections match: ",targets_ID_old == targets_ID_new
print "Old and new active connection strengths match: ",targets_G_old == targets_G_new

#print targets_G_old[0]
#print targets_G_new[0]

###################################
# check weights
###################################

filename_old = directory + "weights_" + str(start_number) + "_.bin"
filename_new = directory + "weights_" + str(start_number) + "_NEW.bin"

(N_RA, trial_number, weights_old) = reading.read_weights(filename_old)
(N_RA, trial_number, weights_new) = reading.read_weights(filename_new)

#print weights_old[0]

print "Old and new weights match: ",np.logical_and.reduce(np.logical_and.reduce(weights_old == weights_new))




#print weights_old[0]
#print weights_new[0]

###################################
# check maturation states
###################################
filename_old = directory + "mature_" + str(start_number) + "_.bin"
filename_new = directory + "mature_" + str(start_number) + "_NEW.bin"

trial_number, gaba_potential_old, firing_rate_short_old, firing_rate_long_old, remodeled_old = reading.read_maturation_info(filename_old)
trial_number, gaba_potential_new, firing_rate_short_new, firing_rate_long_new, remodeled_new= reading.read_maturation_info(filename_new)

#print gaba_potential_new

print "Old and new gaba potentials match: ",gaba_potential_old == gaba_potential_new
print "Old and new short firing rates match: ",firing_rate_short_old == firing_rate_short_new
print "Old and new long firing rates match: ",firing_rate_long_old == firing_rate_long_new
print "Old and new remodeled indicators match: ",remodeled_old == remodeled_new

#print firing_rate_short_old
#print firing_rate_long_old

#print remodeled_new


###################################
# check number of bursts in trials
###################################

filename_old = directory + "num_bursts_in_recent_trials_" + str(start_number) + "_.bin"
filename_new = directory + "num_bursts_in_recent_trials_" + str(start_number) + "_NEW.bin"

num_bursts_in_recent_trials_old = reading.read_num_bursts_in_recent_trials(filename_old)
num_bursts_in_recent_trials_new = reading.read_num_bursts_in_recent_trials(filename_new)

#print num_bursts_in_recent_trials_old[0]
print "Old and new num_bursts_in_recent_trials match: ", np.logical_and.reduce(np.logical_and.reduce(num_bursts_in_recent_trials_old == num_bursts_in_recent_trials_new))

###################################
# check replacement history
###################################

filename_old = directory + "replacement_history_" + str(start_number) + "_.bin"
filename_new = directory + "replacement_history_" + str(start_number) + "_NEW.bin"

time_from_previous_replacement_old = reading.read_replacement_history(filename_old)
time_from_previous_replacement_new = reading.read_replacement_history(filename_new)

print "Old and new time from previous replacement match: ",np.logical_and.reduce(time_from_previous_replacement_old == time_from_previous_replacement_new)


print time_from_previous_replacement_new

###################################
# check last dendritic spike times
###################################

filename_old = directory + "last_dendritic_spike_times_" + str(start_number) + "_.bin"
filename_new = directory + "last_dendritic_spike_times_" + str(start_number) + "_NEW.bin"

last_dend_spike_times_old = reading.read_last_dend_spike_times(filename_old)
last_dend_spike_times_new = reading.read_last_dend_spike_times(filename_new)

print "Old and new last dendritic spike times match: ",np.logical_and.reduce(last_dend_spike_times_old == last_dend_spike_times_new)

#print "last burst times old :", last_dend_spike_times_old
#print "last burst times new :", last_dend_spike_times_new
