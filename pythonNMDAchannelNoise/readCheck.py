# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 16:44:37 2017

@author: jingroup

Script checks that data was read correctly after chain growth was restarted
"""
import reading
import numpy as np

directory = "/home/eugene/Output/networks/networkTest/"
start_number = 5 # starting point

###################################
# check RA -> I connections
###################################

filename_old = directory + "RA_I_connections_" + str(start_number) + "_.bin"
filename_new = directory + "RA_I_connections_" + str(start_number) + "_NEW.bin"

(N_RA, targets_ID_old, targets_G_old) = reading.read_connections(filename_old)
(N_RA, targets_ID_new, targets_G_new) = reading.read_connections(filename_new)

print targets_ID_old == targets_ID_new
print targets_G_old == targets_G_new

###################################
# check I -> RA connections
###################################

filename_old = directory + "I_RA_connections_" + str(start_number) + "_.bin"
filename_new = directory + "I_RA_connections_" + str(start_number) + "_NEW.bin"

(N_RA, targets_ID_old, targets_G_old) = reading.read_connections(filename_old)
(N_RA, targets_ID_new, targets_G_new) = reading.read_connections(filename_new)

print targets_ID_old == targets_ID_new
print targets_G_old == targets_G_new

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

print N_RA
#print targets_ID_new
#print targets_G_new

print targets_ID_old == targets_ID_new

print targets_G_old == targets_G_new

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

print N_RA
#print targets_ID_new
#print targets_G_new

print targets_ID_old == targets_ID_new

print targets_G_old == targets_G_new

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

print np.logical_and.reduce(weights_old == weights_new)




#print weights_old[0]
#print weights_new[0]

###################################
# check maturation states
###################################
filename_old = directory + "mature_" + str(start_number) + "_.bin"
filename_new = directory + "mature_" + str(start_number) + "_NEW.bin"

trial_number, gaba_potential_old, firing_rate_old, remodeled_old, mature_old = reading.read_maturation_info(filename_old)
trial_number, gaba_potential_new, firing_rate_new, remodeled_new, mature_new = reading.read_maturation_info(filename_new)

#print gaba_potential_new

print gaba_potential_old == gaba_potential_new
print firing_rate_old == firing_rate_new
print remodeled_old == remodeled_new
print mature_old == mature_new

###################################
# check number of bursts in trials
###################################

filename_old = directory + "num_bursts_in_recent_trials_" + str(start_number) + "_.bin"
filename_new = directory + "num_bursts_in_recent_trials_" + str(start_number) + "_NEW.bin"

num_bursts_in_recent_trials_old = reading.read_num_bursts_in_recent_trials(filename_old)
num_bursts_in_recent_trials_new = reading.read_num_bursts_in_recent_trials(filename_new)

print num_bursts_in_recent_trials_old == num_bursts_in_recent_trials_new

###################################
# check replacement history
###################################

filename_old = directory + "replacement_history_" + str(start_number) + "_.bin"
filename_new = directory + "replacement_history_" + str(start_number) + "_NEW.bin"

time_from_previous_replacement_old = reading.read_replacement_history(filename_old)
time_from_previous_replacement_new = reading.read_replacement_history(filename_new)

print time_from_previous_replacement_old == time_from_previous_replacement_new

