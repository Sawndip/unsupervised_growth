# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 16:44:37 2017

@author: jingroup

Script reads supersynapses from file
"""
import reading
import numpy as np

###################################
# check RA -> I connections
###################################

filename_old = "/home/eugene/Output/networks/gabaMaturation160417/RA_I_connections.bin"
filename_new = "/home/eugene/Output/networks/gabaMaturation160417/RA_I_connections_NEW.bin"

(N_RA, targets_ID_old, targets_G_old) = reading.read_connections(filename_old)
(N_RA, targets_ID_new, targets_G_new) = reading.read_connections(filename_new)

print targets_ID_old == targets_ID_new
print targets_G_old == targets_G_new

###################################
# check I -> RA connections
###################################

filename_old = "/home/eugene/Output/networks/gabaMaturation160417/I_RA_connections.bin"
filename_new = "/home/eugene/Output/networks/gabaMaturation160417/I_RA_connections_NEW.bin"

(N_RA, targets_ID_old, targets_G_old) = reading.read_connections(filename_old)
(N_RA, targets_ID_new, targets_G_new) = reading.read_connections(filename_new)

print targets_ID_old == targets_ID_new
print targets_G_old == targets_G_new

###################################
# check active connections
###################################

filename_old = "/home/eugene/Output/networks/gabaMaturation160417/RA_RA_active_connections_4600_.bin"
filename_new = "/home/eugene/Output/networks/gabaMaturation160417/RA_RA_active_connections_NEW.bin"

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

filename_old = "/home/eugene/Output/networks/gabaMaturation160417/RA_RA_super_connections.bin"
filename_new = "/home/eugene/Output/networks/gabaMaturation160417/RA_RA_super_connections_NEW.bin"

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

filename_old = "/home/eugene/Output/networks/gabaMaturation160417/weights.bin"
filename_new = "/home/eugene/Output/networks/gabaMaturation160417/weights_NEW.bin"

(N_RA, trial_number, weights_old) = reading.read_weights(filename_old)
(N_RA, trial_number, weights_new) = reading.read_weights(filename_new)

#print weights_old[0]

print np.logical_and.reduce(weights_old == weights_new)




#print weights_old[0]
#print weights_new[0]

###################################
# check maturation states
###################################

filename_old = "/home/eugene/Output/networks/gabaMaturation160417/mature.bin"
filename_new = "/home/eugene/Output/networks/gabaMaturation160417/mature_NEW.bin"

trial_number, gaba_potential_old, firing_rate_old, remodeled_old, mature_old = reading.read_maturation_info(filename_old)
trial_number, gaba_potential_new, firing_rate_new, remodeled_new, mature_new = reading.read_maturation_info(filename_new)

#print gaba_potential_new

print gaba_potential_old == gaba_potential_new
print firing_rate_old == firing_rate_new
print remodeled_old == remodeled_new
print mature_old == mature_new


