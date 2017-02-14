# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 16:44:37 2017

@author: jingroup

Script reads supersynapses from file
"""
import reading

filename_old = "/home/eugene/Output/networks/gabaMaturation300117/RA_RA_super_connections.bin"
filename_new = "/home/eugene/Output/RA_RA_super_connections_NEW.bin"

(N_RA, targets_ID_old, targets_G_old) = reading.read_connections(filename_old)

print N_RA
print targets_ID_old

for i in range(len(targets_ID_old)):
    for j in range(len(targets_ID_old[i])):
        print "Active supersynapse {0} -> {1}".format(i, targets_ID_old[i][j])


(N_RA, targets_ID_new, targets_G_new) = reading.read_connections(filename_new)

print N_RA
print targets_ID_new

print targets_ID_old == targets_ID_new


filename_old = "/home/eugene/Output/networks/gabaMaturation300117/weights.bin"
filename_new = "/home/eugene/Output/weights_NEW.bin"

(N_RA, weights_old) = reading.read_weights(filename_old)
(N_RA, weights_new) = reading.read_weights(filename_new)

print weights_old

print weights_old == weights_new

filename_old = "/home/eugene/Output/networks/gabaMaturation300117/mature.bin"
filename_new = "/home/eugene/Output/mature_NEW.bin"

trial_number, gaba_potential_old, firing_rate_old, remodeled_old, mature_old = reading.read_maturation_info(filename_old)
trial_number, gaba_potential_new, firing_rate_new, remodeled_new, mature_new = reading.read_maturation_info(filename_new)

print gaba_potential_old == gaba_potential_new
print firing_rate_old == firing_rate_new
print remodeled_old == remodeled_new
print mature_old == mature_new
