# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 16:44:37 2017

@author: jingroup

Script reads supersynapses from file
"""
import reading

filename_super = "/home/eugene/Output/networks/test170217/RA_RA_super_connections.bin"
filename_active = "/home/eugene/Output/networks/test170217/RA_RA_active_connections.bin"


(N_RA, targets_ID_super, targets_G_super) = reading.read_connections(filename_super)

print N_RA
print targets_ID_super

print "\nSupersynapses:\n"

for i in range(len(targets_ID_super)):
    for j in range(len(targets_ID_super[i])):
        print "Supersynapse {0} -> {1}".format(i, targets_ID_super[i][j])

(N_RA, targets_ID_active, targets_G_active) = reading.read_connections(filename_active)

print "\nActive synapses:\n"

for i in range(len(targets_ID_active)):
    for j in range(len(targets_ID_active[i])):
        print "Active synapse {0} -> {1}".format(i, targets_ID_active[i][j])

