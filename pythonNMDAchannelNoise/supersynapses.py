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
