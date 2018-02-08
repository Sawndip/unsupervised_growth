# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 10:38:19 2018

@author: jingroup

Script creates pajek .net file for network topology
"""
import reading
import os
import space

dirname = "/home/eugene/Output/networks/chainGrowth/test/"



file_RA_xy = os.path.join(dirname, "RA_xy.bin")
file_I_xy = os.path.join(dirname, "I_xy.bin")

RA2I = os.path.join(dirname, "RA_I_connections.bin")
I2RA = os.path.join(dirname, "I_RA_connections.bin")

file_training = os.path.join(dirname, "training_neurons_clustered.bin")
file_pajek = os.path.join(dirname, "network_topology_clustered.net")

(N_RA, RA_targets, RA_targets_G, _, _) = reading.read_connections(RA2I)
(N_I, I_targets, I_targets_G, _, _) = reading.read_connections(I2RA)

coord_RA = reading.read_coordinates(file_RA_xy)
coord_I = reading.read_coordinates(file_I_xy)

#print targets_ID
#print targets_G
if file_training:
    training_neurons = reading.read_training_neurons(file_training)
else:
    training_neurons = []
    
print "Training neurons: ",training_neurons

with open(file_pajek, 'w') as f:
    f.write("*Vertices {0}\n".format(N_RA+N_I))
    
            
            
    for i in range(N_RA):
        if i in training_neurons:    
            f.write('{0} "{1}" {2} {3} {4} ic Green\n'.format(i+1, i, coord_RA[i][0], coord_RA[i][1], coord_RA[i][2]))
        else:    
            f.write('{0} "{1}" {2} {3} {4} ic Yellow\n'.format(i+1, i, coord_RA[i][0], coord_RA[i][1], coord_RA[i][2]))

    for i in range(N_RA, N_RA+N_I):
        f.write('{0} "{1}" {2} {3} {4} ic Red\n'.format(i+1, i, coord_I[i-N_RA][0], coord_I[i-N_RA][1], coord_I[i-N_RA][2]))

    f.write("*Arcs\n".format(N_RA))
    
    # write targets of HVC(RA) neurons
    for i, targets in enumerate(RA_targets):
        for j, target in enumerate(targets):
            f.write('{0} {1} {2} c Green\n'.format(i+1, N_RA+target+1, RA_targets_G[i][j]))
            
     # write targets of HVC(I) neurons
    for i, targets in enumerate(I_targets):
        for j, target in enumerate(targets):
            f.write('{0} {1} {2} c Red\n'.format(N_RA+i+1, target+1, I_targets_G[i][j]))
            

