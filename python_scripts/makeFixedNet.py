# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 16:19:13 2017

@author: jingroup

Script creates fixed.net text file for pajek from binary files with neuronal coordinates
of HVC(RA) neurons and synapses betweeen them
"""
import reading
import os
import space

dirname = "/home/eugene/Output/networks/sphere_170717_hodgkin/"

#dirname = "/home/eugene/results/noDelays/replacement/sphere/180717_lionx_3/"
arrangement = "sphere" # spatial arrangement of neurons
dim = space.num_coordinates[arrangement] # dimensionality
#dirname = "/home/eugene/Output/networks/140717_huxley/"

trial_number = 21600

file_RA_xy = os.path.join(dirname, "RA_xy_" + str(trial_number) + "_.bin")
file_I_xy = os.path.join(dirname, "I_xy.bin")

RA2I = os.path.join(dirname, "RA_I_connections_" + str(trial_number) + "_.bin")
I2RA = os.path.join(dirname, "I_RA_connections_" + str(trial_number) + "_.bin")

file_training = os.path.join(dirname, "training_neurons.bin")
file_pajek = os.path.join(dirname, "fixed_" + str(trial_number) + "_.net")

(N_RA, RA_targets, RA_targets_G) = reading.read_connections(RA2I)
(N_I, I_targets, I_targets_G) = reading.read_connections(I2RA)

coord_RA = reading.read_coordinates(dim, file_RA_xy)
coord_I = reading.read_coordinates(dim, file_I_xy)

#print targets_ID
#print targets_G
if file_training:
    training_neurons = reading.read_training_neurons(file_training)
else:
    training_neurons = []
    
print "Training neurons: ",training_neurons

with open(file_pajek, 'w') as f:
    f.write("*Vertices {0}\n".format(N_RA+N_I))
    
    if dim == 2:
        for i in range(N_RA):
            if i in training_neurons:    
                f.write('{0} "{1}" {2} {3} ic Green\n'.format(i+1, i, coord_RA[i][0], coord_RA[i][1]))
            else:    
                f.write('{0} "{1}" {2} {3} ic Yellow\n'.format(i+1, i, coord_RA[i][0], coord_RA[i][1]))
                
        for i in range(N_RA, N_RA+N_I):
            f.write('{0} "{1}" {2} {3} ic Red\n'.format(i+1, i, coord_I[i-N_RA][0], coord_I[i-N_RA][1]))
            
            
    if dim == 3:
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
            

