# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 16:19:13 2017

@author: jingroup

Script creates .net text file for pajek from binary files with neuronal coordinates
of HVC(RA) neurons and synapses betweeen them
"""
import reading
import os

dirname = "/home/eugene/Output/networks/sphere_170717_hodgkin/"
dim = 3
#dirname = "/home/eugene/Output/networks/140717_huxley/"

trial_number = 9700

file_xy = os.path.join(dirname, "RA_xy_" + str(trial_number) + "_.bin")
file_super = os.path.join(dirname, "RA_RA_super_connections_" + str(trial_number) + "_.bin")

#file_training = None
file_training = os.path.join(dirname, "training_neurons.bin")
file_pajek = os.path.join(dirname, "super_" + str(trial_number) + "_.net")

(N_RA, targets_ID, targets_G) = reading.read_connections(file_super)
coord = reading.read_coordinates(dim, file_xy)

#print targets_ID
#print targets_G
if file_training:
    training_neurons = reading.read_training_neurons(file_training)
else:
    training_neurons = []
    
print "Training neurons: ",training_neurons

with open(file_pajek, 'w') as f:
    f.write("*Vertices {0}\n".format(N_RA))
    
    if dim == 2:
        SIDE = 100.0 # side of square modeling HVC
        
    
        for i in range(N_RA):
            if i in training_neurons:    
                f.write('{0} "{1}" {2} {3} ic Green\n'.format(i+1, i, coord[i][0] / SIDE, coord[i][1] / SIDE))
            else:    
                f.write('{0} "{1}" {2} {3} ic Yellow\n'.format(i+1, i, coord[i][0] / SIDE, coord[i][1] / SIDE))
            
    if dim == 3:
        for i in range(N_RA):
            if i in training_neurons:    
                f.write('{0} "{1}" {2} {3} {4} ic Green\n'.format(i+1, i, coord[i][0], coord[i][1], coord[i][2]))
            else:    
                f.write('{0} "{1}" {2} {3} {4} ic Yellow\n'.format(i+1, i, coord[i][0], coord[i][1], coord[i][2]))
    
    f.write("*Arcs\n".format(N_RA))
    
    for i, targets in enumerate(targets_ID):
        for j, target in enumerate(targets):
            f.write('{0} {1} {2}\n'.format(i+1, target+1, targets_G[i][j]))
            

