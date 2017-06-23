# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 16:19:13 2017

@author: jingroup

Script creates .net text file for pajek from binary files with neuronal coordinates
of HVC(RA) neurons and synapses betweeen them
"""
import reading
import os

dirname = "/home/eugene/results/noDelays/dispersed/dispersed_1/120617_lionx_2/"

file_xy = os.path.join(dirname, "RA_xy_initial.bin")
file_super = os.path.join(dirname, "RA_RA_super_connections.bin")
file_training = os.path.join(dirname, "training_neurons.bin")
file_pajek = os.path.join(dirname, "super.net")

SIDE = 100.0 # side of square modeling HVC
N_TR = 4 # number of training neurons

(xx, yy) = reading.read_coordinates(file_xy)
(N_RA, targets_ID, targets_G) = reading.read_connections(file_super)

#print targets_ID
#print targets_G

training_neurons = reading.read_training_neurons(file_training)

print "Training neurons: ",training_neurons

with open(file_pajek, 'w') as f:
    f.write("*Vertices {0}\n".format(N_RA))
    
    for i in range(N_RA):
        if i in training_neurons:    
            f.write('{0} "{1}" {2} {3} ic Green\n'.format(i+1, i, xx[i] / SIDE, yy[i] / SIDE))
        else:    
            f.write('{0} "{1}" {2} {3} ic Yellow\n'.format(i+1, i, xx[i] / SIDE, yy[i] / SIDE))
        

    f.write("*Arcs\n".format(N_RA))
    
    for i, targets in enumerate(targets_ID):
        for j, target in enumerate(targets):
            f.write('{0} {1} {2}\n'.format(i+1, target+1, targets_G[i][j]))
            

