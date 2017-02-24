# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 16:19:13 2017

@author: jingroup

Script creates .net text file for pajek from binary files with neuronal coordinates
of HVC(RA) neurons and synapses betweeen them
"""
import reading

file_xy = "/home/eugene/Output/networks/gabaMaturation150217/RA_xy.bin"
file_super = "/home/eugene/Output/networks/gabaMaturation150217/RA_RA_super_connections.bin"

file_pajek = "/home/eugene/Output/networks/gabaMaturation150217/super.net"

SIDE = 100.0 # side of square modeling HVC
N_TR = 4 # number of training neurons

(xx, yy) = reading.read_coordinates(file_xy)
(N_RA, targets_ID, targets_G) = reading.read_connections(file_super)

#print targets_ID
#print targets_G

with open(file_pajek, 'w') as f:
    f.write("*Vertices {0}\n".format(N_RA))
    
    for i in range(N_TR):
        f.write('{0} "{1}" {2} {3} ic Green\n'.format(i+1, i, xx[i] / SIDE, yy[i] / SIDE))
        
    for i in range(N_TR, N_RA):
        f.write('{0} "{1}" {2} {3} ic Yellow\n'.format(i+1, i, xx[i] / SIDE, yy[i] / SIDE))
        

    f.write("*Arcs\n".format(N_RA))
    
    for i, targets in enumerate(targets_ID):
        for j, target in enumerate(targets):
            f.write('{0} {1} {2}\n'.format(i+1, target+1, targets_G[i][j]))
            

