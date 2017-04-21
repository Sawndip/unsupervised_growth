# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 15:30:42 2017

@author: jingroup

Script compares fixed networks
"""

import reading
import numpy as np

###################################
# check HVC(RA) neuron coordinates
###################################
filename_old = "/home/eugene/Output/networks/dispersed/RA_xy.bin"
filename_new = "/home/eugene/Output/networks/clustered/RA_xy.bin"

(xx_old, yy_old) = reading.read_coordinates(filename_old)
(xx_new, yy_new) = reading.read_coordinates(filename_new)

print [i for i,j in zip(xx_old, xx_new) if i!=j]
print xx_old == xx_new
print yy_old == yy_new

###################################
# check HVC(I) neuron coordinates
###################################
filename_old = "/home/eugene/Output/networks/dispersed/I_xy.bin"
filename_new = "/home/eugene/Output/networks/clustered/I_xy.bin"

(xx_old, yy_old) = reading.read_coordinates(filename_old)
(xx_new, yy_new) = reading.read_coordinates(filename_new)

print xx_old == xx_new
print yy_old == yy_new

###################################
# check RA -> I connections
###################################

filename_old = "/home/eugene/Output/networks/dispersed/RA_I_connections.bin"
filename_new = "/home/eugene/Output/networks/clustered/RA_I_connections.bin"

(N_RA, targets_ID_old, targets_G_old) = reading.read_connections(filename_old)
(N_RA, targets_ID_new, targets_G_new) = reading.read_connections(filename_new)

print targets_ID_old[0] == targets_ID_new[0]
print targets_G_old == targets_G_new

###################################
# check I -> RA connections
###################################

filename_old = "/home/eugene/Output/networks/dispersed/I_RA_connections.bin"
filename_new = "/home/eugene/Output/networks/clustered/I_RA_connections.bin"

(N_RA, targets_ID_old, targets_G_old) = reading.read_connections(filename_old)
(N_RA, targets_ID_new, targets_G_new) = reading.read_connections(filename_new)

print targets_ID_old[0] == targets_ID_new[0]
print targets_G_old == targets_G_new

