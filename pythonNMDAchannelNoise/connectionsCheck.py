# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 18:23:26 2017

@author: jingroup

Script checks supersynapses and active connections
"""

import reading
import numpy as np
import reading

directory = "/home/eugene/Output/networks/sphere_170717_hodgkin/"
trial_number = 25000 # starting point
dim = 3 # dimensionality 


filename = directory + "weights_" + str(trial_number) + "_.bin"

(N_RA, trial_number, weights) = reading.read_weights(filename)

print np.where(weights[63] > 0.025)
print weights[63][82]
print weights[63][168]
print weights[63][171]
print weights[63][181]

print weights[63][17]
print weights[63][278]
