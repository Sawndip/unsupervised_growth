# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 10:39:00 2018

@author: jingroup

Script analyzes convergence of inhibitory inputs to HVC-RA neurons
"""
import reading
import matplotlib.pyplot as plt
import os
from collections import Counter
import utils
import numpy as np


trial_number = 0
dirname = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays6/"

fileRA2I = os.path.join(dirname, "RA_I_connections_" + str(trial_number) + ".bin")
fileI2RA = os.path.join(dirname, "I_RA_connections_" + str(trial_number) + ".bin")

fileTraining = os.path.join(dirname, "training_neurons.bin")

(N_RA, targets_id_RA2I, weights_RA2I, _, _) = reading.read_connections(fileRA2I)
(N_I, targets_id_I2RA, weights_I2RA, _, _) = reading.read_connections(fileI2RA)

training_neurons = reading.read_training_neurons(fileTraining)


utils.find_inhibitory_effect_of_hvcra(training_neurons, targets_id_RA2I, weights_RA2I, targets_id_I2RA, weights_I2RA)

### Compare with number of inputs from all interneurons ###
pool_neurons = range(N_RA)
interneurons = range(N_I)

pool_neurons_sorted, num_inhibitory_inputs, input_inhibitory_weight = \
    utils.find_convergence_of_inhibitory_input(interneurons, pool_neurons, targets_id_I2RA, weights_I2RA)
    
#print "Pool neurons that receive inhibition: ",pool_neurons_sorted
print "Number of pool neurons that receive inhibition: ",len(pool_neurons_sorted)
print "Number of inhibitory inputs received: ",num_inhibitory_inputs
#print "Input inhibitory weights: ",input_inhibitory_weight
nbins = 50

f = plt.figure()


ax1 = f.add_subplot(211)
hist, bin_edges = np.histogram(num_inhibitory_inputs, bins=nbins)
#print bin_edges
#print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
ax1.step(bin_centers, hist, where="pre")
ax1.set_xlabel('# of inhibitory inputs')
ax1.set_ylabel('# of neurons')
ax1.set_title('Convergence of inputs from all interneurons')

ax2 = f.add_subplot(212)    
hist, bin_edges = np.histogram(input_inhibitory_weight, bins=nbins)
#print bin_edges
#print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
ax2.step(bin_centers, hist, where="pre")
ax2.set_xlabel('Input inhibitory weight')
ax2.set_ylabel('# of interneurons') 
    
plt.show()