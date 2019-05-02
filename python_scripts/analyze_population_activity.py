# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 20:41:16 2018

@author: jingroup

Script analyzes HVC-RA neuron population activity in the network
"""
import matplotlib.pyplot as plt
import reading
import os
import numpy as np
import utils


dirname = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays/"
trial_number = 350
fileActivityHistory = os.path.join(dirname, "activity_history_" + str(trial_number) + ".bin")
fileTraining = os.path.join(dirname, "training_neurons.bin")

(_, _, activity_history) = reading.read_activity_history(fileActivityHistory) 
training_neurons = reading.read_training_neurons(fileTraining)

np.set_printoptions(threshold=np.nan)

#print activity_history[0][0:10]
#print activity_history[0:100,0:100]

#print activity_history

# check that activity is reasonable
print "Negative numbers in activity: ",np.any(activity_history < 0)
print "Large positive numbers in activity: ",np.any(activity_history > 10)


num_neurons = activity_history.shape[0]
history_size = activity_history.shape[1]

fire_indicators = np.ones_like(activity_history, np.float32)
fire_indicators[activity_history == 0] = 0

#print fire_indicators[156]
#==============================================================================
# if trial_number >= history_size:    
#     num_points = trial_number / window_size
# else:
#     num_points = history_size / window_size
#==============================================================================
    
probability_to_spike = np.empty(num_neurons, np.float32)

for i in range(num_neurons):
    probability_to_spike[i] = float(np.sum(fire_indicators[i])) / float(trial_number)

nbins = 50

plt.figure()
hist, bin_edges = np.histogram(probability_to_spike, bins=nbins)
print bin_edges
print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
plt.step(bin_centers, hist, label="active", where="pre")

plt.ylabel('# of neurons')
plt.xlabel('Population firing rate')
plt.xlim([0, trial_number + 10])

plt.show()
