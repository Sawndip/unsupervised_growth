# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 16:22:09 2018

@author: jingroup

Script analyzes connections that emerge between training HVC-RA neurons
"""

import reading
import matplotlib.pyplot as plt
import os
import numpy as np
import utils



dirname = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays/"

trial_number = 100

fileTraining = os.path.join(dirname, "training_neurons.bin")
fileDend = os.path.join(dirname, "spike_times_dend_" + str(trial_number) + ".bin")
fileSoma = os.path.join(dirname, "spike_times_soma_" + str(trial_number) + ".bin")
fileWeights = os.path.join(dirname, "weights_" + str(trial_number) + ".bin")

training_neurons = sorted(reading.read_training_neurons(fileTraining))

soma_spike_times_training = utils.get_event_times_neurons(fileSoma, training_neurons)

first_somatic_spikes_training = [spikes[0] for spikes in soma_spike_times_training]

    
print first_somatic_spikes_training
spread = max(first_somatic_spikes_training) - min(first_somatic_spikes_training)

print "Spread of training neurons = {0} ms".format(spread)

trialStep = 50
endTrial = 200

input_weights = utils.get_input_weight_time_sequence(dirname, endTrial, trialStep, training_neurons, training_neurons)

print input_weights[:,-1]

ind_max = np.argmax(input_weights[:,-1])

print ind_max

print input_weights.shape[1]

plt.figure()

plt.plot([i*trialStep for i in range(input_weights.shape[1])], input_weights[ind_max])
plt.xlabel('Time(# of trials)')
plt.ylabel('Total input weight')

plt.show()
#(_, _, spike_times_dend, neuron_fired_dend) = reading.read_time_info(fileDend)

#dend_spike_times_training = get_event_times_neurons(fileDend, training_neurons)


#print training_neurons

#print soma_spike_times_training
#print dend_spike_times_training


#print spike_times_soma
#print neuron_fired_soma