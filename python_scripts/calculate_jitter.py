#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 11:50:46 2018

@author: jingroup

Script calculates jitter in first somatic and dendritic spike times of mature neurons
"""
import utils
import reading
import matplotlib.pyplot as plt
import os
import numpy as np

num_trials = 10

dirname = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/noImmatureOut7/"
#dirname = "/home/eugene/results/immature/clusters/4/"

trial_number = 14000
trialStep = 50
N_RA = 1000

training_neurons = [195, 228, 860, 461, 525, 146, 726, 873, 252, 893]
#training_neurons = [23, 127, 179, 196, 221, 228, 285, 291, 438, 493, 592, 669, 726, 807, 822, 824, 850, 886, 909, 964]

fileMature = os.path.join(dirname, "mature_" + str(trial_number - trialStep*num_trials) + ".bin")

(_, _, mature_indicators) = reading.read_mature_indicators(fileMature)

# track spike times of neurons that were mature at the time of the earliest trial
mature_neurons = np.where(mature_indicators == 1)[0]

print "Mature neurons: ",mature_neurons

spike_times_sequence = utils.get_first_spike_time_sequence_relative_to_training(dirname, trial_number, trialStep, N_RA, mature_neurons, training_neurons)

print spike_times_sequence[:, -num_trials:]

mean_spike_time = np.mean(spike_times_sequence[:, -num_trials:], axis=1)
std_spike_time = np.std(spike_times_sequence[:, -num_trials:], axis=1)

print mean_spike_time

plt.figure()
plt.hist(std_spike_time)
plt.xlabel('std of spike time (ms)')
plt.ylabel('Counts')

plt.figure()

mean_spike_time, std_spike_time = zip(*sorted(zip(mean_spike_time, std_spike_time)))
plt.plot(mean_spike_time, std_spike_time)
plt.xlabel('mean spike time (ms)')
plt.ylabel('std of spike time (ms)')

plt.show()