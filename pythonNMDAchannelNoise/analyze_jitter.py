# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 19:37:06 2018

@author: jingroup

Script analyzes network jitter
"""
import reading
import matplotlib.pyplot as plt
import numpy as np

TRAINING_KICK_TIME = 100.0

filename = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/test/maturationTransition2/jitter.bin"
fileMature = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/maturationTransition2/mature_14750.bin"

#filename = "/home/eugene/results/immature/clusters/test/matTrans15/jitter.bin"
#fileMature = "/home/eugene/results/immature/clusters/matTrans15/mature_10500.bin"

N, num_test_trials, \
    probability_soma_spike, average_num_soma_spikes_in_trial, mean_first_soma_spike_time, std_first_soma_spike_time,\
    probability_dend_spike, average_num_dend_spikes_in_trial, mean_first_dend_spike_time, std_first_dend_spike_time = reading.read_jitter(filename) 


(_,_,mature_indicators) = reading.read_mature_indicators(fileMature)

#mature_neurons = np.where(mature_indicators == 1)[0]
mature_neurons = range(N)

mean_first_soma_spike_time_mature = mean_first_soma_spike_time[mature_neurons]
mean_first_dend_spike_time_mature = mean_first_dend_spike_time[mature_neurons]
std_first_soma_spike_time_mature = std_first_soma_spike_time[mature_neurons]
std_first_dend_spike_time_mature = std_first_dend_spike_time[mature_neurons]
probability_soma_spike_mature = probability_soma_spike[mature_neurons]
probability_dend_spike_mature = probability_dend_spike[mature_neurons]
average_num_soma_spikes_in_trial_mature = average_num_soma_spikes_in_trial[mature_neurons]
average_num_dend_spikes_in_trial_mature = average_num_dend_spikes_in_trial[mature_neurons]

indices_sorted_soma_spikes = np.argsort(mean_first_soma_spike_time_mature)
indices_sorted_dend_spikes = np.argsort(mean_first_dend_spike_time_mature)

print N
print num_test_trials
print probability_soma_spike
#print average_num_soma_spikes_in_trial
#print mean_first_soma_spike_time

########## First somatic spike properties vs. time ##################
f = plt.figure()

ax1 = f.add_subplot(311)
ax1.plot(mean_first_soma_spike_time_mature[indices_sorted_soma_spikes] - TRAINING_KICK_TIME, probability_soma_spike_mature[indices_sorted_soma_spikes])
ax1.set_ylabel('p to produce soma spike')
ax1.set_xlim([0, 150])
ax1.set_ylim([-0.1, 1.1])

ax2 = f.add_subplot(312)
ax2.plot(mean_first_soma_spike_time_mature[indices_sorted_soma_spikes] - TRAINING_KICK_TIME, std_first_soma_spike_time_mature[indices_sorted_soma_spikes])
ax2.set_ylabel('jitter in 1st soma spike (ms)')
ax2.set_xlim([0, 150])
ax2.set_ylim([-0.1, 3.0])

ax3 = f.add_subplot(313)
ax3.plot(mean_first_soma_spike_time_mature[indices_sorted_soma_spikes] - TRAINING_KICK_TIME, average_num_soma_spikes_in_trial_mature[indices_sorted_soma_spikes])
ax3.set_ylabel('< # soma spikes >')
ax3.set_xlim([0, 150])
ax3.set_ylim([-0.1, 6.1])

########## First dendritic spike properties vs. time ##################
f = plt.figure()

ax1 = f.add_subplot(311)
ax1.plot(mean_first_dend_spike_time_mature[indices_sorted_dend_spikes] - TRAINING_KICK_TIME, probability_dend_spike_mature[indices_sorted_dend_spikes])
ax1.set_ylabel('p to produce dend spike')
ax1.set_xlim([0, 150])
ax1.set_ylim([-0.1, 1.1])

ax2 = f.add_subplot(312)
ax2.plot(mean_first_dend_spike_time_mature[indices_sorted_dend_spikes] - TRAINING_KICK_TIME, std_first_dend_spike_time_mature[indices_sorted_dend_spikes])
ax2.set_ylabel('jitter in 1st dend spike (ms)')
ax2.set_xlim([0, 50])
ax2.set_ylim([-0.1, 2.0])

ax3 = f.add_subplot(313)
ax3.plot(mean_first_dend_spike_time_mature[indices_sorted_dend_spikes] - TRAINING_KICK_TIME, average_num_dend_spikes_in_trial_mature[indices_sorted_dend_spikes])
ax3.set_ylabel('< # dend spikes >')
ax3.set_xlim([0, 50])
ax3.set_ylim([-0.1, 2.1])

########## Jitter histograms #################
f = plt.figure()
ax1 = f.add_subplot(111)

hist_soma, bin_edges_soma = np.histogram(std_first_soma_spike_time_mature[std_first_soma_spike_time_mature > 0])
hist_dend, bin_edges_dend = np.histogram(std_first_dend_spike_time_mature[std_first_dend_spike_time_mature > 0])
#print bin_edges
#print hist
bin_centers_soma = (bin_edges_soma[:-1] + bin_edges_soma[1:]) / 2.
bin_centers_dend = (bin_edges_dend[:-1] + bin_edges_dend[1:]) / 2.

ax1.step(bin_centers_soma, hist_soma, label="first somatic spike time", where="pre")
#ax1.step(bin_centers_dend, hist_dend, label="first dendritic spike time", where="pre")

ax1.set_xlabel('Jitter (ms)')
ax1.set_ylabel('# of neurons')

ax1.legend()

plt.show()