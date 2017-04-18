# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 12:40:01 2017

@author: jingroup

Script analyzes population activity of HVC(RA) neurons
"""

import reading
import numpy as np
import matplotlib.pyplot as plt


def get_population_firing_rate(filename):
    """
    Calculate population firing rate from data of firing rates of all neurons 
    in the pool
    """
    (_, t, _, _, _, firing_rate) = reading.read_maturation_time_sequence(filename)

    num_timepoints = len(firing_rate[0])
    num_neurons = len(firing_rate)

    firing_rate_total = np.zeros(num_timepoints, dtype=np.float32)

    for i in range(num_timepoints):
        for j in range(num_neurons):
            firing_rate_total[i] += firing_rate[j][i]
            
        firing_rate_total[i] /= float(num_neurons)
     
    return (t, firing_rate_total)           


fileNoisy = "/home/eugene/Output/networks/gabaMaturation120417/maturation_time_sequence.bin"
fileSilent = "/home/eugene/Output/networks/gabaMaturation110417_3/maturation_time_sequence.bin"


t_noisy, firing_rate_noisy = get_population_firing_rate(fileNoisy)
t_silent, firing_rate_silent = get_population_firing_rate(fileSilent)

f = plt.figure()

ax = f.add_subplot(111)

ax.plot(t_noisy, firing_rate_noisy, label="spontaneous activity")
ax.plot(t_silent, firing_rate_silent, label=" no spontaneous activity")


ax.set_xlabel("time (# trials)")
ax.set_ylabel("population firing rate (Hz)")

plt.legend()
plt.show()
#==============================================================================
# files = os.listdir(dataDir)
# 
# num_bursts_total = []
# 
# for f in files:
#     if "spike_times_dend_" in f:
#         filename = os.path.join(dataDir, f)
#         print filename       
#         (trial_number, simulation_time, spike_times, neuron_fired) = reading.read_time_info(filename)
#         num_bursts_total.append(len(spike_times))
# 
# print num_bursts_total
# 
#==============================================================================
