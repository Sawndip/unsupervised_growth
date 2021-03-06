# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 18:20:25 2016

@author: jingroup
"""

"""
Script plots all inhibitory spikes occured in the network
"""

import reading
import matplotlib.pyplot as plt
import numpy as np
import math

N = 550 # number of neurons
TRIAL_DURATION = 500

#fileInterneuron = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays5/spike_times_interneuron_300.bin"
#fileInterneuron = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/test1/test_spike_times_interneuron_10.bin"
#fileInterneuron = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/matTrans2_network2000RA550I/spike_times_interneuron_0.bin"
#fileInterneuron = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/split1/spike_times_interneuron_3500.bin"

fileInterneuron = "/home/eugene/results/immature/clusters/matTrans36/spike_times_interneuron_8400.bin"
#fileInterneuron = "/home/eugene/results/immature/clusters/test/matTrans30/test_spike_times_interneuron_10.bin"

(trial_number, simulation_time, spike_times_interneuron, neuron_fired) = reading.read_time_info(fileInterneuron)

print spike_times_interneuron

spike_times_interneuron = [t for sublist in list(spike_times_interneuron) for t in sublist]
neuron_fired = [ind for sublist in list(neuron_fired) for ind in sublist]

ordered_spikes, ordered_neurons = zip(*sorted(zip(spike_times_interneuron, neuron_fired)))

#for i in xrange(len(ordered_spikes)):
 #   print ordered_spikes[i], ordered_neurons[i]

f = plt.figure()
ax1 = f.add_subplot(111)

for i in range(len(spike_times_interneuron)):
    ax1.vlines(spike_times_interneuron[i], neuron_fired[i]-0.5, neuron_fired[i]+0.5)

#ax.scatter(spike_times, random_ID)
#ax1.set_yticks(random_ID)

#plt.yticks(random_ID, neuron_fired)
ax1.set_xlim([-5, TRIAL_DURATION])
ax1.set_ylabel("neuron ID")
ax1.set_xlabel("Time (ms)")
ax1.set_title("Interneuron spikes")
ax1.set_ylim([0, N])

plt.show()