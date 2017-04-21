# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 14:12:25 2016

@author: jingroup
"""

"""
Script plots all spikes occured in the network
"""

import reading
import matplotlib.pyplot as plt
import numpy as np
import math

N = 300 # number of neurons

fileDend = "/mnt/hodgkin_home/eugene/lionX/gabaMaturation180417_8/spike_times_dend.bin"
fileSoma = "/mnt/hodgkin_home/eugene/lionX/gabaMaturation180417_8/spike_times_soma.bin"

TRIAL_DURATION = 1000

(trial_number, simulation_time, spike_times_dend, neuron_fired_dend) = reading.read_time_info(fileDend)
(trial_number, simulation_time, spike_times_soma, neuron_fired_soma) = reading.read_time_info(fileSoma)

#print "Dedritic spikes: ", neuron_fired_dend
#print "Dendritic spike times: ", spike_times_dend

#print "Somatic spikes: ", neuron_fired_soma
#print "Somatic spike times: ", spike_times_soma




spike_times_soma = [t for sublist in list(spike_times_soma) for t in sublist]
neuron_fired_soma = [ind for sublist in list(neuron_fired_soma) for ind in sublist]

spike_times_dend = [t for sublist in list(spike_times_dend) for t in sublist]
neuron_fired_dend = [ind for sublist in list(neuron_fired_dend) for ind in sublist]

ordered_soma_spikes, ordered_neurons = zip(*sorted(zip(spike_times_soma, neuron_fired_soma)))

for i in xrange(len(ordered_soma_spikes)):
    print ordered_soma_spikes[i], ordered_neurons[i]

f = plt.figure()
ax1 = f.add_subplot(211)

for i in range(len(spike_times_dend)):
    ax1.vlines(spike_times_dend[i], neuron_fired_dend[i]-0.5, neuron_fired_dend[i]+0.5)

#ax.scatter(spike_times, random_ID)
#ax1.set_yticks(random_ID)

#plt.yticks(random_ID, neuron_fired)
ax1.set_xlim([-5, TRIAL_DURATION])
ax1.set_ylabel("neuron ID")
ax1.set_xlabel("relative spike time (ms)")
ax1.set_title("Dendritic spikes")
ax1.set_ylim([0, N])


ax2 = f.add_subplot(212)

for i in range(len(spike_times_soma)):
    ax2.vlines(spike_times_soma[i], neuron_fired_soma[i]-0.5, neuron_fired_soma[i]+0.5)

#ax.scatter(spike_times, random_ID)
#ax1.set_yticks(random_ID)

#plt.yticks(random_ID, neuron_fired)
ax2.set_xlim([-5, TRIAL_DURATION])
ax2.set_ylim([0, N])
ax2.set_ylabel("neuron ID")
ax2.set_xlabel("relative spike time (ms)")
ax2.set_title("Somatic spikes")

plt.show()