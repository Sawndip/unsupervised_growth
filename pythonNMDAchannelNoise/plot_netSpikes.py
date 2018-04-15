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

N = 1000 # number of neurons

#fileDend = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays6/spike_times_dend_1900.bin"
#fileSoma = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays6/spike_times_soma_1900.bin"

fileDend = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/maturationTransition1/spike_times_dend_35.bin"
fileSoma = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/maturationTransition1/spike_times_soma_35.bin"

#fileDend = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/test/noImmatureOut4/test_spike_times_dend_5.bin"
#fileSoma = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/test/noImmatureOut4/test_spike_times_soma_5.bin"


TRIAL_DURATION = 500

(trial_number, simulation_time, spike_times_dend, neuron_fired_dend) = reading.read_time_info(fileDend)
(trial_number, simulation_time, spike_times_soma, neuron_fired_soma) = reading.read_time_info(fileSoma)

print "Number of spiked neurons: ",len(spike_times_soma)

for spikes, neuron_id in zip(spike_times_soma, neuron_fired_soma):
    if len(spikes) > 6:
        print "Neuron {0} produced {1} spikes".format(neuron_id[0], len(spikes))

#print "Dedritic spikes: ", neuron_fired_dend
#print "Dendritic spike times: ", spike_times_dend

#print "Somatic spikes: ", neuron_fired_soma
#print "Somatic spike times: ", spike_times_soma

ordered_soma_spikes_raw, ordered_soma_raw = zip(*sorted(zip(spike_times_soma, neuron_fired_soma)))
ordered_dend_spikes_raw, ordered_dend_raw = zip(*sorted(zip(spike_times_dend, neuron_fired_dend)))

print [ spikes[0] for spikes in ordered_soma_spikes_raw[:] ]
print [ neurons[0] for neurons in ordered_soma_raw[:] ]

spike_times_soma = [t for sublist in list(spike_times_soma) for t in sublist]
neuron_fired_soma = [ind for sublist in list(neuron_fired_soma) for ind in sublist]

spike_times_dend = [t for sublist in list(spike_times_dend) for t in sublist]
neuron_fired_dend = [ind for sublist in list(neuron_fired_dend) for ind in sublist]

#ordered_soma_spikes, ordered_neurons = zip(*sorted(zip(spike_times_soma, neuron_fired_soma)))

ordered_soma_spikes = [t for sublist in list(ordered_soma_spikes_raw) for t in sublist]
ordered_soma_neurons = [t for sublist in list(ordered_soma_raw) for t in sublist]

ordered_id_soma_list = [neuron_id[0] for neuron_id in ordered_soma_raw]
id_soma_map = {i : j for i,j in zip(ordered_id_soma_list, range(len(ordered_id_soma_list)))}

ordered_dend_spikes = [t for sublist in list(ordered_dend_spikes_raw) for t in sublist]
ordered_dend_neurons = [t for sublist in list(ordered_dend_raw) for t in sublist]

ordered_id_dend_list = [neuron_id[0] for neuron_id in ordered_dend_raw]
id_dend_map = {i : j for i,j in zip(ordered_id_dend_list, range(len(ordered_id_dend_list)))}


#print id_soma_map
#for i in xrange(len(ordered_soma_spikes)):
#    print ordered_soma_spikes[i], ordered_neurons[i]

f = plt.figure()
ax1 = f.add_subplot(211)

for i in range(len(spike_times_dend)):
    ax1.vlines(spike_times_dend[i], neuron_fired_dend[i]-0.5, neuron_fired_dend[i]+0.5)

#ax.scatter(spike_times, random_ID)
#ax1.set_yticks(random_ID)

#plt.yticks(random_ID, neuron_fired)
ax1.set_xlim([-5, TRIAL_DURATION])
ax1.set_ylabel("neuron ID")
#ax1.set_xlabel("Time (ms)")
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
ax2.set_xlabel("Time (ms)")
ax2.set_title("Somatic spikes")

f = plt.figure()
ax1 = f.add_subplot(211)

min_soma_spike_time = ordered_soma_spikes[0]
max_soma_spike_time = ordered_soma_spikes[-1] - ordered_soma_spikes[0]

for t, i in zip(ordered_soma_spikes, ordered_soma_neurons):
    if i == 513 or i == 809 or i == 882:
        print i,t
    ax1.vlines(t - min_soma_spike_time, id_soma_map[i]-0.5, id_soma_map[i]+0.5)

#ax.scatter(spike_times, random_ID)

ax1.set_yticks(range(len(ordered_id_soma_list)))

#plt.yticks(range(len(ordered_id_soma_list)), ordered_id_soma_list)
plt.yticks([])

ax1.set_xlim([-5, max_soma_spike_time + 50])
ax1.set_ylabel("neuron ID")
#ax1.set_xlabel("Time (ms)")
ax1.set_title("Ordered somatic spikes")
ax1.set_ylim([-5, len(ordered_id_soma_list)+5])

ax2 = f.add_subplot(212)

min_dend_spike_time = ordered_dend_spikes[0]
max_dend_spike_time = ordered_dend_spikes[-1] - ordered_dend_spikes[0]

for t, i in zip(ordered_dend_spikes, ordered_dend_neurons):
    ax2.vlines(t - min_dend_spike_time, id_dend_map[i]-0.5, id_dend_map[i]+0.5)

#ax.scatter(spike_times, random_ID)

ax2.set_yticks(range(len(ordered_id_dend_list)))

#plt.yticks(range(len(ordered_id_dend_list)), ordered_id_dend_list)
plt.yticks([])

ax2.set_xlim([-5, max_dend_spike_time + 50])
ax2.set_ylabel("neuron ID")
ax2.set_xlabel("Time (ms)")
ax2.set_title("Ordered dendritic spikes")
ax2.set_ylim([-5, len(ordered_id_dend_list)+5])


plt.show()