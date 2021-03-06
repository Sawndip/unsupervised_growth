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

from utils import *



N = 2000 # number of neurons

#fileDend = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays6/spike_times_dend_1900.bin"
#fileSoma = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays6/spike_times_soma_1900.bin"

#fileDend = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/matTrans3_network2000RA550I/spike_times_dend_11550.bin"
#fileSoma = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/matTrans3_network2000RA550I/spike_times_soma_11550.bin"

#fileDend = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/split1/spike_times_dend_7200.bin"
#fileSoma = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/split1/spike_times_soma_7200.bin"


#fileDend = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/chain_1.0/test_spike_times_dend_0.bin"
#fileSoma = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/chain_1.0/test_spike_times_soma_0.bin"


fileDend = "/mnt/hodgkin/eugene/results/immature/clusters/matTrans54/spike_times_dend_6800.bin"
fileSoma = "/mnt/hodgkin/eugene/results/immature/clusters/matTrans54/spike_times_soma_6800.bin"

#fileDend = "/home/eugene/results/immature/clusters/test/matTrans44/test_spike_times_dend_10.bin"
#fileSoma = "/home/eugene/results/immature/clusters/test/matTrans44/test_spike_times_soma_10.bin"

#fileDend = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/test/maturationTransition4/test_spike_times_dend_7.bin"
#fileSoma = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/test/maturationTransition4/test_spike_times_soma_7.bin"


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


ordered_dend_spikes, ordered_dend_neurons = zip(*sorted(zip(spike_times_dend, neuron_fired_dend)))


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


bin_width = 1.0

ordered_dend_spikes = [d - min_dend_spike_time for d in ordered_dend_spikes]

time, burst_density = calculate_burst_density(ordered_dend_spikes, bin_width)

start_time = 0.0
end_time = 120.0
        

print "Mean burst density: ",np.mean(burst_density[(time > start_time) & (time < end_time)])
print "Std burst density: ",np.std(burst_density[(time > start_time) & (time < end_time)])
print "std / mean = ",np.std(burst_density[(time > start_time) & (time < end_time)]) / np.mean(burst_density[(time > start_time) & (time < end_time)])

plt.figure()
plt.plot(time, burst_density)
plt.xlabel('Time (ms)')
plt.ylabel('# of bursts / ms')
#plt.title('Burst density in causal network $G_Emax = 1.5 G_L; G_Imax = 0.5 G_L$')
#plt.title('pruned spread 0.0 ms')
plt.xlim([0,300])
plt.ylim([0,10])


plt.show()