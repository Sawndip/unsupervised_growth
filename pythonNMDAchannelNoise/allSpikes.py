# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import reading
import matplotlib.pyplot as plt
import numpy as np
import math

def std(x):
    mean = np.mean(x)
    sum = 0
    for e in x:
        sum += (e-mean)**2
    return math.sqrt(sum/len(x))

def remap_ID(order):
    already_counted = []
    remapped = []
    ind = 0
    for i in order:
        if i in already_counted:
            ind_in_counted = already_counted.index(i)
            remapped.append(ind_in_counted)
        else:
            remapped.append(ind)
            already_counted.append(i)
            ind += 1
        
    return remapped
    


filename = "/home/eugene/Output/networks/IdealChainTest220217/soma_spikes_trial1.bin"

(trial_number, simulation_time, spike_times, neuron_fired) = reading.read_time_info(filename)

spike_times = [t for sublist in list(spike_times) for t in sublist]
neuron_fired = [ind for sublist in list(neuron_fired) for ind in sublist]


spike_times, neuron_fired = zip(*sorted(zip(spike_times, neuron_fired)))

random_ID = remap_ID(neuron_fired)


f = plt.figure()
ax1 = f.add_subplot(111)

for i in range(len(spike_times)):
    ax1.vlines(spike_times[i], random_ID[i]-0.5, random_ID[i]+0.5)

#ax.scatter(spike_times, random_ID)
ax1.set_yticks(random_ID)

plt.yticks(random_ID, neuron_fired)
ax1.set_xlim([-5, max(spike_times)+5])
ax1.set_ylabel("real neuron ID")
ax1.set_xlabel("relative spike time (ms)")
ax1.set_title("Ordered spikes")



plt.show()    
