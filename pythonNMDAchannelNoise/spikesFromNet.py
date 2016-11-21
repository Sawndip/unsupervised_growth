# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 18:09:57 2016

@author: jingroup
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import reading
import matplotlib.pyplot as plt
import numpy as np
import math

TRIAL_DURATION = 150

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
    


filename = "/home/eugene/Output/soma_spikes.bin"
RARA = "/home/eugene/Output/super163.net"

(trial_number, simulation_time, spike_times, neuron_fired) = reading.read_time_info(filename)

#spike_times = [x - min(spike_times) for x in spike_times]

(N_RA, stronglyConnected) = reading.read_connections_net(RARA)

print "stronglyConnected:", stronglyConnected


spike_times = [t for sublist in list(spike_times) for t in sublist]
neuron_fired = [ind for sublist in list(neuron_fired) for ind in sublist]

first_spike = min(spike_times)

spike_times = map(lambda x: x - first_spike, spike_times)

print "spike_times: ", spike_times
print "neuron_fired: ", neuron_fired

stronglyConnectedFired = [i for i in neuron_fired if i in stronglyConnected]
stronglyConnectedSpikeTimes = [t for ind, t in enumerate(spike_times) if neuron_fired[ind] in stronglyConnected]

print "stronglyConnectedFired: ",stronglyConnectedFired
print "stronglyConnectedSpikeTimes: ",stronglyConnectedSpikeTimes

first_spike_stronglyConnected = min(stronglyConnectedSpikeTimes)

stronglyConnectedSpikeTimes = map(lambda x: x - first_spike_stronglyConnected, stronglyConnectedSpikeTimes)
#==============================================================================
# spikes_to_erase = []
# mean_spike_time = np.mean(spike_times)
# 
# std_spike_times = std(spike_times)
# 
# print "mean spike time: ", mean_spike_time
# print "standard deviation of spike times: ", std_spike_times
# 
# for ind, t in enumerate(spike_times):
#     if np.fabs(t - mean_spike_time) > 5 * std_spike_times:
#         spikes_to_erase.append(ind)
# 
# delete_counter = 0
# 
# for i in spikes_to_erase:
#     del spike_times[ind - delete_counter]
#     del neuron_fired[ind - delete_counter]
#     delete_counter += 1
# 
# print "spike_times after deleting: ", spike_times
# print "neuron_fired after deleting: ", neuron_fired
#==============================================================================

spike_times, neuron_fired = zip(*sorted(zip(spike_times, neuron_fired)))
stronglyConnectedSpikeTimes, stronglyConnectedFired = zip(*sorted(zip(stronglyConnectedSpikeTimes, stronglyConnectedFired)))

print spike_times

print "sorted spike_times: ", spike_times
print "sorted neuron_fired: ", neuron_fired

#random_ID = range(len(neuron_fired))

random_ID = remap_ID(neuron_fired)
random_ID_stronglyConnected = remap_ID(stronglyConnectedFired)

print "random IDs assigned to neurons: ", random_ID

print "Number of strongly connected fired neurons: ",max(random_ID_stronglyConnected) + 1

f1 = plt.figure()
ax1 = f1.add_subplot(211)

for i in range(len(spike_times)):
    ax1.vlines(spike_times[i], random_ID[i]-0.5, random_ID[i]+0.5)

#ax.scatter(spike_times, random_ID)
ax1.set_yticks(random_ID)

plt.yticks(random_ID, neuron_fired)
ax1.set_xlim([-5, max(spike_times)+5])
ax1.set_ylabel("real neuron ID")
ax1.set_xlabel("relative spike time (ms)")
ax1.set_title("Ordered spikes")

ax2 = f1.add_subplot(212)

for i in range(len(stronglyConnectedSpikeTimes)):
    ax2.vlines(stronglyConnectedSpikeTimes[i], random_ID_stronglyConnected[i]-0.5, random_ID_stronglyConnected[i]+0.5)

#ax.scatter(spike_times, random_ID)
ax2.set_yticks(random_ID_stronglyConnected)

plt.yticks(random_ID_stronglyConnected, stronglyConnectedFired)
ax2.set_xlim([-5, max(stronglyConnectedSpikeTimes)+5])
ax2.set_ylabel("real neuron ID")
ax2.set_xlabel("relative spike time (ms)")
ax2.set_title("Ordered spikes of strongly connected neurons")

f2 = plt.figure()
ax3 = f2.add_subplot(111)

for i in range(len(stronglyConnectedSpikeTimes)):
    ax3.vlines(stronglyConnectedSpikeTimes[i], random_ID_stronglyConnected[i]-0.5, random_ID_stronglyConnected[i]+0.5)

#ax.scatter(spike_times, random_ID)
ax3.set_yticks(random_ID_stronglyConnected)
ax3.set_yticklabels([])
#plt.yticks(random_ID_stronglyConnected, stronglyConnectedFired)
ax3.set_ylim([-1, 60])
ax3.set_xlim([-5, TRIAL_DURATION])
ax3.set_xlabel("spike time (ms)")
#ax3.set_title("Ordered spikes of strongly connected neurons")


plt.show()    
