# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 16:49:22 2016

@author: jingroup
"""

"""
Script plots average spike times of neurons in the chain
"""

import reading
import numpy as np
import matplotlib.pyplot as plt
import math

file_spikes = "/home/eugene/lionX/noNMDA_6/chain_test.bin"
RARA = "/home/eugene/lionX/noNMDA_6/super163.net"

SUCCESS_PERCENT = 0.20

def std(x):
    if len(x) == 0:
        return -1
        
    mean = np.mean(x)
    sum = 0
    for e in x:
        sum += (e-mean)**2
    return math.sqrt(sum/len(x))

# read spikes from file
N_RA, num_trials, num_dend_spikes, mean_burst_time_raw, std_burst_time_raw = reading.read_chain_test(file_spikes)


# get rid of all -1 elements

mean_burst_time = []
std_burst_time = []
neuron_id = []

for i in xrange(N_RA):
    if num_dend_spikes[i] > SUCCESS_PERCENT * num_trials:
        neuron_id.append(i)        
        mean_burst_time.append(mean_burst_time_raw[i])
        std_burst_time.append(std_burst_time_raw[i])
        
mean_burst_time, neuron_id, std_burst_time = \
                zip(*sorted(zip(mean_burst_time, neuron_id, std_burst_time)))


print "neuron_id:", neuron_id    
print "mean_burst_time: ", mean_burst_time
print "std_burst_time: ", std_burst_time

f1 = plt.figure()
ax1 = f1.add_subplot(111)

for i in range(len(mean_burst_time)):
    sigma = std_burst_time[i]
    mean = mean_burst_time[i]
        
    ax1.vlines(mean, i-0.5, i+0.5, lw=2)
    ax1.hlines(i, mean - sigma, mean + sigma)
    ax1.vlines(mean-sigma, i-0.25, i+0.25)
    ax1.vlines(mean+sigma, i-0.25, i+0.25)
        

#ax.scatter(spike_times, random_ID)
#ax1.set_yticks(random_ID)

plt.yticks(xrange(len(mean_burst_time)), neuron_id)
ax1.set_xlim([-5, mean_burst_time[-1] + max(std_burst_time)])
ax1.set_ylabel("neuron id")
ax1.set_xlabel("burst time (ms)")
#ax1.set_title("Statistics of spike times")
ax1.grid(True)

f2 = plt.figure()
ax2 = f2.add_subplot(111)

ax2.plot(range(0, len(mean_burst_time)), std_burst_time)
ax2.set_ylabel("standard deviation (ms)")
ax2.set_xlabel("neuron ID")

plt.show()