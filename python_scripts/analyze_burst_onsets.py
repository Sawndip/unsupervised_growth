#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 12:40:25 2017

@author: jingroup

Script analyzes burst onsets in the network
"""

import reading
import matplotlib.pyplot as plt
import numpy as np

def calculate_burst_density(burst_times, bin_width):
    """
    Calculate burst density - number of bursts in time bins
    """
    size = int(burst_times[-1] / bin_width) + 1
    time = np.array([float(i)*bin_width + bin_width/2. for i in range(size)])
    
    num_bursts = np.zeros(size, np.int32)    
    
    for burst_time in burst_times:
        num_bursts[int(burst_time / bin_width)] += 1
        
    burst_density = num_bursts / bin_width
    
    return time, burst_density

CURRENT_INJECTION = 50.0

#filename = "/home/eugene/Output/networks/MLong/simResults/dendritic_tree/causal_1connection/ee0.06_ie0.03_out50active30_1ms_somaSpikes.bin"
filename = "/home/eugene/Output/networks/MLong/simResults/dendritic_tree/random_5500hvcI/ee0.075_ie0.05_burstOnsets_fullTrial_1_somaSpikes.bin"

#filename = "/home/eugene/Output/networks/MLong/simResults/2chain/e0.50_i0.1_rerouted_0ms_dendSpikes.bin"

(trial_number, simulation_time, spike_times_raw, neuron_fired) = reading.read_time_info(filename)

LEAK_CONDUCTANCE = 0.1

N_RA = 20000

print "Number of HVC(RA) neurons = ",N_RA

#print neuron_fired

threshold = 1.0

all_first_spike_times = np.empty(N_RA, np.float32)
all_first_spike_times.fill(-100.0)

num_bursts = 0

for n, time in zip(neuron_fired, spike_times_raw):
    num_bursts += 1
    all_first_spike_times[n[0]] = time[0] - CURRENT_INJECTION

print list(all_first_spike_times[:200])

id_last_fired = np.argmax(all_first_spike_times)

print "Total number of bursts = ",num_bursts

print "Number of silent neurons: ",np.shape(np.where(all_first_spike_times[:(id_last_fired+1)] < 0)[0])[0]
print "id of last fired HVC(RA) = ",id_last_fired
print "Max burst time relative to current injection: ",np.max(all_first_spike_times)

#print burst_times

inh = 0.05 / LEAK_CONDUCTANCE
exc = 0.15 / LEAK_CONDUCTANCE

num_neurons_to_plot = 20000
time_to_plot = 200

#==============================================================================
# # plot unsorted burst times
# f = plt.figure()
# ax1 = f.add_subplot(111)
# 
# 
# print "Max burst time relative to current injection: ",np.max(all_burst_times)
# print "id of neuron that bursts latest: ", np.argmax(all_burst_times)
# num_bursts = []
# 
# 
# 
# #==============================================================================
# # for i in range(len(neuron_fired)):
# #     if neuron_fired[i][0] < num_neurons_to_plot:
# #         num_bursts.append(len(burst_times_raw[i]))
# #         if num_bursts[-1] > threshold:
# #         
# #         for burst_time in burst_times_raw[i]:
# #             ax1.vlines(burst_time - CURRENT_INJECTION, neuron_fired[i][0]-0.5, neuron_fired[i][0]+0.5)
# #     else:
# #         break
# # 
# #==============================================================================
# for i in range(num_neurons_to_plot):
#     if all_burst_times[i] > 0 and all_burst_times[i] < time_to_plot:     
#         ax1.vlines(all_burst_times[i], i-0.5, i+0.5)
#     #else:
#      #   ax1.vlines(all_burst_times[i], i-0.5, i+0.5)
# 
# 
# 
# plt.tick_params(axis='y',which='both',bottom='off',top='off',labelbottom='off')
#     
# ax1.set_ylim([-50, num_neurons_to_plot + 50])      
# ax1.set_xlim([0, time_to_plot])
# ax1.set_ylabel("neuron id")
# ax1.set_xlabel("Time (ms)")
# ax1.set_title("Burst times $G_Imax$ = {0} $G_L$; $G_Emax$ = {1} $G_L$".format(inh, exc))
#==============================================================================


# plot sorted burst times
burst_times_of_first_neurons = all_first_spike_times[:num_neurons_to_plot]

min_burst_time = np.min(burst_times_of_first_neurons[burst_times_of_first_neurons > -10])
burst_times_of_first_neurons = burst_times_of_first_neurons - min_burst_time
burst_times_sorted = np.sort(burst_times_of_first_neurons[burst_times_of_first_neurons > -10])


print burst_times_of_first_neurons[burst_times_of_first_neurons > 42.0 + CURRENT_INJECTION]
print burst_times_sorted[burst_times_sorted > 42.0 + CURRENT_INJECTION]
#spike_times = np.array([spike_time for sublist in spike_times_raw for spike_time in sublist])
#neuron_fired = np.array([neuron_id for sublist in neuron_fired for neuron_id in sublist])

#spike_times = spike_times[neuron_fired < num_neurons_to_plot]

#spike_times = np.sort(spike_times)



f = plt.figure()
ax1 = f.add_subplot(111)

for i in range(burst_times_sorted.shape[0]):
    if burst_times_sorted[i] < time_to_plot:
        ax1.vlines(burst_times_sorted[i], i-0.5, i+0.5)

plt.tick_params(axis='y',which='both',bottom='off',top='off',labelbottom='off')
    
#ax1.set_ylim([-50, num_neurons_to_plot + 50])    
ax1.set_ylim([-50, 2900])   

ax1.set_xlim([0, 150])
#ax1.set_xlim([0, time_to_plot])

ax1.set_ylabel("Neuron id")
ax1.set_xlabel("Time (ms)")
#ax1.set_title("Sorted bursts $G_Imax$ = {0} $G_L$; $G_Emax$ = {1} $G_L$".format(inh, exc))
#ax1.set_title('pruned spread 0.0 ms')

bin_width = 0.75

time, burst_density = calculate_burst_density(burst_times_sorted, bin_width)


plt.figure()
plt.plot(time, burst_density)
plt.xlabel('Time (ms)')
plt.ylabel('# of bursts / ms')
#plt.title('Burst density in causal network $G_Emax = 1.5 G_L; G_Imax = 0.5 G_L$')
#plt.title('pruned spread 0.0 ms')
plt.xlim([0,150])
plt.ylim([0,160])


plt.show()    

