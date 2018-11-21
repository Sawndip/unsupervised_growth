# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 11:36:44 2018

@author: jingroup

Script detects peaks in data
"""

import reading
import matplotlib.pyplot as plt
import numpy as np

import peakutils
from peakutils.plot import plot as pplot

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

#filename = "/home/eugene/Output/networks/MLong/simResults/dendritic_tree/causal_5500hvcI/out50_active360_10ms_5500hvcI/e0.060000_i0.050000_somaSpikes.bin"
#filename = "/home/eugene/Output/networks/MLong/simResults/dendritic_tree/causal_5500hvcI/out50_active210_7ms_5500hvcI/e0.060000_i0.050000_somaSpikes.bin"
#filename = "/home/eugene/Output/networks/MLong/simResults/dendritic_tree/perfectChainRandom/ee0.015_ie0.03_newModel_burstOnsets_0_somaSpikes.bin"
#filename = "/home/eugene/Output/networks/MLong/simResults/dendritic_tree/causal_5500hvcI/out50_active150_5ms_5500hvcI_entire/ee0.060_ie0.05_burstOnsets_0_somaSpikes.bin"
#filename = '/home/eugene/Output/networks/MLong/simResults/dendritic_tree/causal_5500hvcI/out50_sampleWithRewiring_5500hvcI/e0.060000_i0.050000_somaSpikes.bin'
#filename = '/home/eugene/Output/networks/MLong/simResults/dendritic_tree/causal_5500hvcI/out50_active30_1ms_pseudo_growth_synchronous/e0.055000_i0.050000_somaSpikes.bin'
#filename = "/home/eugene/Output/networks/MLong/simResults/dendritic_tree/causal_5500hvcI/out50_active30_1ms_pseudo_growth_bursts/ee0.055_ie0.05_burstOnsets_0_somaSpikes.bin"
#filename = "/mnt/hodgkin/eugene/Output/networks/MLong/simResults/dendritic_tree/parallelChains_5500hvcI/ee0.4_ie0.05_uniform_burstOnsets_0_somaSpikes.bin"
#filename = "/mnt/hodgkin/eugene/results/MLong/prunedChain/ee0.015_ie0.05_uniform_burstOnsets_10_somaSpikes.bin"
#filename = "/mnt/hodgkin/eugene/results/MLong/polychronousChain/ee0.040_ie0.05_uniform_burstOnsets_10_somaSpikes.bin"
  
#filename = "/mnt/hodgkin/eugene/Output/networks/MLong/simResults/dendritic_tree/causal_5500hvcI/out50_active30_1ms_e0.050_5500hvcI/ee0.040_ie0.05_burstOnsets_2_somaSpikes.bin"
#filename = "/mnt/hodgkin/eugene/Output/networks/MLong/simResults/dendritic_tree/random_5500hvcI/ee0.075_ie0.05_burstOnsets_2_somaSpikes.bin"

#filename = "/home/eugene/Output/networks/MLong/simResults/2chain/e0.50_i0.1_rerouted_0ms_dendSpikes.bin"

filename = "/mnt/hodgkin/eugene/results/MLong/parallelChains/ee0.100_ie0.05_2chains_0_somaSpikes.bin"


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

#print list(all_first_spike_times[:200])


min_burst_time = np.min(all_first_spike_times[all_first_spike_times > -10])
all_first_spike_times = all_first_spike_times - min_burst_time
burst_times_sorted = np.sort(all_first_spike_times)


bin_width = 1.0

time, burst_density = calculate_burst_density(burst_times_sorted, bin_width)

indexes = peakutils.indexes(burst_density, thres=0.0, min_dist=3)
indexes_negative = peakutils.indexes(-burst_density, thres=0.0, min_dist=3)
print(indexes)
print(time[indexes], burst_density[indexes])
plt.figure(figsize=(10,6))
#pplot(time, burst_density, indexes)
pplot(time, burst_density, indexes_negative)

plt.title('First estimate')

plt.figure()
plt.plot(time, burst_density)
plt.xlabel('Time (ms)')
plt.ylabel('# of bursts / ms')
#plt.title('Burst density in causal network $G_Emax = 1.5 G_L; G_Imax = 0.5 G_L$')
#plt.title('pruned spread 0.0 ms')
plt.xlim([0,150])
plt.ylim([0,160])


plt.show()    