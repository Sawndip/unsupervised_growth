#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 21 14:05:29 2018

@author: jingroup

Script calculates average burst density of multiple test runs
and estimates oscillations in burst density pattern
"""

import numpy as np
import reading
import os
import matplotlib.pyplot as plt

N_RA = 2000
CURRENT_INJECTION = 100.0
START_TIME = -100.0
END_TIME = 600.0

def get_burst_density_for_run(bin_width, filename):
    """
    Calculates burst density for a given simulation run
    """
    (trial_number, simulation_time, spike_times_raw, neuron_fired) = reading.read_time_info(filename)
    
    print "file = ",filename    
        
    all_first_spike_times = np.empty(N_RA, np.float32)
    all_first_spike_times.fill(-100.0)
    
    num_bursts = 0
    
    for n, time in zip(neuron_fired, spike_times_raw):
        num_bursts += 1
        all_first_spike_times[n[0]] = time[0]
    
    #print list(all_first_spike_times[:200])
    
    id_last_fired = np.argmax(all_first_spike_times)
    
    print "Total number of bursts = ",num_bursts
    
    print "Number of silent neurons: ",np.shape(np.where(all_first_spike_times[:(id_last_fired+1)] < 0)[0])[0]
    print "id of last fired HVC(RA) = ",id_last_fired
    print "Max burst time relative to current injection: ",np.max(all_first_spike_times)

    burst_times_of_fired = all_first_spike_times[all_first_spike_times > 0] 

    burst_times_sorted = np.sort(burst_times_of_fired - CURRENT_INJECTION)
    #print burst_times_sorted
    
    
    
    time, burst_density = calculate_burst_density(burst_times_sorted, bin_width)
    return time, burst_density

def calculate_burst_density(burst_times, bin_width):
    """
    Calculate burst density - number of bursts in time bins
    """
    size = int((END_TIME - START_TIME) / bin_width) + 1
    time = np.array([float(i)*bin_width + START_TIME + bin_width/2. for i in range(size)])
    
    num_bursts = np.zeros(size, np.int32)    
    
    #print min(burst_times)
    
    burst_times = burst_times - START_TIME
    
    #print min(burst_times)
    
    for burst_time in burst_times:
        num_bursts[int(burst_time / bin_width)] += 1
        
    burst_density = num_bursts / bin_width
    
    return time, burst_density

def write_to_file(time, burst_density, filename):
    """
    Write burst density and time to a file
    """
    with open(filename, 'w') as f:
        for t, b in zip(time, burst_density):
            f.write(str(t))
            f.write("\t")
            f.write(str(b))
            f.write("\n")
            

if __name__ == "__main__":
    bin_width = 1.0

    dirname = "/home/eugene/results/immature/clusters/test/matTrans44/"
        
        
    
    files = os.listdir(dirname)
    
    av_burst_density = None
    time = None
    num_runs = 0
    
    for f in files:
        if "soma" in f:
            num_runs += 1
            print os.path.join(dirname, f)
            t, burst_density = get_burst_density_for_run(bin_width, os.path.join(dirname, f))
            
            if time == None:
                time = t
            
            if av_burst_density == None:
                av_burst_density = burst_density
            else:
                av_burst_density = av_burst_density + burst_density
                
            #print av_burst_density[:200]
    
    if num_runs > 0:
        #print list(time)
        
        
        av_burst_density = av_burst_density / float(num_runs)

       # write_to_file(time, av_burst_density, os.path.join(dirname,filename))
        #write_to_file(time, av_burst_density, fileOut)

        #first_nonzero_bin = np.where(av_burst_density > 0)[0][0]        
        #last_nonzero_bin = np.where(av_burst_density > 0)[0][-1]

        start_time = 0.0
        end_time = 120.0

        #print "Mean burst density: ",np.mean(av_burst_density[first_nonzero_bin:last_nonzero_bin])
        #print "Std burst density: ",np.std(av_burst_density[first_nonzero_bin:last_nonzero_bin])
        #print "Std / mean = ", np.std(av_burst_density[first_nonzero_bin:last_nonzero_bin]) / np.mean(av_burst_density[first_nonzero_bin:last_nonzero_bin])
        
        
        print "Mean burst density: ",np.mean(av_burst_density[(time > start_time) & (time < end_time)])
        print "Std burst density: ",np.std(av_burst_density[(time > start_time) & (time < end_time)])
        
        print "Std / mean = ", np.std(av_burst_density[(time > start_time) & (time < end_time)]) / np.mean(av_burst_density[(time > start_time) & (time < end_time)])
        
        
        #indexes = peakutils.indexes(av_burst_density, thres=0.0, min_dist=10)
        #indexes_negative = peakutils.indexes(-av_burst_density, thres=0.0, min_dist=10)
        #print(indexes)
        #print(time[indexes], av_burst_density[indexes])
        #plt.figure(figsize=(10,6))
        #pplot(time, av_burst_density, indexes)
        #pplot(time, av_burst_density, indexes_negative)        
        
        plt.figure()
        plt.plot(time, av_burst_density)
        plt.xlabel('Time (ms)')
        plt.ylabel('# of bursts / ms')
        #plt.title('Burst density in causal network $G_Emax = 1.5 G_L; G_Imax = 0.5 G_L$')
        #plt.title('pruned spread 0.0 ms')
        plt.xlim([-10,200])
        plt.ylim([0,10])
        
        
        plt.show()    
        
#filename = "/home/eugene/Output/networks/MLong/simResults/2chain/e0.50_i0.1_rerouted_0ms_dendSpikes.bin"









