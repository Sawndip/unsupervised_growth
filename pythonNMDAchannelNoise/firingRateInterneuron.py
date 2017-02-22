# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 09:28:35 2017

@author: jingroup

Script analyzes HVC(I) firing rate aligned to dendritic burst of HVC(RA) neuron
"""
import reading
import matplotlib.pyplot as plt
import numpy as np
import os

WINDOW_SIZE = 100 # ms
TIMESTEP = 0.02 # timestep of dynamics in ms
CURRENT_INJECTION = 100 # time of current injection to training neurons in ms
RATE_WINDOW_SIZE = 5.0 # bin in ms in which to calculate neuron firing rate
DENDRITIC_THRESHOLD = 0.0 # threshold for dendritic spike in mV


filename = "/home/eugene/hodgkinData/gabaMaturation300117/membraneTraces/I/I4_trial2.bin"

dataDirRA = "/home/eugene/hodgkinData/gabaMaturation300117/membraneTraces/RA"
dataDirI = "/home/eugene/hodgkinData/gabaMaturation300117/membraneTraces/I"


(t, v, _, _, _, _, _, _, _, flag, _) = reading.read_hhi(filename)

def get_burst_times(filename):
    """
    Extracts dendritic burst indices from file with HVC(RA) data
    """
    (t, _, _, _, _, Vd, _, _, _, _, _, _, _, _, _, _, _, _) = reading.read_hh2(filename)
    
    burst_indices = []
    burst_times = []
    flag = False
        
    for i, v in enumerate(Vd):
        if flag == False and v >= DENDRITIC_THRESHOLD:
            flag = True
        elif flag == True and v < DENDRITIC_THRESHOLD:
            burst_indices.append(i)
            burst_times.append(t[i])
            flag = False
        
    return burst_times, burst_indices

def get_spike_times(filename):
    """
    Extract neuronal spike times
    """    
    (t, _, _, _, _, _, _, _, _, flag, _) = reading.read_hhi(filename)

    flag_previous = flag[0]    

    spike_times = []
    
    for i in range(1, len(flag)):
        if flag[i] - flag_previous < -0.5:
            spike_times.append(t[i])
        flag_previous = flag[i]
        

    return t, spike_times


def firing_rate_sliding_window(t, spike_times):
    """
    Calculate firing rate from spike times using sliding window
    """
    firing_rate = np.zeros(len(t), np.float32)
    #print t
    #print RATE_WINDOW_SIZE / 2.0
    for i in range(len(t)):
        num_spikes_in_window = 0
        for st in spike_times:
            if round(st, 2) >= round(t[i], 2) - RATE_WINDOW_SIZE / 2.0 and round(st, 2) <= round(t[i], 2)+ RATE_WINDOW_SIZE / 2.0:
                #print "t = ",t[i]
                #print "st = ",st
                num_spikes_in_window += 1
                
        if num_spikes_in_window > 0:
            firing_rate[i] = 1000 * num_spikes_in_window / RATE_WINDOW_SIZE

    return firing_rate

def get_spike_triggered_firing_rate(trial_num, fileI, dataDirRA):
    """
    Calculate average firing rate of inhibitory neuron aligned to dendritic bursts
    of HVC(RA) neurons in the chain
    """
    spike_triggered_firing_rate = np.zeros(int(WINDOW_SIZE / TIMESTEP) + 1, np.float32)
    t_spike_triggered = np.zeros(int(WINDOW_SIZE / TIMESTEP) + 1, np.float32)
    
    t, spike_times = get_spike_times(fileI)
    
    # if HVC(I) neuron didn't spike simply return array of zeros
    if len(spike_times) == 0:
        print "Neuron with file {0} did NOT file!".format(fileI)
        return t_spike_triggered, spike_triggered_firing_rate  
    else:
        firing_rate = firing_rate_sliding_window(t, spike_times)
        
        filenames = os.listdir(dataDirRA)    
        
        
        
        num_fired_bursts = 0 # number of HVC(RA) neurons that fired dendritic burst
        
        for n in filenames:                
            filename = os.path.join(dataDirRA, n)
            if "trial" + str(trial_num) + ".bin" in filename: # check only HVC(RA) neurons from the same trial
                #print filename
                burst_times, burst_indices = get_burst_times(filename) # burts time of HVC(RA) neuron
                
                #print burst_times            
    
                if len(burst_times) > 0:
                    num_fired_bursts += 1
                
                for burst_index in burst_indices:
                    window_size_in_indices = int(WINDOW_SIZE / (2*TIMESTEP))
                    
                    #print "burst_index = ",burst_index
                    #print "window_size_in_indices = ",window_size_in_indices
                    
                    for i in range(window_size_in_indices):
                        t_spike_triggered[i] = - WINDOW_SIZE / 2 + i * TIMESTEP  
                        t_spike_triggered[i + window_size_in_indices + 1] = (i + 1) * TIMESTEP
                        spike_triggered_firing_rate[i] += firing_rate[burst_index - window_size_in_indices + i]
                        spike_triggered_firing_rate[i + window_size_in_indices + 1] += firing_rate[burst_index + i + 1]
                        
                    spike_triggered_firing_rate[window_size_in_indices] += firing_rate[burst_index]
                    t_spike_triggered[window_size_in_indices] = 0.0
        
        print num_fired_bursts    
        
        if num_fired_bursts > 0:
            spike_triggered_firing_rate /= float(num_fired_bursts)
        
        return t_spike_triggered, spike_triggered_firing_rate
                
def get_average_spike_triggered_firing_rate(trial_num, dataDirI, dataDirRA):    
    """
    Calculate average spike triggered firing rate for all interneurons
    """
    filenames = os.listdir(dataDirI)
    
    average_spike_triggered_firing_rate = np.zeros(int(WINDOW_SIZE / TIMESTEP) + 1, np.float32)
    
    # loop through all interneurons:
    for n in filenames:                
        filename = os.path.join(dataDirI, n)
        if "trial" + str(trial_num) + ".bin" in filename: # check only HVC(I) neurons from the same trial
            print filename            
            t_spike_triggered, spike_triggered_firing_rate = get_spike_triggered_firing_rate(trial_num, filename, dataDirRA)
            
            average_spike_triggered_firing_rate += spike_triggered_firing_rate
            
    if len(filenames) > 0:
        average_spike_triggered_firing_rate /= len(filenames)
        
    return t_spike_triggered, average_spike_triggered_firing_rate
            
t, spike_times = get_spike_times(filename)
print spike_times
firing_rate = firing_rate_sliding_window(t, spike_times)   
 
 
print firing_rate    

trial_num = 1


t_spike_triggered, average_spike_triggered_firing_rate = get_average_spike_triggered_firing_rate(trial_num, dataDirI, dataDirRA)
    
fig1 = plt.figure()

ax = fig1.add_subplot(111)
ax.plot(t_spike_triggered, average_spike_triggered_firing_rate, 'g', linewidth = 2.0)
ax.set_ylabel("spike_triggered_firing_rate")
ax.set_xlabel("t (ms)")

plt.grid(True)

plt.show()
