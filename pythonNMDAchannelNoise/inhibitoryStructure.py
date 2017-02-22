# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 16:06:35 2017

@author: jingroup

Script calculates how inhibitory inputs are correlated the burst time
"""
import reading
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

dataDir = "/home/eugene/Output/networks/IdealChainTest220217/RA/"

DENDRITIC_THRESHOLD = 0.0 # threshold for burst spike in mV
INHIBITORY_KICK_THRESHOLD = 0.1 # threshold for inhibitory inputs
INHIBITORY_SPIKE_RESOLUTION = 1.0 # resolution in ms wihin which we treat inhibitory spike as the same

def get_burst_times(t, Vd):
    """
    Extracts dendritic burst times from voltage data
    """
    burst_times = []
    flag = False
        
    for i, v in enumerate(Vd):
        if flag == False and v >= DENDRITIC_THRESHOLD:
            flag = True
        elif flag == True and v < DENDRITIC_THRESHOLD:
            burst_times.append(t[i])
            flag = False
        
    return burst_times

def get_inhibitory_input_times(t, Ginh_d):
    """
    Extracts times of inhibitory inputs
    """
    inhibitory_kick_times = []
    inhibitory_kick_strength = [] 
    
    for i, g in enumerate(Ginh_d):
        if i == 0:
            g_previous = g
        else:
            if g - g_previous >= INHIBITORY_KICK_THRESHOLD:
                inhibitory_kick_times.append(t[i])
                inhibitory_kick_strength.append(g - g_previous)
        
            g_previous = g           

    return (inhibitory_kick_times, inhibitory_kick_strength)

def inh_input_burst_difference(t, Vd, Ginh_d):
    """
    Calculates difference between inhibitory input times and burst time
    """
    burst_times = get_burst_times(t, Vd)
    inhibitory_kick_times, inhibitory_kick_strength = get_inhibitory_input_times(t, Ginh_d)
    
    dt = []    
    
    for ti in inhibitory_kick_times:
        for tb in burst_times:
            dt.append(ti - tb)
    
    return dt

def average_inh_input_burst_difference_one_neuron(neuron_id, num_trials):
    """
    Computes average time difference between inhibitory inputs and dendritic bursts 
    for the same neuron over num_trials trials
    """    
    dt_all_trials = []

    num_successful_trials = 0 # number of trials when neuron fired
    # store time differences for all trials in one array 
    for i in range(1, num_trials + 1):
        filename = dataDir + "RA" + str(neuron_id) + "_trial" + str(i) + ".bin"
        (t, Vs, _, _, _, Vd, _, _, _, _, Gexc_d, Ginh_d, _, _, _, _, _, _) = reading.read_hh2(filename)
        
        dt = inh_input_burst_difference(t, Vd, Ginh_d)
        
        if len(dt) > 0:     
            dt_all_trials.extend(dt)
            num_successful_trials += 1
        
    dt_all_trials = sorted(dt_all_trials)
    
    #print dt_all_trials
    
   
    if len(dt_all_trials) == 0:
        print("dt_all_trials is empty for neuron {0}".format(neuron_id))
    
    
    else:
        min_dt = np.floor(min(dt_all_trials))    
        
        if min_dt % 2 == 0:
            min_dt -= 1
        
        bins=np.arange(min_dt, np.ceil(max(dt_all_trials)) + INHIBITORY_SPIKE_RESOLUTION, INHIBITORY_SPIKE_RESOLUTION)    
        
        hist, bins = np.histogram(dt_all_trials, bins=bins)
        hist = hist / float(num_successful_trials)
        
        return hist, bins

def average_inh_input_burst_difference_all_neurons(neurons, num_trials):
    """
    Computes average time difference between inhibitory inputs and dendritic burst
    over trials num_trials and all neurons in neurons
    """
    hist_all_neurons = defaultdict(float)
    
    for n in neurons:
        hist, bins = average_inh_input_burst_difference_one_neuron(n, num_trials)
    
        for i in range(hist.shape[0]):
            if bins[i] in hist_all_neurons:
                hist_all_neurons[bins[i]] += hist[i]
            else:
                hist_all_neurons[bins[i]] = hist[i]
    
    #print hist_all_neurons
    # convert to numpy arrays
    hist = np.array(hist_all_neurons.values())
    bins = np.array(hist_all_neurons.keys())
    
    hist = hist / float(len(neurons)) # compute an average
      
    

    p = bins.argsort() # sort bin array
    
    hist = hist[p]
    bins = bins[p]
    
    bins = np.append(bins, max(bins) + INHIBITORY_SPIKE_RESOLUTION)
    
    #print hist, bins

    plt.figure()    
    
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=INHIBITORY_SPIKE_RESOLUTION)
    plt.title("Histogram for time difference between inhibitory inputs and dendritic burst")
    plt.xlabel("time difference (ms)")
    plt.ylabel("# of inhibitory inputs")
    #for dt in dt

num_trials = 1

#==============================================================================
# # gaba Maturation 300117
# neurons = [71, 186, 187, 84, 44, 175, 219, 238, 224, 70, 288, 117, 99, 276, 23, 24, 165, \
#                                   128, 184, 155, 114, 203, 257, 65, 273, 183, 294, 19, 35, 97, 142, 233, 6, 192, \
#                                   248, 295, 38, 69, 207, 268, 49, 263, 132, 101, 33, 206, 90, 252, 77, 43, 293, 36, \
#                                   5, 180, 282, 65, 34, 267, 208, 66, 146, 179]
#                                   
#==============================================================================

#==============================================================================
# # gabaMaturation 310117
# neurons = [201, 209, 124, 275, 40, 87, 66, 282, 222, 285, 115, 58, 183, 123, 244, 96, 226,\
#            110, 15, 20, 178, 215, 192, 128, 280, 38, 7, 235, 273, 258, 227, 132, 169, 172, \
#                                   243, 100, 188]
#                                   
#==============================================================================

#==============================================================================
# # gabaMaturation 010217
# neurons = [179, 129, 128, 268, 130, 142, 15, 115, 273, 19, 23, 282, 29, 261, 290, \
#            163, 292, 37, 167, 168, 169, 199, 172, 51, 182, 60, 68, 69, 256, 201, 207, 208, 209, \
#            82, 85, 87, 90, 92, 122, 144, 226, 227, 131, 101, 81, 259, 231, 110, 114, 243, 117, 120, \
#            250, 123, 124, 213]
#==============================================================================

num_layers = 11
num_neurons_in_layer = 4
neurons = [n*num_neurons_in_layer + j for n in range(1, num_layers) for j in range(num_neurons_in_layer) ]    

print neurons       
           
average_inh_input_burst_difference_all_neurons(neurons, num_trials)

plt.show()
