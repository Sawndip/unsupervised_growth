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

dataDir = "/home/eugene/Output/networks/RandomChainTest240217/RA/"

DENDRITIC_THRESHOLD = 0.0 # threshold for burst spike in mV
INHIBITORY_KICK_THRESHOLD = 0.1 # threshold for inhibitory inputs
INHIBITORY_SPIKE_RESOLUTION = 1.0 # resolution in ms wihin which we treat inhibitory spike as the same
TIMESTEP = 0.02 # timestep of neuron dynamics in ms
WINDOW_SIZE = 50.0 # window in which to calculate input conductances

def calculate_average_input_conductance_for_neurons(neurons, num_trials, dataDir):
    """
    Calculates average input condictances. Average performed for all neurons in neurons
    and for number of trials num_trials. Files with neuron trials are located in
    dataDir
    """
    num_points_in_window = int(WINDOW_SIZE / TIMESTEP) # number of datapoints inside window
    
    # if even, add one more datapoint    
    if num_points_in_window % 2 == 0:
        num_points_in_window += 1
    
    average_exc_input = np.zeros(num_points_in_window, np.float32)
    average_inh_input = np.zeros(num_points_in_window, np.float32)

    # loop through all neurons    
    for neuron_id in neurons:    
        t_input, exc_input, inh_input = calculate_average_input_conductance_for_neuron(neuron_id, num_trials, dataDir)
        
        average_exc_input += exc_input
        average_inh_input += inh_input
        
    average_exc_input /= float(len(neurons))
    average_inh_input /= float(len(neurons))
    
    return t_input, average_exc_input, average_inh_input
    
        

def calculate_average_input_conductance_for_neuron(neuron_id, num_trials, dataDir):
    """
    Calculates average input conductance for neuron neuron_id. Average is taken
    for num_trials and is aligned to neuron's dendritic burst
    """
    num_points_in_window = int(WINDOW_SIZE / TIMESTEP) # number of datapoints inside window
    
    # if even, add one more datapoint    
    if num_points_in_window % 2 == 0:
        num_points_in_window += 1
        
    half_window = num_points_in_window / 2
    print half_window
    print num_points_in_window
    
    exc_input = np.zeros((num_trials,num_points_in_window), np.float32)
    inh_input = np.zeros((num_trials,num_points_in_window), np.float32)
    
    successful_trials = np.zeros(num_trials, np.bool) # array with bools indicating trials with dendritic bursts    
    
    for n in range(num_trials):
        filename = dataDir + "RA" + str(neuron_id) + "_trial" + str(n+1) + ".bin"
        (t, Vs, _, _, _, Vd, _, _, _, _, Gexc_d, Ginh_d, _, _, _, _, _, _) = reading.read_hh2(filename)
        
        burst_times, burst_indices = get_burst_times(t, Vd) # obtain dendritic burst times        
        
        if len(burst_times) > 0:        
            successful_trials[n] = True
            # loop through all dendritic burst times
            for burst_index in burst_indices:
                for i in range(num_points_in_window):
                    exc_input[n][i] += Gexc_d[burst_index - half_window + i]
                    inh_input[n][i] += Ginh_d[burst_index - half_window + i]
                    
            # average over number of dendritic bursts:
            exc_input[n] = exc_input[n] / float(len(burst_indices))
            inh_input[n] = inh_input[n] / float(len(burst_indices))
        # if no dendritic bursts was produced
        else:
            average_exc_input = np.zeros(num_points_in_window, np.float32)
            average_inh_input = np.zeros(num_points_in_window, np.float32)
            t_input = np.linspace(-half_window * TIMESTEP, half_window * TIMESTEP, num_points_in_window)
            
            return t_input, average_exc_input, average_inh_input
            
    print exc_input.shape
    average_exc_input = np.sum(exc_input[successful_trials], axis=0) / np.sum(successful_trials)
    average_inh_input = np.sum(inh_input[successful_trials], axis=0) / np.sum(successful_trials)
    t_input = np.linspace(-half_window * TIMESTEP, half_window * TIMESTEP, num_points_in_window)
        
    return t_input, average_exc_input, average_inh_input

def get_burst_times(t, Vd):
    """
    Extracts dendritic burst times from voltage data
    """
    burst_times = []
    burst_indices = []
    flag = False
        
    for i, v in enumerate(Vd):
        if flag == False and v >= DENDRITIC_THRESHOLD:
            flag = True
        elif flag == True and v < DENDRITIC_THRESHOLD:
            burst_times.append(t[i])
            burst_indices.append(i)
            flag = False
        
    return burst_times, burst_indices

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
    burst_times, _ = get_burst_times(t, Vd)
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
           
#average_inh_input_burst_difference_all_neurons(neurons, num_trials)

neuron_id = 106
num_trials = 20

#neurons = [276, 23, 99, 24, 155, 128, 165, 114, 184, 203, 257, 183, 294, 273]
neurons = [102, 236, 118, 90, 256, 38, 179, 252, 254, 168, 51, 145, 202, 191, 26, 149, \
            72, 269, 106, 49, 86, 62, 247, 91]

#t_input, average_exc_input, average_inh_input = calculate_average_input_conductance_for_neuron(neuron_id, num_trials, dataDir)
t_input, average_exc_input, average_inh_input = calculate_average_input_conductance_for_neurons(neurons, num_trials, dataDir)

print t_input
print t_input.shape
print average_exc_input.shape

f = plt.figure()

ax1 = f.add_subplot(211)
ax1.set_title("Average conductances of 24 neurons from RandomChainTest240217".format(neuron_id))

#ax1.set_title("Conductances for neuron {0}".format(neuron_id))
ax1.plot(t_input, average_exc_input)
ax1.set_ylabel("$G_{exc, d}$")

ax2 = f.add_subplot(212)

ax2.plot(t_input, average_inh_input)
ax2.set_xlabel("time (ms)")
ax2.set_ylabel("$G_{inh, d}$")


plt.show()
