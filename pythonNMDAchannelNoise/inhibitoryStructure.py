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


dataDir = "/home/eugene/Output/matureTest/gabaMaturation130417/RA/"

DENDRITIC_THRESHOLD = 0.0 # threshold for burst spike in mV
INHIBITORY_KICK_THRESHOLD = 0.1 # threshold for inhibitory inputs
INHIBITORY_SPIKE_RESOLUTION = 1.0 # resolution in ms wihin which we treat inhibitory spike as the same
TIMESTEP = 0.02 # timestep of neuron dynamics in ms
WINDOW_SIZE = 50.0 # window in which to calculate input conductances

def find_silent_neurons(neurons, num_trials, dataDir):
    """
    Finds neurons that did not fire in any of the training trials
    """
    num_successful_trials = [] # number of trials in which dendritic bursts were produced    
    
    for neuron_id in neurons:
        counter = 0 # counter to count successful trials
        for n in range(num_trials):
            filename = dataDir + "RA" + str(neuron_id) + "_trial" + str(n+1) + ".bin"
            (t, Vs, _, _, _, Vd, _, _, _, _, _, _, _, _, _, _, _, _) = reading.read_hh2(filename)
            
            burst_times, burst_indices = get_burst_times(t, Vd) # obtain dendritic burst times
            
            if len(burst_times) > 0:
                counter += 1        
        num_successful_trials.append(counter)
    
    not_fired = [n for i, n in enumerate(neurons) if num_successful_trials[i] == 0]    
    
    return not_fired, num_successful_trials

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
    #print half_window
    #print num_points_in_window
    
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
        
    #print exc_input.shape
    
    t_input = np.linspace(-half_window * TIMESTEP, half_window * TIMESTEP, num_points_in_window)

    # if no dendritic bursts was produced            
    if np.sum(successful_trials) == 0:        
        average_exc_input = np.zeros(num_points_in_window, np.float32)
        average_inh_input = np.zeros(num_points_in_window, np.float32)
        return t_input, average_exc_input, average_inh_input
    
    else:    
        average_exc_input = np.sum(exc_input[successful_trials], axis=0) / np.sum(successful_trials)
        average_inh_input = np.sum(inh_input[successful_trials], axis=0) / np.sum(successful_trials)
        
        return t_input, average_exc_input, average_inh_input

def get_burst_times(t, Vd):
    """
    Extracts dendritic burst times from voltage data
    """
    burst_times = []
    burst_indices = []
    flag = False # indicator of crossing potential threshold
        
    for i, v in enumerate(Vd):
        if flag == False and v >= DENDRITIC_THRESHOLD:
            flag = True
            burst_times.append(t[i])
            burst_indices.append(i)
            
        elif flag == True and v < DENDRITIC_THRESHOLD:
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
    
           
#average_inh_input_burst_difference_all_neurons(neurons, num_trials)


num_trials = 20


#==============================================================================
# # gabaMaturation 300117
# neurons_all = [155, 128, 165, 114, 184, 273, 203, 257, 183, 294, 19, 35, 97, 142, 6, 233, \
#                192, 295, 38, 248, 69, 207, 49, 263]
#==============================================================================

#==============================================================================
# # gabaMaturation 310117
# neurons_all = [96, 226 , 110, 15, 178, 20, 215, 192, 128, 7, 38, 280]
# 
#==============================================================================

#==============================================================================
# # gabaMaturation 010217
# neurons_all = [131, 69, 90, 15, 169, 179, 226, 123, 207, 29, 120, 231, 243, 117, 172, 163, \
#                60, 167, 142]
#==============================================================================

#==============================================================================
# #RandomChainConnections240217
# neurons_all = [225, 198, 204, 114, 193, 78, 47, 177, 206, 22, 110, 264, 42, 169, 155, 299, 235, \
#            118, 89, 255, 90, 178, 251, 253, 98, 167, 51, 145, 201, 75, 190, 148, 72, 268, 194, 105, 49, 86, 246, \
#            184, 27, 109, 216, 173]
#==============================================================================

#==============================================================================
# #RandomChainConnections010317
# neurons_all = [288, 12, 197, 172, 264, 146, 136, 9, 187, 259, 101, 243, 161, 247, 44, 67, 143, 200, 257, 17, \
#                212, 34, 46, 214, 226, 240, 294, 160, 127, 236, 15, 196, 195, 88, 41, 230, 69, 42, 292, 262, \
#                98, 297, 144, 266, 108, 179, 272, 50, 233, 102, 271, 206, 103, 26, 232, 121, 248, 120, 215, 205, \
#                36, 254, 92, 148]
#==============================================================================


#neurons_all = [174, 216, 27, 185, 152, 97, 101, 181, 165, 14, 234, 60, 213, 94, 287, 155, 109, 170, \
#                42, 265, 110, 22, 194, 207, 115, 177]

#==============================================================================
# not_fired, num_successful_trials = find_silent_neurons(neurons_all, num_trials, dataDir)
#         
# print("neurons that did not fire: ",not_fired)
# print("number of successful trials for each neuron:",num_successful_trials)
# 
#==============================================================================

#==============================================================================
# # gabaMaturation 270317 huxley
# #neurons_all = [179, 66 , 11, 123, 173, 129, 148, 287, 199, 174, 285, 298, 144, 20, 161, 165, 205, 89, 17]
# #neurons_all = [199, 174, 285, 298, 144, 20, 161, 165, 205, 89, 17]
# 
#==============================================================================

# gabaMaturation 010417 hodgkin (skip three layers)
#==============================================================================
# neurons_all = [115, 194, 225, 78, 268, 221, 289, 104, 185, 285, 287, 21, 58, 55, 229, \
#                222, 145, 239, 123, 173, 295, 179, 240, 134, 280, 42, 228, 178, 208, 244,\
#                294, 130, 45, 4, 217, 143, 87, 226, 148, 233, 190, 223, 255, 138, 29, 192,\
#                290, 12, 142, 129, 150, 48, 69, 271, 174, 17, 167, 168, 273, 68, 35, 95,\
#                163, 207, 128, 172, 231, 258, 99, 30, 100]
# 
#==============================================================================


#==============================================================================
# #gabaMaturation 280317 huxley (skip three layers)
# #neurons_all = [134, 84, 38, 72, 267,\
# #               34, 137, 77, 20, 188, 200, 136, 173, 13, 206, 5, 118]
# 
# 
# neurons_all = [20, 188, 200, 136, 173, 13, 206, 5, 118]
#==============================================================================

#==============================================================================
# # gabaMaturation040417 huxley (skip three layers)
# neurons_all = [64, 140, 174, 210, 236, \
#                 77, 129, 15, 39, 298, 168, 216, 142, 295, 204, 13, 23, 34, 280, 186, 299, 121, 54, \
#                 269, 292, 105, 9, 35, 57, 251, 100, 69, 260, 182, 136, 237, 134, 26, 66, 157, 286, \
#                 135, 193, 45, 219, 80, 20, 126, 196, 211, 6, 190, 257, 81, 104, 36, 253, 25, 90, 115,\
#                 30, 183, 63, 109, 266, 202, 94, 113, 222, 187, 246, 86, 206, 232, 160, 125, 240, 117, 282,\
#                 152, 19, 259, 198, 128]
#==============================================================================

#gabaMaturation130417 huxley 
neurons_all = [51, 48, 146, 172, 132, 277, 203, 175, 275, 28, 31, 37, 140, 235, 67, 245, 21, 50, 138, 93,\
                 76, 228, 46, 225, 187, 231, 156, 210, 246, 148, 7, 49, 195, 74, 124, 255, 169, 152, 269, 206, \
                 260, 94, 83, 259, 57, 171, 114, 23, 222, 248, 113, 165, 20, 104, 116, 59, 257, 25, 26, 89, 252,\
                 151, 229, 253, 106, 176, 115, 183, 283, 30, 112, 226, 267, 139, 238, 158, 167, 95, 84, 268, 162,\
                 111, 164, 163]

#neurons_all = range(4, 4+4*20)


# find chain layers


# calculate inhibitory conductance of chain neurons

neuron_1 = 51
neuron_2 = 163
#t_input, average_exc_input, average_inh_input = calculate_average_input_conductance_for_neuron(neuron_id, num_trials, dataDir)
t_input, average_exc_input, average_inh_input = calculate_average_input_conductance_for_neurons(neurons_all, num_trials, dataDir)
_, average_exc_input_1, average_inh_input_1 = calculate_average_input_conductance_for_neuron(neuron_1, num_trials, dataDir)
_, average_exc_input_2, average_inh_input_2 = calculate_average_input_conductance_for_neuron(neuron_2, num_trials, dataDir)


f = plt.figure()

ax1 = f.add_subplot(211)
ax1.set_title("Average conductances of {0} neurons from gabaMaturation130417 huxley".format(len(neurons_all)))

#ax1.set_title("Conductances for neuron {0}".format(neuron_id))
ax1.plot(t_input, average_exc_input)
ax1.set_ylabel("$G_{exc, d}$")
_, ymax = ax1.get_ylim()
ax1.set_ylim([0, ymax])

ax2 = f.add_subplot(212)

ax2.plot(t_input, average_inh_input)
ax2.set_xlabel("time (ms)")
ax2.set_ylabel("$G_{inh, d}$")
_, ymax = ax2.get_ylim()
ax2.set_ylim([0, ymax])

f = plt.figure()
ax1 = f.add_subplot(221)

ax1.set_title("Conductances for neuron {0}".format(neuron_1))
ax1.plot(t_input, average_exc_input_1)
ax1.set_ylabel("$G_{exc, d}$")
_, ymax = ax1.get_ylim()
ax1.set_ylim([0, ymax])

ax2 = f.add_subplot(223)

ax2.plot(t_input, average_inh_input_1)
ax2.set_xlabel("time (ms)")
ax2.set_ylabel("$G_{inh, d}$")
_, ymax = ax2.get_ylim()
ax2.set_ylim([0, ymax])

ax1 = f.add_subplot(222)

ax1.set_title("Conductances for neuron {0}".format(neuron_2))
ax1.plot(t_input, average_exc_input_2)
ax1.set_ylabel("$G_{exc, d}$")
_, ymax = ax1.get_ylim()
ax1.set_ylim([0, ymax])

ax2 = f.add_subplot(224)

ax2.plot(t_input, average_inh_input_2)
ax2.set_xlabel("time (ms)")
ax2.set_ylabel("$G_{inh, d}$")
_, ymax = ax2.get_ylim()
ax2.set_ylim([0, ymax])


plt.show()



