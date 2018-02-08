# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 18:20:55 2018

@author: jingroup

Scipt contains useful functions
"""
import numpy as np
import bisect
import os
import reading


def get_input_weight_time_sequence(dirname, end_trial, trialStep, source_neurons, target_neurons):
    """
    Outputs total excitatory conductance input from source_neurons to target_neurons at different
        points during simulation
    """
    num_timepoints = end_trial / trialStep + 1
    
    input_weights = np.zeros((len(target_neurons), num_timepoints), np.float32)

    current_trial = 0
    timepoint = 0

    while current_trial <= end_trial:   
        fileWeights = os.path.join(dirname, "weights_" + str(current_trial) + ".bin")
        (_, _, weights) = reading.read_weights(fileWeights)
        
        for source_id in source_neurons:
            for i, target_id in enumerate(target_neurons):
                input_weights[i][timepoint] += weights[source_id][target_id]
                
        timepoint += 1
        current_trial += trialStep
        
    return input_weights
        
            
            
def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect.bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError

def get_event_times_neurons(filename, neurons):
    """
    Reads event times (spikes in somatic or dendritic compartments )of neurons in trial
    !!! Array neurons has to be sorted!!!
    """
    (_, _, event_times, neuron_fired) = reading.read_time_info(filename)
    
    event_times_neurons = [[] for i in range(len(neurons))]
    
    for neuron_event_times, neuron_fired_id in zip(event_times, neuron_fired):
        try: 
            ind = index(neurons, neuron_fired_id[0])
            event_times_neurons[ind] = neuron_event_times
            
        except ValueError:
            continue
            
    return event_times_neurons
