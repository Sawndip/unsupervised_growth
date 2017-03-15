# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 18:09:22 2016

@author: jingroup

Script analyzes how many indirect connections 
neurons in one group make onto neurons in the next group
"""
import reading
import math
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

training_neurons = [0, 1, 2, 3]

def find_all_indirect_connections_from_previously_fired(neurons, burst_times, neurons_fired, RA2I_targets, I2RA_targets):
    """
    Find how many indirect connections each neuron in neurons receives from all previously 
    fired neurons
    """
    num_indirect_from_previous = []
    for n in neurons:
        previously_fired, num_connections_from_previously_fired = find_indirect_connections_from_previously_fired(n, burst_times, neurons_fired, RA2I_targets, I2RA_targets)
        num_indirect_from_previous.append(num_connections_from_previously_fired)
        #print previously_fired
        #print num_connections_from_previously_fired

    # find minimum length of previously fired neurons
    min_length = 100000
    
    for n in num_indirect_from_previous:
        if len(n) < min_length:
            min_length = len(n)

    # sum and average all numbers of indirect connections
    average_num_from_previous = np.zeros(min_length, np.float32) 

    for num_indirect in num_indirect_from_previous:
        for i in range(min_length):
            average_num_from_previous[i] += num_indirect[i]
            
    # normalize
    for i in range(min_length):
        average_num_from_previous[i] /= float(len(num_indirect_from_previous))
        
    print average_num_from_previous
        

def find_indirect_connections_from_previously_fired(neuron_id, burst_times, neurons_fired, RA2I_targets, I2RA_targets):
    """
    Find how many indirect connections a neuron receives from all previously 
    fired neurons
    """
    print neuron_id
    index_in_fired_array = neurons_fired.index(neuron_id)
    print index_in_fired_array
    
    num_connections_from_previously_fired = []
    previously_fired = neurons_fired[index_in_fired_array-1::-1]
    
    print previously_fired
    
    # loop through all neurons that fired before
    for source in previously_fired:       
        num_indirect_connections = 0 # number of indirect connections from source neuron
        
        #print source
        #print RA2I_targets[source]
        # loop through all inhibitory targets of the neuron that fired before
        for inh_target in RA2I_targets[source]:
            # check if neuron neuron_id is among targets of interneuron
            if neuron_id in I2RA_targets[inh_target]:
                num_indirect_connections += 1
    
        num_connections_from_previously_fired.append(num_indirect_connections)
        
    return previously_fired, num_connections_from_previously_fired

def find_chain_layers(RA2RA_targets):
    """
    Find which neurons are in each chain layer
    """
    chain = [] # chain structure

    # start from training neurons    
    source_layer = set(list(training_neurons))
    
    while len(source_layer) > 0:
        chain.append(list(source_layer))
        new_source_layer = set()
        for i in source_layer:
            new_source_layer.update(RA2RA_targets[i])
        source_layer = new_source_layer
    
    return chain

def get_indirect_inhibition(source, RA2I_targets, I2RA_targets, G):
    """
    Function finds all indirect connections from source neurons
    """
    d = defaultdict(list)    
    
    for s in source:
        for interneuron in RA2I_targets[s]:
            for target in I2RA_targets[interneuron]:
                if target in d:
                    d[target] += G
                else:
                    d[target] = G
                
    return d


def get_indirect_connections_to_target(source, target, RA2I_targets, I2RA_targets):
    """
    Function finds all indirect connections from source group to target group
    """
    d = defaultdict(list)    
    num_connections_from_previous_layer = {}
    
    for t in target:
        num_total_connections = 0
        for s in source:
            num_connections = 0
            for interneuron in RA2I_targets[s]:
                if t in I2RA_targets[interneuron]:
                    num_connections += 1
                
            if num_connections > 0:
                d[t].append((s, num_connections))
                
            num_total_connections += num_connections
            
        num_connections_from_previous_layer[t] = num_total_connections
                
    return d, num_connections_from_previous_layer

def indirect_connections_in_groups(training, RA2I_targets, RA2RA_targets, I2RA_targets, num_layers):
    """
    Function finds indirect connections in num_layers groups of neurons. 
    It starts from training neurons and iterates along the feedforward chain.
    """
    
    source = set(list(training))    
    target = set()
    indirect_connections_from_previous_layer = []
    
    for i in xrange(num_layers):
        # find target neurons
        for s in source:
            target.update(RA2RA_targets[s])
            
        indConnections, num_connections_from_previous_layer = get_indirect_connections_to_target(source, target, RA2I_targets, I2RA_targets)
        print "layer", i, " to layer", i+1
        print "neurons in layer",i,":",source
        print "neurons in layer",i+1,":",target
        print "all connections from previous layer:",indConnections
        print "total number of connections from previous layer",num_connections_from_previous_layer
        
        indirect_connections_from_previous_layer.append(list(num_connections_from_previous_layer.values()))
        source = target
        target = set()
        
        
    return indirect_connections_from_previous_layer
        #print source
        #print target
            
if __name__ == "__main__":
    
    RA2I = "/home/eugene/Output/networks/gabaMaturation130317/RA_I_connections.bin"
    I2RA = "/home/eugene/Output/networks/gabaMaturation130317/I_RA_connections.bin"
    RA2RA = "/home/eugene/Output/networks/gabaMaturation130317/RA_RA_super_connections.bin"
    fileMature = "/home/eugene/Output/networks/gabaMaturation130317/mature.bin"
    fileBursts = "/home/eugene/Output/networks/gabaMaturation130317/spike_times_dend.bin"
    
    (N_RA, RA2I_targets, RA2I_targets_G) = reading.read_connections(RA2I)
    (N_RA, RA2RA_targets, RA2RA_targets_G) = reading.read_connections(RA2RA)
    (N_I, I2RA_targets, I2RA_targets_G) = reading.read_connections(I2RA)
    
    #print RA2I_targets
    #print I2RA_targets
    #print RA2RA_targets
    
    
            
    training = [0, 1, 2, 3]
    num_layers = 3
    
    #indirect_connections_in_groups(training, RA2I_targets, RA2RA_targets, I2RA_targets, num_layers)
    G = 0.2
    
    #print get_indirect_inhibition([0, 1, 2, 3], RA2I_targets, I2RA_targets, G)
    
    num_indirect_connections_from_previous_layer = indirect_connections_in_groups(training, RA2I_targets, RA2RA_targets, I2RA_targets, num_layers)
    
    print num_indirect_connections_from_previous_layer
    
    source = [0, 1, 2, 3]
    targets = [219, 44, 238, 175]
    
    
    print get_indirect_connections_to_target(source, targets, RA2I_targets, I2RA_targets)
    
    f = plt.figure()
    ax = f.add_subplot(111)
    
    for i in range(len(num_indirect_connections_from_previous_layer)):
        x = [i+1]*len(num_indirect_connections_from_previous_layer[i])    
        ax.scatter(x, num_indirect_connections_from_previous_layer[i])
    
    ax.set_title("Number of indirect connections from previous layer")
    ax.set_xlabel("layer id")
    ax.set_ylabel("# indirect connections")
    
    
    
    print find_chain_layers(RA2RA_targets)  
    
    

    (trial_number, simulation_time, burst_times, neurons_fired) = reading.read_time_info(fileBursts)
    
    burst_times = [t for sublist in list(burst_times) for t in sublist]
    neurons_fired = [ind for sublist in list(neurons_fired) for ind in sublist]
    
    
    
    burst_times, neurons_fired = zip(*sorted(zip(burst_times, neurons_fired)))
    
    print burst_times
    print neurons_fired
    
    neurons = [285, 213, 235, 20, 266, 134, 190, 17, 277, 153, 130, 207, 48, 129, 11, 89, 257, 258, 87, 9] 
    
    find_all_indirect_connections_from_previously_fired(neurons, burst_times, neurons_fired, RA2I_targets, I2RA_targets)
    
    
    
    
    plt.show()
    #mature = reading.read_mature(fileMature)
    
    #print "mature = ", mature
    
    #mature_id = [i for i,e in enumerate(mature) if e == 1]
    
    #print "Mature neurons: ", mature_id
    
    #print get_indirect_connections(training, mature_id, RA2I_targets, I2RA_targets)