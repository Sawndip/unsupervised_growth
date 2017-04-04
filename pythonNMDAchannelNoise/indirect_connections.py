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

# default set of training neurons
training_neurons = [0, 1, 2, 3]

# threshold for inhibition that is considered as large enough to drive 
# excitatory immature cells
enough_inhibition = 1.0 

# threshold for inhibition that is considered as too big and would not result
# in target recruitment
too_large_inhibition = 4.0

convergent_inputs_map = {1 : 1, 2 : 2, 3 : 2, 4 : 2}

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

def find_output_connections_to_neurons_not_in_chain(layer_neurons, neurons_in_chain, RA2RA_targets):
    """
    Outputs connections of neurons in the chain layer to the next layer.
    Connections to neurons in the previous layer and to neurons inside the same
    chain layer are ignored
    """
    neurons_in_next_layer = []
    
    for neuron in layer_neurons:
        for target in RA2RA_targets[neuron]:
            if target not in neurons_in_chain and target not in neurons_in_next_layer:
                neurons_in_next_layer.append(target)

    return neurons_in_next_layer

def calculate_effective_inhibition(num_convergent_inputs_to_interneurons, Gie):
    """
    Calculate effective inhibitory input to neurons that receive connections
    from interneurons. Each interneuron receives convergent inputs from source 
    neurons defined in num_convergent_inputs_to_interneurons
    Inhibitory connection strength is Gie
    """
    effective_inhibition = []
    
    for inhibitory_inputs in num_convergent_inputs_to_interneurons:
        if len(inhibitory_inputs) > 0:
            total_inhibitory_input = 0.0            
            for inhibitory_input in inhibitory_inputs:
                total_inhibitory_input += convergent_inputs_map[inhibitory_input] * Gie
                
        else:
            total_inhibitory_input = 0.0
            
        effective_inhibition.append(total_inhibitory_input)
        
    return effective_inhibition

def create_chain(num_layers):
    """
    Creates artificial chain
    """            
    chain = [list(training_neurons)]
    
    RA2I_targets = [] # connections from HVC(RA) to HVC(I) neurons
    I2RA_targets = [] # connections from HVC(I) to HVC(RA) neurons
    
    num_in_layer = len(training_neurons) # number of neurons in the chain layer    

    current_neuron = 0

    #RA2I_targets.append([0])
    
    for i in range(num_in_layer):
        RA2I_targets.append([i])
        I2RA_targets.append([i + num_in_layer])
        #I2RA_targets[-1].append(i + 2*num_in_layer)
    #I2RA_targets[-1].extend([j+2*num_in_layer for j in range(num_in_layer)])
    
    for i in range(num_layers - 1):
        current_neuron += num_in_layer
        chain.append(range(current_neuron, current_neuron + num_in_layer))

        #RA2I_targets.append([i+1])
    
        for j in range(num_in_layer):
            RA2I_targets.append([j + current_neuron])
        
            I2RA_targets.append([j + current_neuron + 1 + 0*num_in_layer])
        #I2RA_targets.append([j+(i+2)*num_in_layer for j in range(num_in_layer)])
        #I2RA_targets[-1].extend([j+(i+3)*num_in_layer for j in range(num_in_layer)])
            
    print chain
    print RA2I_targets
    print I2RA_targets
  
    return (chain, RA2I_targets, I2RA_targets)  
  
def find_chain_formation_mechanism(chain, RA2I_targets, I2RA_targets, starting_layer, Gie):
    """
    Count three types of inhibition:
    1.) Inhibition received from the previous layer
    2.) Inhibition received from the layers prior to previous (layers 2 and 3 before)
    3.) Inhibition received from the neurons in the same group
    """
    # number of several previous layers that may contribute to recruitment mechanism
    num_several_previous_layers = 2
    
    # number of neurons recruited due to inhibitory input from the previous chain layer
    num_recruited_due_to_previous = 0    
    
    # number of neurons recruited due to inhibitory input from chain layers before previous
    num_recruited_due_to_layers_before_previous = 0    
    
    # number of neurons recruited due to mixture of inhibitory input from all previous chain layer
    num_recruited_due_to_mixture_previous = 0        
    
    # number of neurons recruited due to inhibitory input from the same chain layer
    num_recruited_due_to_same = 0  
    
    # number of neurons recruited due to mixture of three mechanisms
    num_recruited_due_to_mixture_of_all_mechanisms = 0
    
    # number of neurons recruited due to unknown mechanism
    num_recruited_due_to_unknown_mechanism = 0
    
    # loop through the chain starting at layer starting_layer
    for i in range(starting_layer-1, len(chain)):
        print "Neurons in the layer: ",chain[i]
        # evaluate different mechanisms
        # first check if inhibitory input from layers before prvious is strong enough
        # to recruit neurons
        inhibition_from_layers_before_previous = [0.0] * len(chain[i])            
            
        for j in range(2, num_several_previous_layers+2):   
            # make sure that we don't go beyond the chain
            if i - j >= 0:
                connections_from_some_previous_layer, num_connections_from_some_previous_layer, num_convergent_inputs_to_interneurons_connected_to_targets = \
                            find_indirect_connections_to_target(chain[i-j], chain[i], RA2I_targets, I2RA_targets)
                
                inhibition_from_some_layer_before_previous = calculate_effective_inhibition(num_convergent_inputs_to_interneurons_connected_to_targets, Gie)
                
                for k in range(len(inhibition_from_some_layer_before_previous)):
                    inhibition_from_layers_before_previous[k] += inhibition_from_some_layer_before_previous[k]
        
       
        print "Inhibition from layers before previous: ",inhibition_from_layers_before_previous
        
    
        # neurons that are likely to be recruited due to other reasons
        neurons_recruited_for_other_reasons = []        
        inhibition_for_neurons_recruited_for_other_reasons = []
        
        for k, inhibition in enumerate(inhibition_from_layers_before_previous):
            if inhibition < enough_inhibition:
                neurons_recruited_for_other_reasons.append(chain[i][k])
                inhibition_for_neurons_recruited_for_other_reasons.append(inhibition)
                    
        num_recruited_due_to_layers_before_previous += len(chain[i]) - len(neurons_recruited_for_other_reasons)
        
        if len(neurons_recruited_for_other_reasons) > 0:
            print "Neurons recruited for other reasons: ",neurons_recruited_for_other_reasons
                    
            # for neurons recruited due to other reasons
            # check if they receive strong input previous layers
            connections_from_previous_layer, num_connections_from_previous_layer, num_convergent_inputs_to_interneurons_connected_to_targets = \
                            find_indirect_connections_to_target(chain[i-1], neurons_recruited_for_other_reasons, RA2I_targets, I2RA_targets)
                
            inhibition_from_previous_layer = calculate_effective_inhibition(num_convergent_inputs_to_interneurons_connected_to_targets, Gie)
        
            print "Inhibition from previous layer: ",inhibition_from_previous_layer
            
            # update list of neurons recruited due to other mechanism            
            neurons_to_check_for_same_layer_connections = []            
            inhibition_for_neurons_to_check = []
            
            for k in range(len(neurons_recruited_for_other_reasons)):
                # if inhibition from several previous is strong enough on its own
                if inhibition_from_previous_layer[k] >= enough_inhibition:
                    num_recruited_due_to_previous += 1
                # if combination of inhibition from previous and inhibition from several previous
                # is enough
                elif inhibition_for_neurons_recruited_for_other_reasons[k] + inhibition_from_previous_layer[k] >= enough_inhibition:
                    num_recruited_due_to_mixture_previous += 1
                # if inhibition is still not strong enough, add neurons to the new list
                else:    
                    #print "Different mechanism for neuron ",neurons_recruited_for_other_reasons[i]
                    neurons_to_check_for_same_layer_connections.append(neurons_recruited_for_other_reasons[k])
                    inhibition_for_neurons_to_check.append(inhibition_for_neurons_recruited_for_other_reasons[k] + inhibition_from_previous_layer[k])
              
            # check if remaining neurons receive strong inhibition from the same layer
            if len(neurons_to_check_for_same_layer_connections) > 0:
                print "Neurons that remain to be checked for inhibitory connections within layer: ",neurons_to_check_for_same_layer_connections
                connections_between_neurons, num_connections_between_neurons,  num_convergent_inputs_to_interneurons_connected_to_targets = \
                        find_indirect_connections_between_neurons(chain[i], RA2I_targets, I2RA_targets)
                        
                all_inhibition_from_same_layer = calculate_effective_inhibition(num_convergent_inputs_to_interneurons_connected_to_targets, Gie)
                
                indices_of_neurons_to_check = [chain[i].index(neuron_to_check) for neuron_to_check in neurons_to_check_for_same_layer_connections]
                
                inhibition_from_same_layer = [all_inhibition_from_same_layer[index] for index in indices_of_neurons_to_check]           
                
                print "Inhibition from same layer: ",inhibition_from_same_layer
            
                for k in range(len(neurons_to_check_for_same_layer_connections)):
                    if inhibition_from_same_layer[k] >= enough_inhibition:
                        num_recruited_due_to_same += 1
                    elif inhibition_from_same_layer[k] +  inhibition_for_neurons_to_check[k] >= enough_inhibition:
                        num_recruited_due_to_mixture_of_all_mechanisms += 1
                    else:
                        num_recruited_due_to_unknown_mechanism  += 1
                        print "Mechanism for recruting neuron {0} is not known!".format(neurons_to_check_for_same_layer_connections[k])

    print "Number of neurons recruited due to inhibition from the previous layer: ",num_recruited_due_to_previous
    print "Number of neurons recruited due to inhibition from layers before previous: ",num_recruited_due_to_layers_before_previous
    print "Number of neurons recruited due to mixture of inhibition from all previous layers: ", num_recruited_due_to_mixture_previous
    print "Number of neurons recruited due to inhibition from same layer: ", num_recruited_due_to_same
    print "Number of neurons recruited due to mixture of all mechanisms: ", num_recruited_due_to_mixture_of_all_mechanisms
    print "Number of neurons recruited due to unknown mechanism: ", num_recruited_due_to_unknown_mechanism
    
    return (num_recruited_due_to_previous, num_recruited_due_to_layers_before_previous, num_recruited_due_to_mixture_previous, \
        num_recruited_due_to_same, num_recruited_due_to_mixture_of_all_mechanisms, num_recruited_due_to_unknown_mechanism)
    
def find_chain_layers(RA2RA_targets):
    """
    Find which neurons are in each chain layer
    """
    chain = [list(training_neurons)] # chain structure
    
    # start from training neurons    
    source_layer = list(training_neurons)
    neurons_in_chain = set(list(training_neurons))
    
    neurons_in_next_layer = find_output_connections_to_neurons_not_in_chain(source_layer, neurons_in_chain, RA2RA_targets)

    print neurons_in_next_layer    
    
    while len(neurons_in_next_layer) > 0:
        neurons_in_chain.update(neurons_in_next_layer)
        chain.append(neurons_in_next_layer)
        source_layer = neurons_in_next_layer
        neurons_in_next_layer = find_output_connections_to_neurons_not_in_chain(source_layer, neurons_in_chain, RA2RA_targets)
        print neurons_in_next_layer   
        
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

def find_indirect_connections_between_neurons(neurons, RA2I_targets, I2RA_targets):
    """
    Functions finds all indirect connections between neurons in neurons.
    Self connections are ignored
    """
    connections_between_neurons = []
    num_connections_between_neurons = []

    num_convergent_inputs_to_interneurons_connected_to_targets = []
    
    
    for target in neurons:
        connections_to_target = {}
        connections_to_target[target] = []
        num_total_connections = 0
        interneurons_connected_to_target = set()    
        
        for source in neurons:
            if source != target:
                num_connections_from_source = 0
            
                for interneuron in RA2I_targets[source]:
                    if target in I2RA_targets[interneuron]:
                        interneurons_connected_to_target.add(interneuron)
                        num_connections_from_source += 1
                
                if num_connections_from_source > 0:
                    connections_to_target[target].append((source, num_connections_from_source))
                
                num_total_connections += num_connections_from_source
                 
        source_except_target = [n for n in neurons if n != target]                
        
        #print source_except_target        
        
        num_convergent_inputs_to_interneurons = find_num_connections_to_interneurons(source_except_target, interneurons_connected_to_target, RA2I_targets)     
        #print num_connections_to_interneurons
        num_convergent_inputs_to_interneurons_connected_to_targets.append(num_convergent_inputs_to_interneurons)
        
        connections_between_neurons.append(connections_to_target)
        num_connections_between_neurons.append(num_total_connections)
        
    return connections_between_neurons, num_connections_between_neurons,  num_convergent_inputs_to_interneurons_connected_to_targets
                
        

def find_num_connections_to_interneurons(source, interneurons, RA2I_targets):
    """
    Function finds number of connections that source neurons make on interneurons
    """    
    num_connections_to_interneurons = []    
    
    for interneuron in interneurons:
        num_connections = 0
        for s in source:
            if interneuron in RA2I_targets[s]:
                num_connections += 1
            
        num_connections_to_interneurons.append(num_connections)
        
    return num_connections_to_interneurons

def find_indirect_connections_to_target(source, target, RA2I_targets, I2RA_targets):
    """
    Function finds all indirect connections from source group to target group
    """
    connections_from_previous_layer = []    
    num_connections_from_previous_layer = []
    
    num_convergent_inputs_to_interneurons_connected_to_targets = []   
    
    for t in target:
        connections_to_target_from_previous_layer = {}
        connections_to_target_from_previous_layer[t] = []
        num_total_connections = 0
        interneurons_connected_to_target = set()
        
        for s in source:
            num_connections_from_source = 0
            
            for interneuron in RA2I_targets[s]:
                if t in I2RA_targets[interneuron]:
                    interneurons_connected_to_target.add(interneuron)
                    num_connections_from_source += 1
                
            if num_connections_from_source > 0:
                
                connections_to_target_from_previous_layer[t].append((s, num_connections_from_source))
                
                
            num_total_connections += num_connections_from_source
            #print "num_connections_from_source = ",num_connections_from_source
            #print "num_total_connections = ",num_total_connections

        
        num_convergent_inputs_to_interneurons = find_num_connections_to_interneurons(source, interneurons_connected_to_target, RA2I_targets)     
        #print num_connections_to_interneurons
        num_convergent_inputs_to_interneurons_connected_to_targets.append(num_convergent_inputs_to_interneurons)
        
        connections_from_previous_layer.append(connections_to_target_from_previous_layer)
        num_connections_from_previous_layer.append(num_total_connections)
                
    return connections_from_previous_layer, num_connections_from_previous_layer, num_convergent_inputs_to_interneurons_connected_to_targets

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

def autolabel(rects, name):
    # attach some text labels
    for ii,rect in enumerate(rects):
        height = rect.get_height()
        plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, '%s'% (name[ii]),
                ha='center', va='bottom')
    
           
if __name__ == "__main__":
    
    RA2I = "/home/eugene/Output/networks/gabaMaturation280317/RA_I_connections.bin"
    I2RA = "/home/eugene/Output/networks/gabaMaturation280317/I_RA_connections.bin"
    RA2RA = "/home/eugene/Output/networks/gabaMaturation280317/RA_RA_super_connections.bin"
    fileMature = "/home/eugene/Output/networks/gabaMaturation280317/mature.bin"
    fileBursts = "/home/eugene/Output/networks/gabaMaturation280317/spike_times_dend.bin"
    
    (N_RA, RA2I_targets, RA2I_targets_G) = reading.read_connections(RA2I)
    (N_RA, RA2RA_targets, RA2RA_targets_G) = reading.read_connections(RA2RA)
    (N_I, I2RA_targets, I2RA_targets_G) = reading.read_connections(I2RA)
      
    
    for targets_G in I2RA_targets_G:
        if len(targets_G) > 0:
            # inhibitory input strength
            Gie = targets_G[0]
            break
    
    print Gie
    
    
    
    #print I2RA_targets
    #print I2RA_targets
    #print RA2RA_targets
            
    training = [0, 1, 2, 3]
    num_layers = 5
    
    #indirect_connections_in_groups(training, RA2I_targets, RA2RA_targets, I2RA_targets, num_layers)
    
    #print get_indirect_inhibition([0, 1, 2, 3], RA2I_targets, I2RA_targets, G)
    
    
    
    #(chain, RA2I_targets, I2RA_targets) = create_chain(num_layers)   
    
    
    chain = find_chain_layers(RA2RA_targets)    
    print chain
    
    starting_layer = 2 
    num_mechanisms = 6
    
    
    freq = find_chain_formation_mechanism(chain, RA2I_targets, I2RA_targets, starting_layer, Gie)
    
    name = ["previous layer", "before previous", "mixed previous", "same layer", "mixed three", "unknown"]    
      
    
    index = np.arange(1, num_mechanisms+1)
    bar_width = 0.5 
    rects = plt.bar(index, freq, bar_width)
    plt.xlim([0, num_mechanisms + 1])
    plt.ylim([0, max(freq)+1])
    plt.ylabel("number of recruited neurons")
   
    plt.title("Number of neurons recruited by different mechanisms")
    plt.xticks([])
    autolabel(rects, name)
    
    plt.show()       
    
    #num_indirect_connections_from_previous_layer = indirect_connections_in_groups(training, RA2I_targets, RA2RA_targets, I2RA_targets, num_layers)
    
    #print num_indirect_connections_from_previous_layer
    
    source = [0, 1, 2, 3]
    targets = [148, 173, 129, 287]
    