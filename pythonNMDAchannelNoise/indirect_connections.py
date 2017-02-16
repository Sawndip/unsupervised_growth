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
    
    RA2I = "/home/eugene/Output/RA_I_connections.bin"
    I2RA = "/home/eugene/Output/I_RA_connections.bin"
    RA2RA = "/home/eugene/Output/RA_RA_super_connections.bin"
    fileMature = "/home/eugene/Output/mature.bin"

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
    
    plt.show()
    #mature = reading.read_mature(fileMature)
    
    #print "mature = ", mature
    
    #mature_id = [i for i,e in enumerate(mature) if e == 1]
    
    #print "Mature neurons: ", mature_id
    
    #print get_indirect_connections(training, mature_id, RA2I_targets, I2RA_targets)