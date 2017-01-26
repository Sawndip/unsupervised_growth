# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 18:09:22 2016

@author: jingroup

Script analyzes how many indirect connections 
neurons in one group make onto neurons in the next group
"""

RA2I = "/home/eugene/Output/RA_I_connections.bin"
I2RA = "/home/eugene/Output/I_RA_connections.bin"
RA2RA = "/home/eugene/Output/RA_RA_super_connections.bin"
fileMature = "/home/eugene/Output/mature.bin"

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
    
    for t in target:
        for s in source:
            num_connections = 0
            for interneuron in RA2I_targets[s]:
                if t in I2RA_targets[interneuron]:
                    num_connections += 1
                
            if num_connections > 0:
                d[t].append((s, num_connections))
                
    return d

def indirect_connections_in_groups(training, RA2I_targets, RA2RA_targets, I2RA_targets, num_layers):
    """
    Function finds indirect connections in num_layers groups of neurons. 
    It starts from training neurons and iterates along the feedforward chain.
    """
    
    source = set(list(training))    
    target = set()

    for i in xrange(num_layers):
        # find target neurons
        for s in source:
            target.update(RA2RA_targets[s])
            
        indConnections = get_indirect_connections(source, target, RA2I_targets, I2RA_targets)
        print "layer ", i
        print indConnections
        
        source = target
        target = set()
        
        print source
        print target
            

(N_RA, RA2I_targets, RA2I_targets_G) = reading.read_connections(RA2I)
(N_RA, RA2RA_targets, RA2RA_targets_G) = reading.read_connections(RA2RA)
(N_I, I2RA_targets, I2RA_targets_G) = reading.read_connections(I2RA)

#print RA2I_targets
#print I2RA_targets
#print RA2RA_targets


        
training = [0, 1, 2, 3]
targets = [13, 16, 18, 25, 28, 46, 47, 52, 56, 62, 67, 71, 75, 84, 102, 103, 105, 113, 119, 121, 141, 147, 151, 156, 157, 171, 175, 176, 186, 187, 189, 218, 219, 220, 238, 242, 245, 246, 247, 262, 264, 265, 269, 278, 281, 289, 291, 298]
num_layers = 2

#indirect_connections_in_groups(training, RA2I_targets, RA2RA_targets, I2RA_targets, num_layers)
G = 0.2

print get_indirect_inhibition([0, 1, 2, 3], RA2I_targets, I2RA_targets, G)

#mature = reading.read_mature(fileMature)

#print "mature = ", mature

#mature_id = [i for i,e in enumerate(mature) if e == 1]

#print "Mature neurons: ", mature_id

#print get_indirect_connections(training, mature_id, RA2I_targets, I2RA_targets)