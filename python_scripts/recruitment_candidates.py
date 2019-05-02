# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 12:42:25 2017

@author: jingroup

Scripts shows canditate neurons for recruitment
"""

import reading
import indirect_connections as connections
import os


if __name__ == "__main__":
    
    networkDir = "/home/eugene/Output/networks/sphere"
    
    RA2I = os.path.join(networkDir, "RA_I_connections_initial.bin")
    I2RA = os.path.join(networkDir, "I_RA_connections_initial.bin")
    
    training_neurons = [131, 95, 277, 179]
    map_inputs_to_spikes = {0 : 0, 1 : 1, 2 : 2, 3 : 2, 4 : 3} # map of number of convergent inputs from HVC(RA) to interneuron spikes
    
    inh_enough = 1.0 # strength of inhibition that is strong enough for recruitment
    Gie = 0.10 # inhibition strength
    
    (N_RA, RA2I_targets, RA2I_targets_G) = reading.read_connections(RA2I)
    (N_I, I2RA_targets, I2RA_targets_G) = reading.read_connections(I2RA)
    
    for i in range(N_RA):
        if i not in training_neurons:
            con, num_connections, num_convergent_inputs_to_interneurons = \
                connections.find_indirect_connections_to_target(training_neurons, [i], RA2I_targets, I2RA_targets)
            
            if num_connections[0] > 0:            
                effective_inh = 0
                
                for num_conv_inputs in num_convergent_inputs_to_interneurons[0]:
                    effective_inh += map_inputs_to_spikes[num_conv_inputs]
                
                if effective_inh * Gie >= inh_enough:
                    print "Target neuron %s" % i
                    print "Connections from training neurons: ", con
                    print "Number of convergent inputs to interneurons connected to target: ",num_convergent_inputs_to_interneurons
               
                    print "Effective inhibitory input: ",effective_inh * Gie