# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 19:29:31 2017

@author: jingroup

Script read and writes training neurons 
"""

import reading
import struct
import numpy as np


def write_training_neurons(training_neurons, filename):
    """
    Creates new file with training neurons
    """
    with open(filename, "wb") as f:
        f.write(struct.pack('<i', len(training_neurons)))
    
        for neuron in training_neurons:
            f.write(struct.pack('<i', neuron))
            

if __name__ == "__main__":
    #filename = "/home/eugene/Output/networks/MLong/generateNetwork/dendritic_tree/parallelChains_20_5500hvcI/training_neurons.bin"
    #
    filename = "/home/eugene/Output/networks/chainGrowth/network2000/training_neurons_20.bin"
    #filename = "/home/eugene/results/immature/clusters/5/training_neurons.bin"
    
    N_RA = 2000
    N_TR = 20    
    
    #training_neurons = range(N_TR)
    training_neurons = np.random.choice(N_RA, N_TR, replace=False)
    
    #print "Training neurons: ",training_neurons
    write_training_neurons(training_neurons, filename)    
    
    #training_neurons_read = reading.read_training_neurons(filename)
    
    
    #print training_neurons_read

