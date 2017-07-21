# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 19:29:31 2017

@author: jingroup

Script read and writes training neurons 
"""

import reading
import struct


def write_training_neurons(training_neurons, filename):
    """
    Creates new file with training neurons
    """
    with open(filename, "wb") as f:
        f.write(struct.pack('<i', len(training_neurons)))
    
        for neuron in training_neurons:
            f.write(struct.pack('<i', neuron))
            

if __name__ == "__main__":
    
    filename = "/home/eugene/Output/networks/sphere/training_neurons.bin"
    
    training_neurons = [131, 95, 277, 179]
    
    write_training_neurons(training_neurons, filename)    
    
    training_neurons_read = reading.read_training_neurons(filename)
    
    
    print training_neurons_read

