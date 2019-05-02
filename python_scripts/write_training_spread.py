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
            

if __name__ == "__main__":
    
    filename = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/noImmatureOut2/training_spread.bin"
    
    spread_times = [2.24469, -1.90563, 1.19244, -2.28566, 1.38835, -1.51658, -0.684362, -1.67624, -2.41211, -2.10426]
    
    with open(filename, "wb") as f:
        f.write(struct.pack('<i', len(spread_times)))
        
        
        for t in spread_times:
            f.write(struct.pack('<d', t))
    
    
Neuron 146 spread time = -2.24469
Neuron 195 spread time = -1.90563
Neuron 228 spread time = 1.19244
Neuron 252 spread time = -2.28566
Neuron 461 spread time = 1.38835
Neuron 525 spread time = -1.51658
Neuron 726 spread time = -0.684362
Neuron 860 spread time = -1.67624
Neuron 873 spread time = -2.41211
Neuron 893 spread time = -2.10426

