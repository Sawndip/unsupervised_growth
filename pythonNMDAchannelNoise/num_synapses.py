# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 15:47:33 2017

@author: jingroup

Script checks number of active and supersynapses
"""

import matplotlib.pyplot as plt
import reading

filename = "/home/eugene/Output/num_synapses.bin"

trial_number, num_active_synapses, num_supersynapses = reading.read_num_synapses(filename)
        
print trial_number
print num_active_synapses
print num_supersynapses

