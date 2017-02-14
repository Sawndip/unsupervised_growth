# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 15:47:33 2017

@author: jingroup

Script checks statistics of synaptic weight distribution
"""

import matplotlib.pyplot as plt
import reading

filename = "/home/eugene/Output/networks/gabaMaturation100217/weight_statistics.bin"

trial_number, mean, std = reading.read_weight_statistics(filename)
        
print trial_number
print mean
print std


f = plt.figure()

ax1 = f.add_subplot(211)
ax1.plot(trial_number, mean)
ax1.set_ylabel("mean synaptic weight")

ax2 = f.add_subplot(212)
ax2.plot(trial_number, std)
ax2.set_xlabel("time (# trials)")
ax2.set_ylabel("std of synaptic weight")

