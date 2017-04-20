# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 15:47:33 2017

@author: jingroup

Script checks number of active and supersynapses
"""

import matplotlib.pyplot as plt
import reading

filename = "/mnt/hodgkin_home/eugene/lionX/gabaMaturation180417_3/num_synapses.bin"

trial_number, num_active_synapses, num_supersynapses = reading.read_num_synapses(filename)
        
print trial_number
print num_active_synapses
print num_supersynapses


f = plt.figure()

ax1 = f.add_subplot(211)
ax1.plot(trial_number, num_active_synapses)
ax1.set_ylabel("num_active_synapses")

ax2 = f.add_subplot(212)
ax2.plot(trial_number, num_supersynapses)
ax2.set_xlabel("time (# trials)")
ax2.set_ylabel("num_supersynapses")
