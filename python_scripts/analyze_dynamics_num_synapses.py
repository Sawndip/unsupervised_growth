#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 17:20:11 2018

@author: jingroup

Script analyzes dynamics of active and super synapses 
"""
import reading
import os
import matplotlib.pyplot as plt

dirname = "/home/eugene/results/immature/clusters/matTrans55/"

trials = [2600, 5000, 10000, 13400, 17000, 21200] # matTrans55
#trials = [2600, 5000, 8800, 12200, 15000, 19000] # matTrans56
#trials = [3000, 5000, 9200, 12600, 16000, 20400] # matTrans57

num_active = []
num_super = []

for trial in trials:
    fileActiveSynapses = os.path.join(dirname, "RA_RA_active_connections_" + str(trial) + ".bin")
    fileSuperSynapses = os.path.join(dirname, "RA_RA_super_connections_" + str(trial) + ".bin")
    
    (N, _, active_synapses) = reading.read_synapses(fileActiveSynapses)
    (_, _, super_synapses) = reading.read_synapses(fileSuperSynapses)

    total_a = 0
    total_s = 0

    for i in range(N):
        total_a += len(active_synapses[i])
        total_s += len(super_synapses[i])
        
    num_active.append(total_a)
    num_super.append(total_s)


plt.figure()
plt.plot(trials, num_active, '-bo', label='active')
plt.plot(trials, num_super, '-go', label='super')

plt.ylabel('# of synapses')
plt.xlabel('time (trials)')
plt.legend()
plt.show()