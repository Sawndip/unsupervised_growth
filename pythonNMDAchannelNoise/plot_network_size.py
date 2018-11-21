# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 21:57:54 2018

@author: jingroup

Script plots network size vs trial number
"""
import matplotlib.pyplot as plt
import reading
import os
import numpy as np

def get_networkSize_vs_time(dirname):
    """
    Get network size and time from files in the directory
    """
    files = os.listdir(dirname)
    
    time = []
    network_size = []
    
    
    
    for f in files:
        if "RA_RA_super_connections" in f:    
            time.append(int(f.split("_")[-1].split(".bin")[0]))
            
            (N_RA, _, super_synapses) = reading.read_synapses(os.path.join(dirname, f))
            
            has_input_supersynapse = np.zeros(N_RA, np.bool)
    
            #print has_input_supersynapse        
            
            for i in range(N_RA):
                for target in super_synapses[i]:
                    has_input_supersynapse[target] = True
                    
            network_size.append(len(np.where(has_input_supersynapse == True)[0]))
            
    return time, network_size



dataDir = "/mnt/hodgkin/eugene/results/immature/clusters/full/matTrans62"
dataDir_2 = "/mnt/hodgkin/eugene/results/immature/clusters/full/matTrans63"


time, network_size = get_networkSize_vs_time(dataDir)
time, network_size =zip(*sorted( zip(time, network_size)))

time_2, network_size_2 = get_networkSize_vs_time(dataDir_2)
time_2, network_size_2 =zip(*sorted( zip(time_2, network_size_2)))

plt.figure()
plt.plot(time, network_size, label='neurogenesis')
plt.plot(time_2, network_size_2, label='no neurogenesis')
plt.xlabel('Time (# trials)')
plt.ylabel('Network size')
plt.legend()
plt.show()
