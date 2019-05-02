# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 22:01:49 2018

@author: jingroup

Script analyzes inhibitory inputs to population of HVC-RA neurons
"""
import matplotlib.pyplot as plt
import reading
import os
import numpy as np

dirname = "/home/eugene/Output/networks/chainGrowth/matureTest/RA"
fileMature = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays3/mature_5300.bin"

files = os.listdir(dirname)

(_, _, mature_indicators) = reading.read_mature_indicators(fileMature)

mature_neurons = set(np.where(mature_indicators == 1)[0])

counter = 0
counter_mature = 0

for f in files:
    t, Vs, Vd, Gexc_d, Ginh_d, n, h, r, c, Ca = reading.read_hh2_buffer_full(os.path.join(dirname, f))
    
    if counter == 0:
        Ginh = np.copy(Ginh_d)
    else:
        Ginh += Ginh_d
 
    neuron_id = int(f.split(".bin")[0][2:])
    
    print neuron_id
     
    if neuron_id in mature_neurons:
        if counter_mature == 0:
            Ginh_mature = np.copy(Ginh_d)
        else:
            Ginh_mature += Ginh_d
            
        counter_mature += 1
    
    counter += 1
    
        
Ginh = Ginh / float(len(files))
Ginh_mature = Ginh_mature / float(counter_mature) 
    
tmin = 0
tmax = 500


# conductances
f = plt.figure()

ax1 = f.add_subplot(111)
ax1.plot(t, Ginh, label='all neurons')
ax1.plot(t, Ginh_mature, label='mature neurons')
ax1.set_ylabel("< Ginh > (mS/cm^2)")
ax1.set_xlim([tmin, tmax])

ax1.set_xlabel("Time (ms)")
ax1.legend()

plt.show()
