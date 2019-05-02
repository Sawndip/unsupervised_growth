#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 01:07:17 2018

@author: jingroup

Script plots results of inhibitory and excitatory inputs test for modeled HVC
"""
import reading
import matplotlib.pyplot as plt
import numpy as np
import os

dirname = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/inhAndExcTest/"

filename = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/inhAndExcTest/Ginh_0.00_Gee_20.00.bin"
fileTraining = "/home/eugene/Output/networks/chainGrowth/test/training_neurons_random.bin"

          
            

#print probability_for_relevant_spike
#print mean_relevant_spike_time

#print mean_spike_time
#print std_spike_time

training = reading.read_training_neurons(fileTraining)

#print probability_for_nonrelevant_spike

nbins = 50

plt.figure()

files = os.listdir(dirname)

indicatorString = "Ginh_3.00"

relevant_files = [f for f in files if indicatorString in f]
print relevant_files
Gee = []

for f in relevant_files:
    Gee.append(float((f.split(".bin")[0]).split("_")[-1]))
    
Gee.sort()
print Gee

sorted_files = [indicatorString + "_Gee_" + str(G) + "0.bin" for G in Gee]


for i, f in enumerate(sorted_files):
    probability_for_relevant_spike, mean_relevant_spike_time, std_relevant_spike_time, \
        probability_for_nonrelevant_spike, mean_spike_time, std_spike_time = reading.read_inhAndExc_test(os.path.join(dirname, f))
  
    label = "Gee = " + str(Gee[i])
    plt.hist(np.delete(probability_for_relevant_spike, training), label=label, bins=50)
#plt.hist(np.delete(probability_for_nonrelevant_spike, training), label='spontaneous activity', bins=50)

plt.legend()
plt.xlabel('probability to spike')
plt.ylabel('# of neurons')
plt.xlim([-0.1, 1.1])

plt.show()