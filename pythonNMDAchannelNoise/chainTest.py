# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 16:45:56 2016

@author: jingroup

Script reads chain test file
"""

import reading

filename = "/home/eugene/Output/chain_test.bin"

N_RA, num_trials, num_dend_spikes, mean_burst_time, std_burst_time = reading.read_chain_test(filename)

print "N_RA = ",N_RA
print "num_trials = ",num_trials
print "num_dend_spikes = ",num_dend_spikes

print "mean_burst_time = ",mean_burst_time
print "std_burst_time = ",std_burst_time
