# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 19:19:50 2017

@author: jingroup

Script reads file with average trial execution time information
"""

import reading

filename = "/home/eugene/Output/exec_time_1.bin"

(N_RA, N_I, np, timestep, network_update_frequency, average_execution_time) = reading.read_execution_time(filename)

print "N_RA = ",N_RA
print "N_I = ",N_I
print "np = ",np
print "timestep = ",timestep
print "network_update_frequency = ",network_update_frequency
print "average_execution_time = ",average_execution_time
