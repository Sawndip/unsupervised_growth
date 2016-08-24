# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 16:46:09 2016

@author: jingroup
"""

"""
Script read time information of running simulation
"""

import reading as read

filename = "/home/eugene/Output/timeInfo.bin"

trial_number, internal_time, network_time = read.read_simTime_info(filename)

print "trial_number = ", trial_number
print "internal_time = ", internal_time
print "network_time = ", network_time

