# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:12:47 2017

@author: jingroup

Script reads replacement history
"""

import reading
import numpy as np

directory = "/mnt/hodgkin_home/eugene/Output/networks/networkTest/"
start_number = 0 # starting point

filename = directory + "replacement_history_" + str(start_number) + "_.bin"

time_from_previous_replacement = reading.read_replacement_history(filename)

print len(time_from_previous_replacement)
print time_from_previous_replacement


