# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 11:32:11 2016

@author: jingroup
"""

import reading

fileMature = "/home/eugene/Output/mature13.bin"

mature = reading.read_mature(fileMature)

print "mature = ", mature

mature_id = [i for i,e in enumerate(mature) if e == 1]

print "Mature neurons: ", mature_id