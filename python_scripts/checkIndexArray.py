# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 15:44:53 2017

@author: jingroup

Script checks global index array
"""
import reading

filename = "/home/eugene/hodgkinData/gabaMaturation300117/membraneTraces/global_index_array.bin"

(N_RA, Id_RA_global) = reading.read_global_index_array(filename)

print N_RA
print Id_RA_global 

