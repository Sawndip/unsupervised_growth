# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:12:47 2017

@author: jingroup

Script reads replacement history
"""

import reading
import numpy as np

directory = "/mnt/hodgkin/eugene/results/immature/clusters/matTrans62/"
start_number = 52000 # starting point

filename = directory + "replacement_history_" + str(start_number) + ".bin"

time_from_previous_replacement = reading.read_replacement_history(filename)

coord_HVCRA = reading.read_coordinates(directory + "RA_xy_" + str(start_number) + ".bin")


print len(time_from_previous_replacement)
print time_from_previous_replacement


print time_from_previous_replacement[2]

print time_from_previous_replacement[2][1334]
print time_from_previous_replacement[2][234]


print coord_HVCRA[1334]
print coord_HVCRA[234]

print coord_HVCRA[234] == coord_HVCRA[1334]

#for i, t in enumerate(time_from_previous_replacement[2]):
#    if t < start_number:
#        print "Neuron {0} with last replacement trial {1}".format(i, t)

