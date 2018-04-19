#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 15:39:22 2018

@author: jingroup

Script read maturation properties from file
"""
import reading
import numpy as np

filename = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/maturationTransition3/maturation_properties_150.bin"


(N, trial_number, mature_indicators, maturation_rate, Erest, GCa) = reading.read_maturation_properties(filename)


mature_neurons = np.where(mature_indicators == 1)
immature_neurons =  np.where(mature_indicators == 0)

print list(mature_neurons[0])
print N
print trial_number
print list(mature_indicators[mature_neurons])
print list(maturation_rate[mature_neurons])
#print list(maturation_rate[immature_neurons])

print list(Erest[mature_neurons])
#print list(Erest[immature_neurons])
print np.min(Erest[immature_neurons])

print np.argmin(Erest[immature_neurons])


print list(GCa[mature_neurons])
print np.max(GCa[immature_neurons])
print np.argmax(GCa[immature_neurons])
