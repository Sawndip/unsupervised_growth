#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 21:26:57 2018

@author: jingroup

Script reads replace neurons
"""
import reading

filename = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/maturationTransition3/replaced_neurons.bin"

replaced_neurons = reading.read_replaced_neurons(filename)

print replaced_neurons

