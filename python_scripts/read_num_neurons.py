#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 12:30:28 2018

@author: jingroup

Script checks number of neuron in a file
"""
import reading

filename = "/home/eugene/results/immature/clusters/matTrans28/num_neurons.bin"

N_RA, N_I = reading.read_num_neurons(filename)

print N_RA, N_I

