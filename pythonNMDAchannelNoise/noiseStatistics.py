# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 10:53:03 2016

@author: jingroup
"""

"""
Script calculates noise statistics of the 2 compartment RA neuron in HVC
"""

import reading
import numpy as np
import math

filename = "/home/eugene/Output/RAneurons/RA97.bin"

TRIAL_START = 200
TRIAL_END = 800

def std(x):
    mean = np.mean(x)

    s = 0

    for e in x:
        s = s + (e - mean) ** 2
    
    return math.sqrt(s / len(x))

(t, Vs, Is, n, h, Vd, Id, r, c, Ca, Gexc_d, Ginh_d, Gexc_s, Ginh_s, s_s, s_d, Ei, flag, Nsoma, Ndend) = reading.read_hh2(filename)


vs = [v for ind, v in enumerate(Vs) if t[ind] > TRIAL_START and t[ind] < TRIAL_END]
vd = [v for ind, v in enumerate(Vd) if t[ind] > TRIAL_START and t[ind] < TRIAL_END]

#print vs


sigma_vs = std(vs)
sigma_vd = std(vd)
sigma_sd = std(s_d)

print "Mean fraction of open NMDA channels in dendritic compartment: ", np.mean(s_d)
print "Standard deviation in fraction of open NMDA channels in dendritic compartment: ", sigma_sd

print "Mean Vs: ", np.mean(Vs), "mV"
print "Standard deviation in Vs: ", sigma_vs, "mV"

print "Mean Vd: ", np.mean(Vd), "mV"
print "Standard deviation in Vd: ", sigma_vd, "mV"