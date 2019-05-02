# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 14:02:27 2016

@author: jingroup

Script runs simulations to measure trial execution time for different number of
processes; network_update_frequency; diff number of neurons
"""

import subprocess

dirname = '/home/eugene/Output/timing/netUpdateScaling_2/'

#np_max = np.arange(0, 24, 1)
N_RA = 300
N_I = 100
network_update_frequency = [0.02, 0.04, 0.08, 0.1, 0.16, 0.32, 0.64, 1.0, 1.28]

np = 24

for freq in network_update_frequency:
    print "freq = ",freq
    name = 'freq' +  str(freq) + '.bin'
    params = ['mpirun', '-np', str(np), '/home/eugene/projects/chaingrowth/ParallelNMDAchannelNoise/timeIt/out', 
            str(N_RA), str(N_I), str(freq), name, dirname]
    subprocess.call(params)

