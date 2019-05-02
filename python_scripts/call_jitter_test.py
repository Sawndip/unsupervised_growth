# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 09:48:59 2018

@author: jingroup
"""

import subprocess
import os

dirname = '/mnt/hodgkin/eugene/results/immature/clusters/'

np = 8

simnames = ["matTrans68", "matTrans62", "matTrans69", "matTrans66", "matTrans64"]
trials = [37400, 58400, 27800, 71400, 57200]

num_test_trials = 10

for simname, trial in zip(simnames, trials):
    outDir = dirname + "test/" + simname + "/trial" + str(trial) + "/"
    dataDir =  dirname + simname + "/"
    
    params = ['mpirun', '-np', str(np), '/home/eugene/projects/chaingrowth/ParallelNMDAchannelNoise/chainTest/main', 
              dataDir, str(trial), str(num_test_trials), outDir]

    
    print params

    subprocess.call(params)






