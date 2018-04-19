# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 11:32:11 2016

@author: jingroup
"""
import numpy as np
import reading
import os

trial_number = 10500

#dirname = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/maturationTransition2/"
dirname = "/home/eugene/results/immature/clusters/matTrans15/"

fileMature = os.path.join(dirname, "mature_" + str(trial_number) + ".bin")

(_, _, mature_indicators) = reading.read_mature_indicators(fileMature)


mature_neurons = np.where(mature_indicators == 1)[0]
print "Mature neurons: ", list(mature_neurons)