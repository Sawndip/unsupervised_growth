# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 10:54:21 2018

@author: jingroup
"""
import reading
import numpy as np
import os


def findNeuronalReplacementTimes(dataDir):
    """
    Function finds neuronal replacement times
    """

    files = os.listdir(dataDir)
    
    trials = []
    for f in files:
        if "replacement_history" in f:
            trials.append(f.split("replacement_history_"))
    
    
    print trials
    #currentTrial = 1000 
    
    
    
    #filename = dataDir + "replacement_history_" + str(currentTrial) + ".bin"
    
    #if  os.path.getsize(filename) > 0:        
      #  N_RA, _, time_from_previous_replacement = reading.read_replacement_history(filename)
        
        #print time_from_previous_replacement


if __name__ == "__main__":
    dataDir = "/mnt/hodgkin/eugene/results/immature/clusters/matTrans62/"
    findNeuronalReplacementTimes(dataDir)