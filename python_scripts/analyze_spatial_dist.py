# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 14:35:53 2017

@author: jingroup

Script analyzes spatial distributions
"""
import space
import reading
import math
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import os

def get_fraction_chain_neurons(dirname):
    """
    Get fraction of storngly connected neurons in the chain
    """
    RARA = os.path.join(dirname, "RA_RA_super_connections.bin")
    
    strongly_connected = set()
    
    (N_RA, RA_super_targets, _) = reading.read_connections(RARA)

    
    for i in xrange(N_RA):
        if len(RA_super_targets[i]) > 0:
            strongly_connected.update(RA_super_targets[i])
            
    return len(strongly_connected) / float(N_RA)
                

def find_latest_checkpoint(dirname):
    """
    Find trial number of the latest checkpoint created in simulation
    """
    files = os.listdir(dirname)

    latest_checkpoint = -1

    for f in files:
        if "replacement_history_" in f:
            trial_num = int(f.split("_")[2])
            
            #print f.split("_")[2]
            if trial_num > latest_checkpoint:
                latest_checkpoint = trial_num
                
            
    return latest_checkpoint
    
def get_statistics_of_grown_conn(arrangement, dirname):
    """
    Get mean and std of grown connections between HVC(RA) neurons
    
    Input: arrangement: spatial arrangement of neurons: "square", "sphere", "cube"
           dirname: directory with data
    """
    latest_checkpoint = find_latest_checkpoint(dirname) 

    print "latest_checkpoint for directory {0} = {1}".format(dirname, latest_checkpoint)
    
    
    RARA = os.path.join(dirname, "RA_RA_super_connections_" + str(latest_checkpoint) + "_.bin")
    RA_xy = os.path.join(dirname, "RA_xy_" + str(latest_checkpoint) + "_.bin")
    I2RA = os.path.join(dirname, "I_RA_connections_" + str(latest_checkpoint) + "_.bin")
    
    (N_I, _, I_targets_G) = reading.read_connections(I2RA)
    
    Gie = -1.0    
    
    for i in xrange(N_I):
        if len(I_targets_G[i]) > 0:
            Gie = I_targets_G[i][0]
            break
    
       
    
    coord = reading.read_coordinates(space.num_coordinates[arrangement], RA_xy) # coordinates of HVC(RA) neurons

    dist_RA2RA = [] # array with distances between supersynapses    

    (N_RA, RA_super_targets, _) = reading.read_connections(RARA)

   
    for i in xrange(N_RA):
        if len(RA_super_targets[i]) > 0:
            for target_ind in RA_super_targets[i]:
                dist_RA2RA.append(space.distance_function[arrangement](coord[i], coord[target_ind]))
    
    #mean  = np.mean(dist_RA2RA)
    #std = np.std(dist_RA2RA)
    
    return dist_RA2RA, Gie


if __name__ == "__main__":
    dirname = "/home/eugene/results/delays/5ms/sphere/"
    arrangement = "sphere"
    #plot_spatial_dist_fixed(dirname)
    directories = os.listdir(dirname)    
    
    Gie = []    # strength of inhibitory connection
    fraction_chain_neurons = [] # fraction of strongly connected neurons in the chain
    
    distances_between_chain = [] # distances between chain HVC(RA) neurons
    mean = [] # mean distance between HVC(RA) neurons
    std = [] # standard deviation of the distance between HVC(RA) neurons
    
    dirs_to_exclude = ["matureTest", "figures", "180717_lionx_5", "180717_lionx_6", "180717_lionx_7"]    
    
    for dir in directories:    
        if dir not in dirs_to_exclude:
            
            #m, s, gie = get_statistics_of_grown_conn(os.path.join(dirname, dir), networkDir)
            
            dist, gie = get_statistics_of_grown_conn(arrangement, os.path.join(dirname, dir))
            
            
            if gie < 0:
                print "Error! No inhibitory connections are present!"
            else:
                #mean.append(m)
                #std.append(s)
                
                Gie.append(gie)
                distances_between_chain.append(dist)            
                fraction_chain_neurons.append(get_fraction_chain_neurons(os.path.join(dirname, dir)))
        
    #Gie, mean, std, fraction_chain_neurons = zip(*sorted(zip(Gie, mean, std, fraction_chain_neurons)))
    Gie, distances_between_chain, fraction_chain_neurons = zip(*sorted(zip(Gie, distances_between_chain, fraction_chain_neurons)))
    
    #print distances_between_chain
    #print mean
    #print std    
    
    f = plt.figure()
    
    ax = f.add_subplot(111)
    #for i in range(len(Gie)):
    ax.boxplot(distances_between_chain, labels=Gie)
    #ax.errorbar(Gie, mean,yerr=std, fmt='o')
    #ax.set_xlim([0, 0.75])
    ax.set_ylabel("Distance between strongly connected HVC(RA)")
    ax.set_xlabel("$Gie (mS/cm^2)$")
    
    f = plt.figure()
    
    ax = f.add_subplot(111)
    
    ax.scatter(Gie, fraction_chain_neurons)

    ax.set_ylim([0, max(fraction_chain_neurons) + 0.1])
    
    ax.set_xlim([0, max(Gie) + 0.05])
    ax.set_ylabel("fraction of neurons in the chain")
    ax.set_xlabel("$Gie (mS/cm^2)$")
    
    plt.show()
