# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 14:35:53 2017

@author: jingroup

Script analyzes spatial distributions
"""
import reading
import math
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import os

SQUARE_SIDE = 100
SIDE = SQUARE_SIDE * math.sqrt(2)


def d(x1, y1, x2, y2):
    """ Computes distance between two points in space"""
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
        
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
                
    
def get_statistics_of_grown_conn(dirname, networkDir):
    """
    Get mean and std of grown connections between HVC(RA) neurons
    """
    RARA = os.path.join(dirname, "RA_RA_super_connections.bin")
    RA_xy = os.path.join(networkDir, "RA_xy.bin")
    I2RA = os.path.join(dirname, "I_RA_connections.bin")
    
    (N_I, _, I_targets_G) = reading.read_connections(I2RA)
    
    Gie = -1.0    
    
    for i in xrange(N_I):
        if len(I_targets_G[i]) > 0:
            Gie = I_targets_G[i][0]
            break
    
       
    
    (x_RA, y_RA) = reading.read_coordinates(RA_xy) # coordinates of HVC(RA) neurons

    dist_RA2RA = [] # array with distances between supersynapses    

    (N_RA, RA_super_targets, _) = reading.read_connections(RARA)

    
    for i in xrange(N_RA):
        if len(RA_super_targets[i]) > 0:
            for target_ind in RA_super_targets[i]:
                dist_RA2RA.append(d(x_RA[i], y_RA[i], x_RA[target_ind], y_RA[target_ind])/SIDE)
    
    #mean  = np.mean(dist_RA2RA)
    #std = np.std(dist_RA2RA)
    
    return dist_RA2RA, Gie

def plot_spatial_dist_fixed(dirname, networkDir):
    """
    Plot spatial distributions of fixed connections between HVC(RA) and HVC(I) neurons
    using files in directory dirname
    """

    RA_xy = os.path.join(networkDir, "RA_xy.bin")
    I_xy = os.path.join(networkDir, "I_xy.bin")
    RA2I = os.path.join(dirname, "RA_I_connections.bin")
    I2RA = os.path.join(dirname, "I_RA_connections.bin")

    (x_I, y_I) = reading.read_coordinates(I_xy)
    (x_RA, y_RA) = reading.read_coordinates(RA_xy)

    (N_RA, RA_targets, RA_targets_G) = reading.read_connections(RA2I)
    (N_I, I_targets, I_targets_G) = reading.read_connections(I2RA)

    # Here we calculate mean distance between RA and their I targets
    
    dist_RA2I = []
    
    for i in xrange(N_RA):
        if len(RA_targets[i]) > 1:
            for target_ind in RA_targets[i]:
                dist_RA2I.append(d(x_RA[i], y_RA[i], x_I[target_ind], y_I[target_ind])/SIDE)
    
    
    # Here we calculate mean distance between I and their RA targets
    
    dist_I2RA = []
    
    for i in xrange(N_I):
        if len(I_targets[i]) > 1:
            for target_ind in I_targets[i]:
                dist_I2RA.append(d(x_I[i], y_I[i], x_RA[target_ind], y_RA[target_ind])/SIDE)
                
    numBins = 30
    
    f1 = plt.figure()
    ax = f1.add_subplot(111)
    ax.hist(dist_RA2I, numBins)
    ax.set_xlabel("Normalized distance")
    ax.set_ylabel("Number of connections")
    ax.set_title("Spatial distribution of connections from RA onto I neurons")
    ax.set_xlim([0,1])
    
    numBins = 25
    
    f2 = plt.figure()
    ax = f2.add_subplot(111)
    ax.hist(dist_I2RA, numBins)
    ax.set_xlabel("Normalized distance")
    ax.set_ylabel("Number of connections")
    ax.set_title("Spatial distribution of connections from I onto RA neurons")
    ax.set_xlim([0,1])
    
    plt.show()


if __name__ == "__main__":
    dirname = "/mnt/hodgkin/eugene/lionX/clustered/r0.2/"
    networkDir = "/mnt/hodgkin/eugene/lionX/clustered/network/"
    #plot_spatial_dist_fixed(dirname)
    directories = os.listdir(dirname)    
    
    Gie = []    # strength of inhibitory connection
    fraction_chain_neurons = [] # fraction of strongly connected neurons in the chain
    
    distances_between_chain = [] # distances between chain HVC(RA) neurons
    mean = [] # mean distance between HVC(RA) neurons
    std = [] # standard deviation of the distance between HVC(RA) neurons
    
    for dir in directories:    
        
        
        #m, s, gie = get_statistics_of_grown_conn(os.path.join(dirname, dir), networkDir)
        
        dist, gie = get_statistics_of_grown_conn(os.path.join(dirname, dir), networkDir)
        
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
    ax.set_xlim([0, 0.75])
    ax.set_ylabel("fraction of neurons in the chain")
    ax.set_xlabel("$Gie (mS/cm^2)$")
    
    plt.show()
