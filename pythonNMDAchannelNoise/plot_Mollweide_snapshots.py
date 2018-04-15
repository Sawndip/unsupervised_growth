# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 12:31:07 2018

@author: jingroup

Script plots location of fired neurons using Mollweide projection of sphere
"""
import matplotlib.pyplot as plt
import math
import numpy as np
import utils
import reading
import os

R = 1.0 # radius of sphere

num_iter = 10000
tolerance = 1e-7

window = 10
num_rows = 3
num_cols = 3

#dirname = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/noImmatureOut8/"
dirname = "/home/eugene/results/immature/clusters/11/"

trial_number = 13600

fileCoordRA = os.path.join(dirname, "RA_xy_" + str(trial_number) + ".bin")
fileTraining = os.path.join(dirname, "training_neurons.bin")
#fileSoma = os.path.join(dirname, "spike_times_soma_" + str(trial_number) + ".bin")
#fileSoma = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/test/noImmatureOut8/test_spike_times_soma_5.bin"
fileSoma = "/home/eugene/results/immature/clusters/test/11/test_spike_times_soma_5.bin"

coord = reading.read_coordinates(fileCoordRA)
training_neurons = reading.read_training_neurons(fileTraining)

(_, _, spike_times_soma, neuron_fired_soma) = reading.read_time_info(fileSoma)

print "Number of spiked neurons: ",len(spike_times_soma)

for spikes, neuron_id in zip(spike_times_soma, neuron_fired_soma):
    if len(spikes) > 6:
        print "Neuron {0} produced {1} spikes".format(neuron_id[0], len(spikes))

#print "Dedritic spikes: ", neuron_fired_dend
#print "Dendritic spike times: ", spike_times_dend

#print "Somatic spikes: ", neuron_fired_soma
#print "Somatic spike times: ", spike_times_soma

#print ordered_soma_spikes_raw
#print ordered_soma_raw

first_spike_times = [spikes[0] for spikes in list(spike_times_soma)]
id_fired = [ids[0] for ids in list(neuron_fired_soma)]

ordered_first_spike_times, ordered_id_fired = zip(*sorted(zip(first_spike_times, id_fired)))

### find first spike of one of the training neurons
start_id = -1

for i, neuron_id in enumerate(ordered_id_fired):
    if neuron_id in training_neurons:
        start_id = i
        break
    
if start_id < 0:
    print "No training neurons fired!"
else:
    
    # calculate meridians
    num_meridians = 10
    longitude_meridians = np.linspace(-np.pi, np.pi, num_meridians)
    meridians = []
    
    for i in range(num_meridians):    
        meridians.append(utils.get_meridian(longitude_meridians[i], 1000, tolerance, num_iter))
        
    # plot parallels
    num_parallels = 10
    latitude_parallels = np.linspace(-np.pi/2.0, np.pi/2.0, num_parallels)
    parallels = []
    
    for i in range(num_parallels):    
        parallels.append(utils.get_parallel(latitude_parallels[i], 1000, tolerance, num_iter))
        
        
        
    
    start_time = ordered_first_spike_times[start_id]
    current_id = start_id
    
    f, axarr = plt.subplots(num_rows, num_cols)
    
    for i in range(num_rows):
        for j in range(num_cols):
            neurons_to_plot = []            
            
            while current_id < len(ordered_id_fired):
                if ordered_first_spike_times[current_id] < start_time + window:
                    neurons_to_plot.append(ordered_id_fired[current_id])
                    current_id += 1
                else:
                    break
            
            if len(neurons_to_plot) > 0:
                ### calculate latitude and longitude of active HVC-RA neurons ####
                longitude, latitude = utils.calculate_longAndLat(coord[neurons_to_plot])
                
                print min(longitude) / np.pi
                print max(longitude) / np.pi
                print min(latitude) / (np.pi/2.)
                print max(latitude) / (np.pi/2.)
                
                Mollweide_RA = np.empty(shape=(len(neurons_to_plot),2), dtype=np.float32) # array with Mollweide coordinates of active HVC(RA) neurons
                Mollweide_RA.fill(np.nan)
                
                
                for k in range(len(neurons_to_plot)):
                    Mollweide_RA[k,0], Mollweide_RA[k,1] = utils.Mollweide_projection(longitude[k], latitude[k], tolerance, num_iter)
                
                
                colors = []
                s = []
                
                
                
               
                axarr[i][j].scatter(Mollweide_RA[:,0], Mollweide_RA[:,1])
                
    
            for k in range(num_parallels):
                axarr[i][j].plot(parallels[k][0], parallels[k][1], c='k')
          
            for k in range(num_meridians):
                axarr[i][j].plot(meridians[k][0], meridians[k][1], c='k')
            
            axarr[i][j].set_xticks([])
            axarr[i][j].set_yticks([])
            axarr[i][j].set_title("{0} - {1} ms".format(start_time, start_time + window))
            
            start_time += window
          
    #ax.scatter(coord_I[:,0], coord_I[:,1], coord_I[:,2], c='r')
    #ax.scatter(coord_RA[:,0], coord_RA[:,1], coord_RA[:,2], c='b')
    
    plt.show()

