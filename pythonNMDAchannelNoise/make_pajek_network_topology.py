# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 10:38:19 2018

@author: jingroup

Script creates pajek .net file for network topology
"""
import reading
import os
import utils
import numpy as np

def write_allCoords(training_neurons, coord_HVCRA, coord_HVCI, filename):
    """
    Write coordinates of HVC(RA) and HVC(I) neurons to file
    
    Input: 
            training_neurons - list with ids of training_neurons
            coord_HVCRA - numpy array with coordinates of HVC(RA) neurons
            coord_HVCRI - numpy array with coordinates of HVC(I) neurons
            f - filestream to which data is written
    """
    N_RA = coord_HVCRA.shape[0] # number of HVC(RA) neurons            
    N_I = coord_HVCI.shape[0] # number of HVC(I) neurons            
        
    with open(filename, 'w') as f:
        f.write("*Vertices {0}\n".format(N_RA+N_I))
         
       
        for i in range(N_RA):
            if i in training_neurons:    
                f.write('{0} "{1}" {2} {3} {4} ic Green\n'.format(i+1, i, coord_HVCRA[i][0], coord_HVCRA[i][1], coord_HVCRA[i][2]))
            else:    
                f.write('{0} "{1}" {2} {3} {4} ic Red\n'.format(i+1, i, coord_HVCRA[i][0], coord_HVCRA[i][1], coord_HVCRA[i][2]))
    
        for i in range(N_RA, N_RA+N_I):
            f.write('{0} "{1}" {2} {3} {4} ic Blue\n'.format(i+1, i, coord_HVCI[i-N_RA][0], coord_HVCI[i-N_RA][1], coord_HVCI[i-N_RA][2]))


def write_pajek_network_topology(dirname):
    """
    Create .net file with locations and connections between HVC-RA and HVC-I
    neurons
    """
    file_RA_xy = os.path.join(dirname, "RA_xy.bin")
    file_I_xy = os.path.join(dirname, "I_xy.bin")
    
    RA2I = os.path.join(dirname, "RA_I_connections.bin")
    I2RA = os.path.join(dirname, "I_RA_connections.bin")
    
    file_training = os.path.join(dirname, "training_neurons_clustered.bin")
    file_pajek = os.path.join(dirname, "network_topology_clustered.net")
    
    (N_RA, RA_targets, RA_targets_G, _, _) = reading.read_connections(RA2I)
    (N_I, I_targets, I_targets_G, _, _) = reading.read_connections(I2RA)
    
    coord_RA = reading.read_coordinates(file_RA_xy)
    coord_I = reading.read_coordinates(file_I_xy)
    
    #print targets_ID
    #print targets_G
    if file_training:
        training_neurons = reading.read_training_neurons(file_training)
    else:
        training_neurons = []
        
    print "Training neurons: ",training_neurons
    
    with open(file_pajek, 'w') as f:
        f.write("*Vertices {0}\n".format(N_RA+N_I))
        
                
                
        for i in range(N_RA):
            if i in training_neurons:    
                f.write('{0} "{1}" {2} {3} {4} ic Green\n'.format(i+1, i, coord_RA[i][0], coord_RA[i][1], coord_RA[i][2]))
            else:    
                f.write('{0} "{1}" {2} {3} {4} ic Yellow\n'.format(i+1, i, coord_RA[i][0], coord_RA[i][1], coord_RA[i][2]))
    
        for i in range(N_RA, N_RA+N_I):
            f.write('{0} "{1}" {2} {3} {4} ic Red\n'.format(i+1, i, coord_I[i-N_RA][0], coord_I[i-N_RA][1], coord_I[i-N_RA][2]))
    
        f.write("*Arcs\n".format(N_RA))
        
        # write targets of HVC(RA) neurons
        for i, targets in enumerate(RA_targets):
            for j, target in enumerate(targets):
                f.write('{0} {1} {2} c Green\n'.format(i+1, N_RA+target+1, RA_targets_G[i][j]))
                
         # write targets of HVC(I) neurons
        for i, targets in enumerate(I_targets):
            for j, target in enumerate(targets):
                f.write('{0} {1} {2} c Red\n'.format(N_RA+i+1, target+1, I_targets_G[i][j]))

def write_pajek_hvcRA_coord(dirname, trial_number):
    """
    Create .net file with locations of HVC-RA neurons in array.
    Mature neurons are highlighted
    """        
    file_RA_xy = os.path.join(dirname, "RA_xy_" + str(trial_number) + ".bin")
    
    file_training = os.path.join(dirname, "training_neurons.bin")
    file_pajek = os.path.join(dirname, "network_" + str(trial_number) + ".net")
    fileMature = os.path.join(dirname, "mature_" + str(trial_number) + ".bin")
    
    coord_RA = reading.read_coordinates(file_RA_xy)
    training_neurons = reading.read_training_neurons(file_training)
    (_, _, mature_indicators) = reading.read_mature_indicators(fileMature)
    
    mature_neurons = np.where(mature_indicators == 1)[0]   
    num_neurons = coord_RA.shape[0]
    # sort array with neurons and training neurons #
    training_neurons.sort()
    mature_neurons.sort()
    
    with open(file_pajek, 'w') as f:
        f.write("*Vertices {0}\n".format(num_neurons))
        
        for i in range(num_neurons):
            if i in training_neurons:    
                f.write('{0} "{1}" {2} {3} {4} ic Green\n'.format(i+1, i, coord_RA[i][0], coord_RA[i][1], coord_RA[i][2]))
            elif i in mature_neurons:  
                f.write('{0} "{1}" {2} {3} {4} ic Black\n'.format(i+1, i, coord_RA[i][0], coord_RA[i][1], coord_RA[i][2]))
            else:
                f.write('{0} "{1}" {2} {3} {4} ic Yellow\n'.format(i+1, i, coord_RA[i][0], coord_RA[i][1], coord_RA[i][2]))

def write_pajek_neurons_connected_by_supersynapses(dirname, trial_number):
    """
    Create .net file with locations and supersynaptic connections for HVC-RA neurons connected by supersynapses
    """
    file_RA_xy = os.path.join(dirname, "RA_xy_" + str(trial_number) + ".bin")
    
    file_training = os.path.join(dirname, "training_neurons.bin")
    file_pajek = os.path.join(dirname, "network_" + str(trial_number) + ".net")
    fileSuperSynapses = os.path.join(dirname, "RA_RA_super_connections_" + str(trial_number) + ".bin")
    fileWeights = os.path.join(dirname, "weights_" + str(trial_number) + ".bin")

    coord_RA = reading.read_coordinates(file_RA_xy)
    training_neurons = reading.read_training_neurons(file_training)
    (N_RA, _, super_synapses) = reading.read_synapses(fileSuperSynapses)
    (N_RA, _, weights) = reading.read_weights(fileWeights) 

    network_neurons = set(training_neurons)
        

    for i in range(N_RA):
        for target in super_synapses[i]:
            network_neurons.add(target)
    
    network_neurons = sorted(list(network_neurons))
    
    
    num_neurons = len(network_neurons)
    # sort array with neurons and training neurons #
    training_neurons.sort()
    
    with open(file_pajek, 'w') as f:
        f.write("*Vertices {0}\n".format(num_neurons))
        
        for i, neuron_id in enumerate(network_neurons):
            if neuron_id in training_neurons:    
                f.write('{0} "{1}" {2} {3} {4} ic Green\n'.format(i+1, neuron_id, coord_RA[neuron_id][0], coord_RA[neuron_id][1], coord_RA[neuron_id][2]))
            else:    
                f.write('{0} "{1}" {2} {3} {4} ic Yellow\n'.format(i+1, neuron_id, coord_RA[neuron_id][0], coord_RA[neuron_id][1], coord_RA[neuron_id][2]))
        
        
        f.write("*Arcs\n")
        
        # write targets of HVC(RA) neurons
        for i, source_id in enumerate(network_neurons):
            for target_id in super_synapses[source_id]:
                try:
                    ind = utils.index(network_neurons, target_id)                 
                    f.write('{0} {1} {2} c Green\n'.format(i+1, ind+1, weights[source_id][target_id]))
                except ValueError:
                    continue

 
def write_pajek_neurons(dirname, trial_number):
    """
    Create .net file with locations and connections between mature HVC-RA neurons in array
    """
    file_RA_xy = os.path.join(dirname, "RA_xy_" + str(trial_number) + ".bin")
    
    file_training = os.path.join(dirname, "training_neurons.bin")
    file_pajek = os.path.join(dirname, "network_" + str(trial_number) + ".net")
    fileMature = os.path.join(dirname, "mature_" + str(trial_number) + ".bin")
    fileSuperSynapses = os.path.join(dirname, "RA_RA_super_connections_" + str(trial_number) + ".bin")
    fileWeights = os.path.join(dirname, "weights_" + str(trial_number) + ".bin")
    
    coord_RA = reading.read_coordinates(file_RA_xy)
    training_neurons = reading.read_training_neurons(file_training)
    (N_RA, _, weights) = reading.read_weights(fileWeights) 
    (_, _, mature_indicators) = reading.read_mature_indicators(fileMature)
    (_, _, super_synapses) = reading.read_synapses(fileSuperSynapses)

    mature_neurons = np.where(mature_indicators == 1)[0]   
    #print list(mature_neurons)
    #mature_neurons = range(N_RA)
    num_neurons = len(mature_neurons)
    # sort array with neurons and training neurons #
    training_neurons.sort()
    mature_neurons.sort()
    
    with open(file_pajek, 'w') as f:
        f.write("*Vertices {0}\n".format(num_neurons))
        
        for i, neuron_id in enumerate(mature_neurons):
            if neuron_id in training_neurons:    
                f.write('{0} "{1}" {2} {3} {4} ic Green\n'.format(i+1, neuron_id, coord_RA[neuron_id][0], coord_RA[neuron_id][1], coord_RA[neuron_id][2]))
            else:    
                f.write('{0} "{1}" {2} {3} {4} ic Yellow\n'.format(i+1, neuron_id, coord_RA[neuron_id][0], coord_RA[neuron_id][1], coord_RA[neuron_id][2]))
        
        
        f.write("*Arcs\n".format(num_neurons))
        
        # write targets of HVC(RA) neurons
        for i, source_id in enumerate(mature_neurons):
            for target_id in super_synapses[source_id]:
                try:
                    ind = utils.index(mature_neurons, target_id)                 
                    f.write('{0} {1} {2} c Green\n'.format(i+1, ind+1, weights[source_id][target_id]))
                except ValueError:
                    continue

def write_pajek_network_subset(dirname, trial_number, N, fileSpikes):
    """
    Create .net file with locations and connections between mature HVC-RA neurons in array
        first N mature neurons that spiked are plotted
    """
    file_RA_xy = os.path.join(dirname, "RA_xy_" + str(trial_number) + ".bin")
    
    file_training = os.path.join(dirname, "training_neurons.bin")
    file_pajek = os.path.join(dirname, "network_subset_" + str(trial_number) + ".net")
    fileMature = os.path.join(dirname, "mature_" + str(trial_number) + ".bin")
    fileSuperSynapses = os.path.join(dirname, "RA_RA_super_connections_" + str(trial_number) + ".bin")
    fileWeights = os.path.join(dirname, "weights_" + str(trial_number) + ".bin")
    
    coord_RA = reading.read_coordinates(file_RA_xy)
    training_neurons = reading.read_training_neurons(file_training)
    (N_RA, _, weights) = reading.read_weights(fileWeights) 
    (_, _, mature_indicators) = reading.read_mature_indicators(fileMature)
    (_, _, super_synapses) = reading.read_synapses(fileSuperSynapses)

   
    #print list(mature_neurons)
    #mature_neurons = range(N_RA)
   
    # sort array with neurons and training neurons #
    training_neurons.sort()
   
    
   
#fileDend = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/test/noImmatureOut4/test_spike_times_dend_5.bin"
#fileSoma = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/test/noImmatureOut4/test_spike_times_soma_5.bin"



    (_, _, spike_times_soma, neuron_fired_soma) = reading.read_time_info(fileSpikes)

    ordered_soma_spikes_raw, ordered_soma_raw = zip(*sorted(zip(spike_times_soma, neuron_fired_soma)))

    first_mature_spiked = []

    for spikes, neuron_ids in zip(ordered_soma_spikes_raw, ordered_soma_raw):
        if len(first_mature_spiked) >= N:
            break
        
        if mature_indicators[neuron_ids[0]] == 1:
            first_mature_spiked.append(neuron_ids[0])

    first_mature_spiked.sort()
    
    num_neurons = len(first_mature_spiked)
    
    with open(file_pajek, 'w') as f:
        f.write("*Vertices {0}\n".format(num_neurons))
        
        for i, neuron_id in enumerate(first_mature_spiked):
            if neuron_id in training_neurons:    
                f.write('{0} "{1}" {2} {3} {4} ic Green\n'.format(i+1, neuron_id, coord_RA[neuron_id][0], coord_RA[neuron_id][1], coord_RA[neuron_id][2]))
            else:    
                f.write('{0} "{1}" {2} {3} {4} ic Yellow\n'.format(i+1, neuron_id, coord_RA[neuron_id][0], coord_RA[neuron_id][1], coord_RA[neuron_id][2]))
        
        
        f.write("*Arcs\n".format(num_neurons))
        
        # write targets of HVC(RA) neurons
        for i, source_id in enumerate(first_mature_spiked):
            for target_id in super_synapses[source_id]:
                try:
                    ind = utils.index(first_mature_spiked, target_id)                 
                    f.write('{0} {1} {2} c Green\n'.format(i+1, ind+1, weights[source_id][target_id]))
                except ValueError:
                    continue
  
if __name__ == "__main__":
    #dirname = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/maturationTransition4/"

    dirname = "/mnt/hodgkin/eugene/results/immature/clusters/matTrans62/"
    trial_number = 23800
   

    #dirname = "/mnt/hodgkin/eugene/results/immature/clusters/matTrans62/"
    #trial_number = 23800

   
   # fileSpikes = "/home/eugene/results/immature/clusters/test/matTrans29/test_spike_times_soma_10.bin"

    dirname = "/mnt/hodgkin/eugene/Output/networks/chainGrowth/network200RA55I"
    fileTraining = "/mnt/hodgkin/eugene/Output/networks/chainGrowth/network200RA55I/training_neurons_random.bin"
    
    
    training_neurons = set(reading.read_training_neurons(fileTraining))

    
    coord_HVCRA = reading.read_coordinates(os.path.join(dirname, "RA_xy.bin"))
    coord_HVCI = reading.read_coordinates(os.path.join(dirname, "I_xy.bin"))    
    
    #print training_neurons
    #print coord_RA
    
    
    #write_pajek_neurons_connected_by_supersynapses(dirname, trial_number)
    write_allCoords(training_neurons, coord_HVCRA, coord_HVCI, os.path.join(dirname, "pajek.net"))
    #write_pajek_neurons(dirname, trial_number)
    #write_pajek_network_subset(dirname, trial_number, N, fileSpikes)
    #write_pajek_hvcRA_coord(dirname, trial_number)