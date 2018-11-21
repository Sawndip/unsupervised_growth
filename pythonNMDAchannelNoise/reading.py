# -*- coding: utf-8 -*-
"""
Created on Tue Mar 01 00:09:53 2016

@author: Eugene
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 12:42:25 2016

@author: Eugene
"""

import struct
import numpy as np

SIZE_OF_INT = 4
SIZE_OF_DOUBLE = 8


def read_connections(filename):
    """
    Read connections of neurons from binary file
    """
    with open(filename, 'rb') as file:
        data = file.read()
    
        N = struct.unpack("<i", data[:SIZE_OF_INT])[0]
        
        #print N     
        
        targets_ID = []
        weights = []
        syn_lengths = []
        axonal_delays = []
        
        ind = SIZE_OF_INT    
        
        for i in range(N):
            neurons_ID = struct.unpack("<i", data[ind:(ind+SIZE_OF_INT)])[0]
            number_of_targets = struct.unpack("<i", data[(ind+SIZE_OF_INT):(ind+2*SIZE_OF_INT)])[0]
            
            temp_ID = []
            temp_G = []    
            temp_lengths = []    
            temp_delays = []    
            
            ind += 2*SIZE_OF_INT            
            
            for j in range(number_of_targets):
                ID = struct.unpack("<i", data[ind:(ind + SIZE_OF_INT)])[0]
                #print ID
                G = struct.unpack("<d", data[(ind + SIZE_OF_INT):(ind + SIZE_OF_INT + SIZE_OF_DOUBLE)])[0]
                length = struct.unpack("<d", data[(ind + SIZE_OF_INT + SIZE_OF_DOUBLE):(ind + SIZE_OF_INT + 2*SIZE_OF_DOUBLE)])[0]
                delay = struct.unpack("<d", data[(ind + SIZE_OF_INT + 2*SIZE_OF_DOUBLE):(ind + SIZE_OF_INT + 3*SIZE_OF_DOUBLE)])[0]
                                                
                #print G            
                temp_ID.append(ID)
                temp_G.append(G)
                temp_lengths.append(length)
                temp_delays.append(delay)
            
                ind += SIZE_OF_INT + 3*SIZE_OF_DOUBLE
            
            targets_ID.append(temp_ID)
            weights.append(temp_G)
            syn_lengths.append(temp_lengths)
            axonal_delays.append(temp_delays)
            
            
        return (N, targets_ID, weights, syn_lengths, axonal_delays)
    
def read_coordinates(filename):
    """
    Read coordinates of neurons from binary file
    """
    with open(filename, 'rb') as file:
        data = file.read()
   
        N = struct.unpack('<i', data[:SIZE_OF_INT])[0] # get number of neurons
        dimensionality = struct.unpack('<i', data[SIZE_OF_INT:2*SIZE_OF_INT])[0] # get number of neurons
        model_interneuron_distance = struct.unpack('<d', data[2*SIZE_OF_INT:(2*SIZE_OF_INT+SIZE_OF_DOUBLE)])[0] # get average distance between interneurons in the model
        
       
        ind = 2*SIZE_OF_INT + SIZE_OF_DOUBLE # position in file    
        
        if dimensionality == 2:    
             coord = np.empty(shape=(N,2), dtype=np.float32) # coordinates
             coord.fill(np.nan)
             
             for i in range(N):
                 coord[i][0] = struct.unpack('<d', data[ind:(ind + SIZE_OF_DOUBLE)])[0]
                 coord[i][1] = struct.unpack('<d', data[(ind + SIZE_OF_DOUBLE):(ind + 2*SIZE_OF_DOUBLE)])[0]
                
                 ind += 2 * SIZE_OF_DOUBLE
            
             return coord
            
        if dimensionality == 3:    
             coord = np.empty(shape=(N,3), dtype=np.float32) # coordinates
             coord.fill(np.nan)
        
             for i in range(N):
                 coord[i][0] = struct.unpack('<d', data[ind:(ind + SIZE_OF_DOUBLE)])[0]
                 coord[i][1] = struct.unpack('<d', data[(ind + SIZE_OF_DOUBLE):(ind + 2*SIZE_OF_DOUBLE)])[0]
                 coord[i][2] = struct.unpack('<d', data[(ind + 2*SIZE_OF_DOUBLE):(ind + 3*SIZE_OF_DOUBLE)])[0]
                
                
                 ind += 3 * SIZE_OF_DOUBLE
            
             return coord
            
        else:
            print "Dimensionality %s is not supported!",dimensionality
            return -1
    

def read_global_index_array(filename):
    """
    Read global array with HVC(RA) neuronal ids
    """
    with open(filename, 'rb') as file:
        data = file.read()
    file.close()        
    N_RA = struct.unpack('<i', data[:SIZE_OF_INT])[0]
    Id_RA_global = np.zeros(N_RA)
    
    for i in range(N_RA):
        Id_RA_global[i] = struct.unpack('<i', data[(SIZE_OF_INT * (i+1)):(SIZE_OF_INT * (i+2))])[0]
        
    return (N_RA, Id_RA_global) 

def read_weights(filename):
    """
    Read all synaptic connections of RA to RA neurons
    """
    with open(filename, 'rb') as file:
        data = file.read()
           
        N_RA = struct.unpack('<i', data[:SIZE_OF_INT])[0]
        trial_number = struct.unpack('<i', data[SIZE_OF_INT:2*SIZE_OF_INT])[0]
        weights = np.zeros((N_RA, N_RA), np.float32)
        
        pos = 2*SIZE_OF_INT     
        
        #print "N_RA = ", N_RA
        #print "trial_number = ", trial_number
        
        #print "num bytes = ", len(data[2*SIZE_OF_INT:]) 
        #print "num_datapoints = ",len(data[2*SIZE_OF_INT:]) / SIZE_OF_DOUBLE    
        
        for i in range(N_RA):
            for j in range(N_RA):
                weights[i][j] = struct.unpack('<d', data[pos:(pos + SIZE_OF_DOUBLE)])[0]
                pos += SIZE_OF_DOUBLE
        
        return (N_RA, trial_number, weights)  
    
def read_maturation_properties(filename):
    """
    Read maturation properties of HVC-RA neurons
    """
    with open(filename, 'rb') as file:
        data = file.read()
           
        N = struct.unpack('<i', data[:SIZE_OF_INT])[0]
        trial_number = struct.unpack('<i', data[SIZE_OF_INT:2*SIZE_OF_INT])[0]
        
        mature_indicators = np.zeros(N, np.int32)
        maturation_rate = np.zeros(N, np.int32)
        Erest = np.zeros(N, np.float32)
        GCa = np.zeros(N, np.float32)
        
        ind = 2*SIZE_OF_INT     
        
        for i in range(N):
            mature_indicators[i] = struct.unpack('<i', data[ind:(ind+SIZE_OF_INT)])[0]
            maturation_rate[i] = struct.unpack('<i', data[(ind+SIZE_OF_INT):(ind+2*SIZE_OF_INT)])[0]
            Erest[i] = struct.unpack('<d', data[(ind+2*SIZE_OF_INT):(ind+2*SIZE_OF_INT+SIZE_OF_DOUBLE)])[0]
            GCa[i] = struct.unpack('<d', data[(ind+2*SIZE_OF_INT+SIZE_OF_DOUBLE):(ind+2*SIZE_OF_INT+2*SIZE_OF_DOUBLE)])[0]
            
            ind += 2*SIZE_OF_INT+2*SIZE_OF_DOUBLE
        
        return (N, trial_number, mature_indicators, maturation_rate, Erest, GCa)    

def read_mature_indicators(filename):
    """
    Read maturation indicators of HVC-RA neurons
    """
    with open(filename, 'rb') as file:
        data = file.read()
           
        N = struct.unpack('<i', data[:SIZE_OF_INT])[0]
        trial_number = struct.unpack('<i', data[SIZE_OF_INT:2*SIZE_OF_INT])[0]
        
        mature_indicators = np.zeros(N, np.int32)
        
        ind = 2*SIZE_OF_INT     
        
        for i in range(N):
            mature_indicators[i] = struct.unpack('<i', data[ind:(ind + SIZE_OF_INT)])[0]
            ind += SIZE_OF_INT
        
        return (N, trial_number, mature_indicators) 
        
def read_remodeled_indicators(filename):
    """
    Read axon-remodeling indicators of HVC-RA neurons
    """
    with open(filename, 'rb') as file:
        data = file.read()
           
        N = struct.unpack('<i', data[:SIZE_OF_INT])[0]
        trial_number = struct.unpack('<i', data[SIZE_OF_INT:2*SIZE_OF_INT])[0]
        
        remodeled_indicators = np.zeros(N, np.int32)
        
        ind = 2*SIZE_OF_INT     
        
        for i in range(N):
            remodeled_indicators[i] = struct.unpack('<i', data[ind:(ind + SIZE_OF_INT)])[0]
            ind += SIZE_OF_INT
        
        return (N, trial_number, remodeled_indicators) 

def read_replacement_history(filename):
    """
    Read replacement history of HVC-RA neurons
    """
    with open(filename, 'rb') as file:
        data = file.read()
           
        N = struct.unpack('<i', data[:SIZE_OF_INT])[0]
        trial_number = struct.unpack('<i', data[SIZE_OF_INT:2*SIZE_OF_INT])[0]
        
        replacement_history = np.zeros(N, np.int32)
        
        ind = 2*SIZE_OF_INT     
        
        for i in range(N):
            replacement_history[i] = struct.unpack('<i', data[ind:(ind + SIZE_OF_INT)])[0]
            ind += SIZE_OF_INT
        
        return (N, trial_number, replacement_history) 

def read_activity_history(filename):
    """
    Read activity history of HVC-RA neurons
    """
    with open(filename, 'rb') as file:
        data = file.read()
           
        N = struct.unpack('<i', data[:SIZE_OF_INT])[0]
        rate_window = struct.unpack('<i', data[SIZE_OF_INT:2*SIZE_OF_INT])[0]
        trial_number = struct.unpack('<i', data[2*SIZE_OF_INT:3*SIZE_OF_INT])[0]
        
        activity_history = np.zeros((N, rate_window), np.int32)
        
        ind = 3*SIZE_OF_INT     
        
        for i in range(N):
            for j in range(rate_window):
                activity_history[i][j] = struct.unpack('<i', data[ind:(ind + SIZE_OF_INT)])[0]
                ind += SIZE_OF_INT
        
        return (N, trial_number, activity_history)  

def read_axonal_delays(filename):
    """
    Read axonal delays between neurons
    """
    with open(filename, 'rb') as file:
        data = file.read()
           
        N = struct.unpack('<i', data[:SIZE_OF_INT])[0]
        trial_number = struct.unpack('<i', data[SIZE_OF_INT:2*SIZE_OF_INT])[0]
        axonal_delays = [[] for i in range(N)]
        
        pos = 2*SIZE_OF_INT     
        
        #print "N_RA = ", N_RA
        #print "trial_number = ", trial_number
        
        #print "num bytes = ", len(data[2*SIZE_OF_INT:]) 
        #print "num_datapoints = ",len(data[2*SIZE_OF_INT:]) / SIZE_OF_DOUBLE    
        
        for i in range(N):
            num_targets = struct.unpack('<i', data[pos:(pos+SIZE_OF_INT)])[0]
            pos += SIZE_OF_INT
            
            for j in range(num_targets):
                axonal_delays[i].append(struct.unpack('<d', data[pos:(pos + SIZE_OF_DOUBLE)])[0])
                pos += SIZE_OF_DOUBLE
        
        return (N, trial_number, axonal_delays)  

def read_jitter(filename):
    """
    Read results of test chain simulation with jitter
    """
    with open(filename, 'rb') as file:
        data = file.read()
           
        N = struct.unpack('<i', data[:SIZE_OF_INT])[0]
        num_test_trials = struct.unpack('<i', data[SIZE_OF_INT:2*SIZE_OF_INT])[0]
        
        probability_soma_spike = np.empty(N, np.float32)        
        average_num_soma_spikes_in_trial = np.empty(N, np.float32)        
        mean_first_soma_spike_time = np.empty(N, np.float32)        
        std_first_soma_spike_time = np.empty(N, np.float32)        
        
        probability_dend_spike = np.empty(N, np.float32)        
        average_num_dend_spikes_in_trial = np.empty(N, np.float32)        
        mean_first_dend_spike_time = np.empty(N, np.float32)        
        std_first_dend_spike_time = np.empty(N, np.float32)        
                
        
        pos = 2*SIZE_OF_INT     
        
        #print "N_RA = ", N_RA
        #print "trial_number = ", trial_number
        
        #print "num bytes = ", len(data[2*SIZE_OF_INT:]) 
        #print "num_datapoints = ",len(data[2*SIZE_OF_INT:]) / SIZE_OF_DOUBLE    
        
        for i in range(N):
            probability_soma_spike[i] = struct.unpack('<d', data[pos:(pos+SIZE_OF_DOUBLE)])[0]
            average_num_soma_spikes_in_trial[i] = struct.unpack('<d', data[(pos+SIZE_OF_DOUBLE):(pos+2*SIZE_OF_DOUBLE)])[0]
            mean_first_soma_spike_time[i] = struct.unpack('<d', data[(pos+2*SIZE_OF_DOUBLE):(pos+3*SIZE_OF_DOUBLE)])[0]
            std_first_soma_spike_time[i] = struct.unpack('<d', data[(pos+3*SIZE_OF_DOUBLE):(pos+4*SIZE_OF_DOUBLE)])[0]
            
            probability_dend_spike[i] = struct.unpack('<d', data[(pos+4*SIZE_OF_DOUBLE):(pos+5*SIZE_OF_DOUBLE)])[0]
            average_num_dend_spikes_in_trial[i] = struct.unpack('<d', data[(pos+5*SIZE_OF_DOUBLE):(pos+6*SIZE_OF_DOUBLE)])[0]
            mean_first_dend_spike_time[i] = struct.unpack('<d', data[(pos+6*SIZE_OF_DOUBLE):(pos+7*SIZE_OF_DOUBLE)])[0]
            std_first_dend_spike_time[i] = struct.unpack('<d', data[(pos+7*SIZE_OF_DOUBLE):(pos+8*SIZE_OF_DOUBLE)])[0]
        
            pos += 8*SIZE_OF_DOUBLE
            
        return (N, num_test_trials, \
                    probability_soma_spike, average_num_soma_spikes_in_trial, mean_first_soma_spike_time, std_first_soma_spike_time, \
                    probability_dend_spike, average_num_dend_spikes_in_trial, mean_first_dend_spike_time, std_first_dend_spike_time)  

def read_hhi(filename):
    """
    Functions reads binary file with HHI_final neuron data
    """
    with open(filename, mode = "rb") as file:
        fileContent = file.read()

    file.close()	#	close file

    dataPointsNumber = struct.unpack("<i", fileContent[:SIZE_OF_INT])[0] #	unpacking number of data points
    Nspikes = struct.unpack("<i", fileContent[SIZE_OF_INT:2*SIZE_OF_INT])[0]
    timeStep = struct.unpack("<d", fileContent[2*SIZE_OF_INT:(2*SIZE_OF_INT + SIZE_OF_DOUBLE)])[0]   
    
    
    ind = 2*SIZE_OF_INT + SIZE_OF_DOUBLE    
    
    t = np.zeros(dataPointsNumber)
    v = np.zeros(dataPointsNumber)
    I = np.zeros(dataPointsNumber)
    n = np.zeros(dataPointsNumber)
    m = np.zeros(dataPointsNumber)
    h = np.zeros(dataPointsNumber)
    w = np.zeros(dataPointsNumber)
    Ge = np.zeros(dataPointsNumber)
    Gi = np.zeros(dataPointsNumber)
    flag = np.zeros(dataPointsNumber)
    
    
    for i in range(dataPointsNumber):
        t[i] = struct.unpack("<d", fileContent[ind:(ind+SIZE_OF_DOUBLE)])[0]
        v[i] = struct.unpack("<d", fileContent[(ind+SIZE_OF_DOUBLE):(ind+2*SIZE_OF_DOUBLE)])[0]
        I[i] = struct.unpack("<d", fileContent[(ind+2*SIZE_OF_DOUBLE):(ind+3*SIZE_OF_DOUBLE)])[0]
        n[i] = struct.unpack("<d", fileContent[(ind+3*SIZE_OF_DOUBLE):(ind+4*SIZE_OF_DOUBLE)])[0]
        m[i] = struct.unpack("<d", fileContent[(ind+4*SIZE_OF_DOUBLE):(ind+5*SIZE_OF_DOUBLE)])[0]
        h[i] = struct.unpack("<d", fileContent[(ind+5*SIZE_OF_DOUBLE):(ind+6*SIZE_OF_DOUBLE)])[0]
        w[i] = struct.unpack("<d", fileContent[(ind+6*SIZE_OF_DOUBLE):(ind+7*SIZE_OF_DOUBLE)])[0]
        Ge[i] = struct.unpack("<d", fileContent[(ind+7*SIZE_OF_DOUBLE):(ind+8*SIZE_OF_DOUBLE)])[0]
        Gi[i] = struct.unpack("<d", fileContent[(ind+8*SIZE_OF_DOUBLE):(ind+9*SIZE_OF_DOUBLE)])[0]
        flag[i] = struct.unpack("<i", fileContent[(ind+9*SIZE_OF_DOUBLE):(ind+9*SIZE_OF_DOUBLE+SIZE_OF_INT)])[0]
        
        ind += 9*SIZE_OF_DOUBLE+SIZE_OF_INT
    
    return (t, v, I, n, m, h, w, Ge, Gi, flag, Nspikes)

def read_hh2(fileName):
    """
    Functions reads binary file with HH2_final neuron data
    """
   
    with open(fileName, mode = "rb") as file:
	fileContent = file.read()

    file.close()	#	close file

    dataPointsNumber = struct.unpack("<i", fileContent[:SIZE_OF_INT])[0] #	unpacking number of data points
    Nsoma = struct.unpack("<i", fileContent[SIZE_OF_INT:2*SIZE_OF_INT])[0]
    #print dataPointsNumber[0]*11
    Ndend = struct.unpack("<i", fileContent[2*SIZE_OF_INT:3*SIZE_OF_INT])[0]
    #timeStep = struct.unpack("<d", fileContent[3*SIZE_OF_INT:(3*SIZE_OF_INT+SIZE_OF_DOUBLE)])[0]
    
    t = np.zeros(dataPointsNumber)
    Vs = np.zeros(dataPointsNumber)
    Is = np.zeros(dataPointsNumber)
    n = np.zeros(dataPointsNumber)
    h = np.zeros(dataPointsNumber)
    Vd = np.zeros(dataPointsNumber)
    Id = np.zeros(dataPointsNumber)
    r = np.zeros(dataPointsNumber)
    c = np.zeros(dataPointsNumber)
    Ca = np.zeros(dataPointsNumber)
    Gexc_d = np.zeros(dataPointsNumber)
    Ginh_d = np.zeros(dataPointsNumber)
    Gexc_s = np.zeros(dataPointsNumber)
    Ginh_s = np.zeros(dataPointsNumber)
    
    flag = np.zeros(dataPointsNumber)
    Ei = np.zeros(dataPointsNumber)

    ind = 3*SIZE_OF_INT + SIZE_OF_DOUBLE
    
    for i in range(dataPointsNumber):
        t[i] = struct.unpack("<d", fileContent[ind:(ind+SIZE_OF_DOUBLE)])[0]        
        Vs[i] = struct.unpack("<d", fileContent[(ind+SIZE_OF_DOUBLE):(ind+2*SIZE_OF_DOUBLE)])[0]        
        Is[i] = struct.unpack("<d", fileContent[(ind+2*SIZE_OF_DOUBLE):(ind+3*SIZE_OF_DOUBLE)])[0]        
        n[i] = struct.unpack("<d", fileContent[(ind+3*SIZE_OF_DOUBLE):(ind+4*SIZE_OF_DOUBLE)])[0]        
        h[i] = struct.unpack("<d", fileContent[(ind+4*SIZE_OF_DOUBLE):(ind+5*SIZE_OF_DOUBLE)])[0]        
        Vd[i] = struct.unpack("<d", fileContent[(ind+5*SIZE_OF_DOUBLE):(ind+6*SIZE_OF_DOUBLE)])[0]        
        Id[i] = struct.unpack("<d", fileContent[(ind+6*SIZE_OF_DOUBLE):(ind+7*SIZE_OF_DOUBLE)])[0]        
        r[i] = struct.unpack("<d", fileContent[(ind+7*SIZE_OF_DOUBLE):(ind+8*SIZE_OF_DOUBLE)])[0]        
        c[i] = struct.unpack("<d", fileContent[(ind+8*SIZE_OF_DOUBLE):(ind+9*SIZE_OF_DOUBLE)])[0]        
        Ca[i] = struct.unpack("<d", fileContent[(ind+9*SIZE_OF_DOUBLE):(ind+10*SIZE_OF_DOUBLE)])[0]        
        Gexc_d[i] = struct.unpack("<d", fileContent[(ind+10*SIZE_OF_DOUBLE):(ind+11*SIZE_OF_DOUBLE)])[0]        
        Ginh_d[i] = struct.unpack("<d", fileContent[(ind+11*SIZE_OF_DOUBLE):(ind+12*SIZE_OF_DOUBLE)])[0]        
        Gexc_s[i] = struct.unpack("<d", fileContent[(ind+12*SIZE_OF_DOUBLE):(ind+13*SIZE_OF_DOUBLE)])[0]
        Ginh_s[i] = struct.unpack("<d", fileContent[(ind+13*SIZE_OF_DOUBLE):(ind+14*SIZE_OF_DOUBLE)])[0]
       
        Ei[i] = struct.unpack("<d", fileContent[(ind+14*SIZE_OF_DOUBLE):(ind+15*SIZE_OF_DOUBLE)])[0]
       
 #print Ei[i]
	
#	if i // 100 == 0:
#		print t[i]
        flag[i] = struct.unpack("<i", fileContent[(ind+15*SIZE_OF_DOUBLE):(ind+15*SIZE_OF_DOUBLE+SIZE_OF_INT)])[0]
        
        

        ind += 15*SIZE_OF_DOUBLE + SIZE_OF_INT
    
    return (t, Vs, Is, n, h, Vd, Id, r, c, Ca, Gexc_d, Ginh_d, Gexc_s, Ginh_s, Ei, flag, Nsoma, Ndend)


def read_hh2_buffer_full(filename):
    """
    Reads full output of HH2_buffer neuron
    """
    with open(filename, mode = "rb") as f:
        data = f.read()
        
        # calculate number of datapoints
        num_datapoints = len(data) / (10 * SIZE_OF_DOUBLE)
        
        t = np.empty(num_datapoints, np.float32)
        Vs = np.empty(num_datapoints, np.float32)
        Vd = np.empty(num_datapoints, np.float32)
        Gexc_d = np.empty(num_datapoints, np.float32)
        Ginh_d = np.empty(num_datapoints, np.float32)
        
        n = np.empty(num_datapoints, np.float32)
        h = np.empty(num_datapoints, np.float32)
        r = np.empty(num_datapoints, np.float32)
        c = np.empty(num_datapoints, np.float32)
        Ca = np.empty(num_datapoints, np.float32)
        
        #print num_datapoints
        
        start_ind = 0
        
        for i in range(num_datapoints):
            t[i] = struct.unpack("<d", data[start_ind:(start_ind+SIZE_OF_DOUBLE)])[0]
            Vs[i] = struct.unpack("<d", data[(start_ind+SIZE_OF_DOUBLE):(start_ind+2*SIZE_OF_DOUBLE)])[0]
            Vd[i] = struct.unpack("<d", data[(start_ind+2*SIZE_OF_DOUBLE):(start_ind+3*SIZE_OF_DOUBLE)])[0]
            Gexc_d[i] = struct.unpack("<d", data[(start_ind+3*SIZE_OF_DOUBLE):(start_ind+4*SIZE_OF_DOUBLE)])[0]
            Ginh_d[i] = struct.unpack("<d", data[(start_ind+4*SIZE_OF_DOUBLE):(start_ind+5*SIZE_OF_DOUBLE)])[0]
        
            n[i] = struct.unpack("<d", data[(start_ind+5*SIZE_OF_DOUBLE):(start_ind+6*SIZE_OF_DOUBLE)])[0]
            h[i] = struct.unpack("<d", data[(start_ind+6*SIZE_OF_DOUBLE):(start_ind+7*SIZE_OF_DOUBLE)])[0]
            r[i] = struct.unpack("<d", data[(start_ind+7*SIZE_OF_DOUBLE):(start_ind+8*SIZE_OF_DOUBLE)])[0]
            c[i] = struct.unpack("<d", data[(start_ind+8*SIZE_OF_DOUBLE):(start_ind+9*SIZE_OF_DOUBLE)])[0]
            Ca[i] = struct.unpack("<d", data[(start_ind+9*SIZE_OF_DOUBLE):(start_ind+10*SIZE_OF_DOUBLE)])[0]
        
            start_ind += 10*SIZE_OF_DOUBLE
            
        return t, Vs, Vd, Gexc_d, Ginh_d, n, h, r, c, Ca

def get_RA2RA_graph(filename):
    """
    Function returns edges and weights for RA to RA connections
    """
    (N_RA, targets_RA_RA_ID, targets_RA_RA_G) = read_connections(filename)
    
    edges = []
    weights = []

    for i in range(N_RA):    
        weights.extend(targets_RA_RA_G[i])
        temp = [(i,j) for j in targets_RA_RA_ID[i]]
        edges.extend(temp)
    
    #weights = map(lambda x: x*10, weights)
        
    return (N_RA, edges, weights)
    
def get_RA2I_graph(filename):
    """
    Function returns edges and weights for RA to I connections
    """
    (N_RA, targets_RA_I_ID, targets_RA_I_G) = read_connections(filename)
    
    edges = []
    weights = []

    for i in range(N_RA):    
        weights.extend(targets_RA_I_G[i])
        temp = [(i,j+N_RA) for j in targets_RA_I_ID[i]]
        edges.extend(temp)
    
    #weights = map(lambda x: x*10, weights)
        
    return (N_RA, edges, weights)

def get_I2RA_graph(N_RA, filename):
    """
    Function returns edges and weights for I to RA connections
    """
    (N_I, targets_I_RA_ID, targets_I_RA_G) = read_connections(filename)
    
    edges = []
    weights = []

    for i in range(N_I):    
        weights.extend(targets_I_RA_G[i])
        temp = [(i+N_RA,j) for j in targets_I_RA_ID[i]]
        edges.extend(temp)
    
    #weights = map(lambda x: x*10, weights)
        
    return (N_I, edges, weights)  


def read_execution_time(filename):
    """
    Reads average execution time of a trial from file
    """
    with open(filename, "rb") as file:
        data = file.read()
        file.close()
    
    N_RA = struct.unpack("<i", data[:SIZE_OF_INT])[0] # num of RA neurons
    N_I = struct.unpack("<i", data[SIZE_OF_INT:2*SIZE_OF_INT])[0] # num of I neurons
    np = struct.unpack("<i", data[2*SIZE_OF_INT:3*SIZE_OF_INT])[0] # num of processses
    
    timestep = struct.unpack("<d", data[3*SIZE_OF_INT:(3*SIZE_OF_INT+SIZE_OF_DOUBLE)])[0] # simulation timestep
    network_update_frequency = struct.unpack("<d", data[(3*SIZE_OF_INT+SIZE_OF_DOUBLE):(3*SIZE_OF_INT+2*SIZE_OF_DOUBLE)])[0]    
    average_execution_time = struct.unpack("<d", data[(3*SIZE_OF_INT+2*SIZE_OF_DOUBLE):(3*SIZE_OF_INT+3*SIZE_OF_DOUBLE)])[0]    
    
    return (N_RA, N_I, np, timestep, network_update_frequency, average_execution_time)
 
def read_time_info(filename):
    with open(filename, "rb") as file:
        data = file.read()
        file.close()
    
    trial_number = struct.unpack("<i", data[:SIZE_OF_INT])[0]
    simulation_time = struct.unpack("<d", data[SIZE_OF_INT:(SIZE_OF_INT+SIZE_OF_DOUBLE)])[0]
    N = struct.unpack("<i", data[(SIZE_OF_INT+SIZE_OF_DOUBLE):(2*SIZE_OF_INT+SIZE_OF_DOUBLE)])[0] 
    #print len(data[(2*SIZE_OF_INT+SIZE_OF_DOUBLE):])
    ind = 2*SIZE_OF_INT+SIZE_OF_DOUBLE    
    spike_times = []
    neuron_fired = []
    
    for i in range(N):
        spike_array_size = struct.unpack("<i", data[ind:(ind+SIZE_OF_INT)])[0]
        single_neuron_spikes = []        
        single_neuron_fired = []        
        ind += SIZE_OF_INT
        for j in range(spike_array_size):
            temp = struct.unpack("<d", data[ind:(ind + SIZE_OF_DOUBLE)])[0]
            single_neuron_spikes.append(temp)
            single_neuron_fired.append(i)
            ind += SIZE_OF_DOUBLE
        
        if len(single_neuron_spikes) > 0:
            spike_times.append(single_neuron_spikes)
            neuron_fired.append(single_neuron_fired)
    #spike_times = struct.unpack("<{0}d".format(N_RA), data[(2*SIZE_OF_INT+SIZE_OF_DOUBLE):])
    return (trial_number, simulation_time, spike_times, neuron_fired)
    
def read_maturation_time_sequence(filename):
    """
    Reads neuron states of target neurons such as remodeled and mature indicators; gaba reverse potential
     and firing rate as functions of time measured in trial_numbers
    """
    with open(filename, "rb") as file:
        data = file.read()
       
        num_neurons = struct.unpack("<i", data[:SIZE_OF_INT])[0] # number of neurons     
        
        ind = SIZE_OF_INT    
        
        t = [] # time in trial numbers    
        remodeled = [[] for i in xrange(num_neurons)] # indicator for neuron to be remodeled
        gaba_potential = [[] for i in xrange(num_neurons)] # reverse GABA potential    
        firing_rate = [[] for i in xrange(num_neurons)] # neuron firing rate
        
        num_datapoints = (len(data) - SIZE_OF_INT) / (SIZE_OF_INT + num_neurons * (1*SIZE_OF_INT + 2*SIZE_OF_DOUBLE)) # number of datapoints in file
        
        for i in xrange(num_datapoints):
            t.append(struct.unpack("<i", data[ind:(ind + SIZE_OF_INT)])[0])
    
            for j in xrange(num_neurons):        
                remodeled[j].append(struct.unpack("<i", data[(ind + SIZE_OF_INT ):(ind + 2*SIZE_OF_INT)])[0])
                
                gaba_potential[j].append(struct.unpack("<d", data[(ind + 2*SIZE_OF_INT):(ind + 2*SIZE_OF_INT + SIZE_OF_DOUBLE)])[0])
                firing_rate[j].append(struct.unpack("<d", data[(ind + 2*SIZE_OF_INT + SIZE_OF_DOUBLE):(ind + 2*SIZE_OF_INT + 2*SIZE_OF_DOUBLE)])[0])
                
                ind += SIZE_OF_INT + 2*SIZE_OF_DOUBLE
            
            ind += SIZE_OF_INT
        
        return (t, remodeled, gaba_potential, firing_rate)
        

def read_synaptic_weights_time_sequence(filename):
    """
    Reads time dynamics of synaptic weights from source to target neurons
    """
    with open(filename, "rb") as file:
        data = file.read()
       
        num_neurons = struct.unpack("<i", data[:SIZE_OF_INT])[0] # number of neurons
        
        ind = SIZE_OF_INT    
            
        num_datapoints = (len(data) - SIZE_OF_INT) \
                                    / (num_neurons * num_neurons * SIZE_OF_DOUBLE + SIZE_OF_INT)
    
        t = []    
        weights = np.empty((num_datapoints, num_neurons, num_neurons), np.float32)
        
            
        for k in xrange(num_datapoints):    
            t.append(struct.unpack("<i", data[ind:(ind + SIZE_OF_INT)])[0])
            
            for i in xrange(num_neurons):
                for j in xrange(num_neurons):        
                    weights[k][i][j] = struct.unpack("<d", data[(ind + SIZE_OF_INT):(ind + SIZE_OF_INT + SIZE_OF_DOUBLE)])[0]
                    
                    ind += SIZE_OF_DOUBLE
            ind += SIZE_OF_INT
            
        return t, weights

def read_weight_statistics(filename):
    """
    Read mean synaptc weights and standard deviation of synaptic weights
    as a function of time
    """
    with open(filename, "rb") as file:
        data = file.read()
        file.close()
    
    trial_number = []
    mean = []
    std = []
    
    num_datapoints = len(data) / (SIZE_OF_INT + 2*SIZE_OF_DOUBLE)
    
    ind = 0
    
    for i in xrange(num_datapoints):
        
        trial_number.append(struct.unpack("<i", data[ind:(ind+SIZE_OF_INT)])[0])
        mean.append(struct.unpack("<d", data[(ind+SIZE_OF_INT):(ind + SIZE_OF_INT + SIZE_OF_DOUBLE)])[0])
        std.append(struct.unpack("<d", data[(ind + SIZE_OF_INT + SIZE_OF_DOUBLE):(ind + SIZE_OF_INT + 2*SIZE_OF_DOUBLE)])[0])
        
        ind += SIZE_OF_INT + 2*SIZE_OF_DOUBLE
   
    return (trial_number, mean, std)

def read_replaced_neurons(filename):
    """
    Read replaced neurons from file
    """
    with open(filename, "rb") as file:
        data = file.read()
        
        replaced_neurons = {}
        
        
        ind = 0
        
        while ind < len(data):    
            trial_number = struct.unpack("<i", data[ind:(ind+SIZE_OF_INT)])[0]
            num_replaced = struct.unpack("<i", data[(ind+SIZE_OF_INT):(ind+2*SIZE_OF_INT)])[0]
            
            #print trial_number
            #print num_replaced
            
            replaced = []
            
            
            for i in range(num_replaced):
                replaced.append(struct.unpack("<i", data[(ind+(i+2)*SIZE_OF_INT):(ind+(i+3)*SIZE_OF_INT)])[0])
            
            #print replaced
            
            replaced_neurons[trial_number] = list(replaced)
            
            ind += 2*SIZE_OF_INT + num_replaced*SIZE_OF_INT
       
        return replaced_neurons

def read_num_neurons(filename):
    """
    Read number of neurons in HVC network file
    """
    with open(filename, "rb") as file:
        data = file.read()
        
        N_RA = struct.unpack("<i", data[0:SIZE_OF_INT])[0]
        N_I = struct.unpack("<i", data[SIZE_OF_INT:2*SIZE_OF_INT])[0]
        
   
        return (N_RA, N_I)

def read_num_synapses(filename):
    """
    Read number of active and supersynapses from file
    """
    with open(filename, "rb") as file:
        data = file.read()
        file.close()
    
    trial_number = []
    num_active = []
    num_super = []
    
    num_datapoints = len(data) / (3 * SIZE_OF_INT)
    
    ind = 0
    
    for i in xrange(num_datapoints):
        
        trial_number.append(struct.unpack("<i", data[ind:(ind+SIZE_OF_INT)])[0])
        num_active.append(struct.unpack("<i", data[(ind+SIZE_OF_INT):(ind + 2*SIZE_OF_INT)])[0])
        num_super.append(struct.unpack("<i", data[(ind + 2*SIZE_OF_INT):(ind + 3*SIZE_OF_INT)])[0])
        
        ind += 3*SIZE_OF_INT
   
    return (trial_number, num_active, num_super)

def read_synapses(filename):
    """
    Read active or super
    """
    with open(filename, "rb") as file:
        data = file.read()
        
        N = struct.unpack("<i", data[:SIZE_OF_INT])[0] # number of neurons
        trial_number = struct.unpack("<i", data[SIZE_OF_INT:2*SIZE_OF_INT])[0] # trial number
        
        targets = [[] for i in range(N)]
        
        ind = 2*SIZE_OF_INT        
        
        for i in range(N):
            num_targets = struct.unpack("<i", data[ind:(ind+SIZE_OF_INT)])[0] # number of targets for neuron
            ind += SIZE_OF_INT
            
            for j in range(num_targets):
                targets[i].append(struct.unpack("<i", data[ind:(ind+SIZE_OF_INT)])[0])
                ind += SIZE_OF_INT
                
        return (N, trial_number, targets)
   
def read_synaptic_info(filename):
    with open(filename, "rb") as file:
        data = file.read()
        file.close()
        
    num_datapoints = len(data) / (3*SIZE_OF_INT) # number of datapoints in file

    trial_num = [] # trial number
    num_active = [] # numver of active synapses
    num_super = [] # number os supersynapses

    ind = 0 
    
    for i in xrange(num_datapoints):
        trial_num.append(struct.unpack("<i", data[ind:(ind + SIZE_OF_INT)])[0])
        num_active.append(struct.unpack("<i", data[(ind + SIZE_OF_INT):(ind + 2*SIZE_OF_INT)])[0])
        num_super.append(struct.unpack("<i", data[(ind + 2*SIZE_OF_INT):(ind + 3*SIZE_OF_INT)])[0])
        ind += 3*SIZE_OF_INT
    
    return (trial_num, num_active, num_super)
    
def read_maturation_info(filename):
    with open(filename, "rb") as file:
        data = file.read()
        file.close()
    
    N_RA = struct.unpack("<i", data[:SIZE_OF_INT])[0]
    trial_number = struct.unpack("<i", data[SIZE_OF_INT:2*SIZE_OF_INT])[0]
    
    gaba_potential = []    
    firing_rate_short = []
    firing_rate_long = []
  
    remodeled = []    
    
    ind = 2*SIZE_OF_INT
    
    for i in xrange(N_RA):
        gaba_potential.append(struct.unpack("<d", data[ind:(ind + SIZE_OF_DOUBLE)])[0])
        firing_rate_short.append(struct.unpack("<d", data[(ind + SIZE_OF_DOUBLE):(ind + 2*SIZE_OF_DOUBLE)])[0])
        firing_rate_long.append(struct.unpack("<d", data[(ind + 2*SIZE_OF_DOUBLE):(ind + 3*SIZE_OF_DOUBLE)])[0])

        remodeled.append(struct.unpack("<i", data[(ind + 3*SIZE_OF_DOUBLE):(ind + 3*SIZE_OF_DOUBLE + SIZE_OF_INT)])[0])
       
        ind += 3*SIZE_OF_DOUBLE + SIZE_OF_INT
    
    return trial_number, gaba_potential, firing_rate_short, firing_rate_long, remodeled

def read_simTime_info(filename):
    with open(filename, "rb") as file:
        data = file.read()
        file.close()
    
    trial_number = struct.unpack("<i", data[:SIZE_OF_INT])[0]
    internal_time = struct.unpack("<d", data[SIZE_OF_INT:(SIZE_OF_INT + SIZE_OF_DOUBLE)])[0]
    network_time = struct.unpack("<d", data[(SIZE_OF_INT + SIZE_OF_DOUBLE):(SIZE_OF_INT + 2*SIZE_OF_DOUBLE)])[0]
        
    return trial_number, internal_time, network_time
    
def read_kick_response(filename):
    with open(filename, "rb") as file:
        data = file.read()
        file.close()
    
    N_RA = struct.unpack("<i", data[:SIZE_OF_INT])[0]  
    mu_soma = struct.unpack("<d", data[SIZE_OF_INT:(SIZE_OF_INT + SIZE_OF_DOUBLE)])[0] 
    sigma_soma = struct.unpack("<d", data[(SIZE_OF_INT + SIZE_OF_DOUBLE):(SIZE_OF_INT + 2*SIZE_OF_DOUBLE)])[0] 
    mu_dend = struct.unpack("<d", data[(SIZE_OF_INT + 2*SIZE_OF_DOUBLE):(SIZE_OF_INT + 3*SIZE_OF_DOUBLE)])[0] 
    sigma_dend = struct.unpack("<d", data[(SIZE_OF_INT + 3*SIZE_OF_DOUBLE):(SIZE_OF_INT + 4*SIZE_OF_DOUBLE)])[0] 
    Ei = struct.unpack("<d", data[(SIZE_OF_INT + 4*SIZE_OF_DOUBLE):(SIZE_OF_INT + 5*SIZE_OF_DOUBLE)])[0] 
    inh_kick = struct.unpack("<d", data[(SIZE_OF_INT + 5*SIZE_OF_DOUBLE):(SIZE_OF_INT + 6*SIZE_OF_DOUBLE)])[0] 
    NMDA_kick = struct.unpack("<d", data[(SIZE_OF_INT + 6*SIZE_OF_DOUBLE):(SIZE_OF_INT + 7*SIZE_OF_DOUBLE)])[0] 

    kick_number = struct.unpack("<i", data[(SIZE_OF_INT + 7*SIZE_OF_DOUBLE):(2*SIZE_OF_INT + 7*SIZE_OF_DOUBLE)])[0]
    num_soma_spikes = struct.unpack("<i", data[(2*SIZE_OF_INT + 7*SIZE_OF_DOUBLE):(3*SIZE_OF_INT + 7*SIZE_OF_DOUBLE)])[0]
    num_dend_spikes = struct.unpack("<i", data[(3*SIZE_OF_INT + 7*SIZE_OF_DOUBLE):(4*SIZE_OF_INT + 7*SIZE_OF_DOUBLE)])[0]
    average_dend_spike_time = struct.unpack("<d", data[(4*SIZE_OF_INT + 7*SIZE_OF_DOUBLE):(4*SIZE_OF_INT + 8*SIZE_OF_DOUBLE)])[0]
    std_dend_spike_time = struct.unpack("<d", data[(4*SIZE_OF_INT + 8*SIZE_OF_DOUBLE):(4*SIZE_OF_INT + 9*SIZE_OF_DOUBLE)])[0]
        
    return N_RA, mu_soma, sigma_soma, mu_dend, sigma_dend, Ei, inh_kick, NMDA_kick, kick_number, num_soma_spikes, \
            num_dend_spikes, average_dend_spike_time, std_dend_spike_time

def read_noise_check(filename):
    with open(filename, "rb") as file:
        data = file.read()
        file.close()
    
    N_RA = struct.unpack("<i", data[:SIZE_OF_INT])[0]  
    mu_soma = struct.unpack("<d", data[SIZE_OF_INT:(SIZE_OF_INT + SIZE_OF_DOUBLE)])[0] 
    sigma_soma = struct.unpack("<d", data[(SIZE_OF_INT + SIZE_OF_DOUBLE):(SIZE_OF_INT + 2*SIZE_OF_DOUBLE)])[0] 
    mu_dend = struct.unpack("<d", data[(SIZE_OF_INT + 2*SIZE_OF_DOUBLE):(SIZE_OF_INT + 3*SIZE_OF_DOUBLE)])[0] 
    sigma_dend = struct.unpack("<d", data[(SIZE_OF_INT + 3*SIZE_OF_DOUBLE):(SIZE_OF_INT + 4*SIZE_OF_DOUBLE)])[0] 
    
    num_soma_spikes = struct.unpack("<i", data[(SIZE_OF_INT + 4*SIZE_OF_DOUBLE):(2*SIZE_OF_INT + 4*SIZE_OF_DOUBLE)])[0]
    num_dend_spikes = struct.unpack("<i", data[(2*SIZE_OF_INT + 4*SIZE_OF_DOUBLE):(3*SIZE_OF_INT + 4*SIZE_OF_DOUBLE)])[0]
    
    total_sim_time = struct.unpack("<d", data[(3*SIZE_OF_INT + 4*SIZE_OF_DOUBLE):(3*SIZE_OF_INT + 5*SIZE_OF_DOUBLE)])[0]        
    mean_vs = struct.unpack("<d", data[(3*SIZE_OF_INT + 5*SIZE_OF_DOUBLE):(3*SIZE_OF_INT + 6*SIZE_OF_DOUBLE)])[0]
    std_vs = struct.unpack("<d", data[(3*SIZE_OF_INT + 6*SIZE_OF_DOUBLE):(3*SIZE_OF_INT + 7*SIZE_OF_DOUBLE)])[0]
    mean_vd = struct.unpack("<d", data[(3*SIZE_OF_INT + 7*SIZE_OF_DOUBLE):(3*SIZE_OF_INT + 8*SIZE_OF_DOUBLE)])[0]
    std_vd = struct.unpack("<d", data[(3*SIZE_OF_INT + 8*SIZE_OF_DOUBLE):(3*SIZE_OF_INT + 9*SIZE_OF_DOUBLE)])[0]
    
    
    return N_RA, mu_soma, sigma_soma, mu_dend, sigma_dend, num_soma_spikes, num_dend_spikes, total_sim_time,\
            mean_vs, std_vs, mean_vd, std_vd
    
    
def read_chain_test(filename):
    with open(filename, "rb") as file:
        data = file.read()
        file.close()
    
    N_RA = struct.unpack("<i", data[:SIZE_OF_INT])[0]  
    num_trials = struct.unpack("<i", data[(SIZE_OF_INT):(2*SIZE_OF_INT)])[0]
    
    firing_robustness = []    
    average_num_dend_spikes_in_trial = []
    average_num_soma_spikes_in_trial = []
    mean_burst_time = []
    std_burst_time = []    
    
    ind = 2*SIZE_OF_INT
    
    for i in xrange(N_RA):
        firing_robustness.append(struct.unpack("<d", data[ind:(ind + SIZE_OF_DOUBLE)])[0])
        
        average_num_dend_spikes_in_trial.append(struct.unpack("<d", data[(ind + SIZE_OF_DOUBLE):(ind + 2*SIZE_OF_DOUBLE)])[0])
        average_num_soma_spikes_in_trial.append(struct.unpack("<d", data[(ind + 2*SIZE_OF_DOUBLE):(ind + 3*SIZE_OF_DOUBLE)])[0])
        
        mean_burst_time.append(struct.unpack("<d", data[(ind + 3*SIZE_OF_DOUBLE):(ind + 4*SIZE_OF_DOUBLE)])[0])
        std_burst_time.append(struct.unpack("<d", data[(ind + 4*SIZE_OF_DOUBLE):(ind + 5 * SIZE_OF_DOUBLE)])[0])
    
        ind = ind + 5*SIZE_OF_DOUBLE
    return N_RA, num_trials, firing_robustness, average_num_dend_spikes_in_trial, average_num_soma_spikes_in_trial,mean_burst_time, std_burst_time

def read_last_dend_spike_times(filename):
    """
    Reads last dendritic spike times of each HVC(RA) neuron
    """
    with open(filename, "rb") as file:
        data = file.read()
        file.close()
    
    N_RA = struct.unpack("<i", data[:SIZE_OF_INT])[0]  
    trial_number = struct.unpack("<i", data[SIZE_OF_INT:2*SIZE_OF_INT])[0]  

    last_dend_spike_times = np.empty(N_RA, dtype=np.float32)
    
    for i in range(N_RA):
        last_dend_spike_times[i] = struct.unpack("<d", data[(2*SIZE_OF_INT + i*SIZE_OF_DOUBLE):(2*SIZE_OF_INT + (i+1)*SIZE_OF_DOUBLE)])[0]
            
    return last_dend_spike_times
    
    
def read_num_bursts_in_recent_trials(filename):
    """
    Read number of bursts produced by each HVC(RA) neuron during last 
    RATE_WINDOW_LONG trials
    """
    with open(filename, "rb") as file:
        data = file.read()
        file.close()
    
    N_RA = struct.unpack("<i", data[:SIZE_OF_INT])[0]  
    rate_window = struct.unpack("<i", data[SIZE_OF_INT:2*SIZE_OF_INT])[0]  
    trial_number = struct.unpack("<i", data[2*SIZE_OF_INT:3*SIZE_OF_INT])[0]  

    num_bursts_in_recent_trials = np.empty(shape=(N_RA,rate_window), dtype=np.int32)
    
    for i in range(N_RA):
        for j in range(rate_window):
            num_bursts_in_recent_trials[i][j] = struct.unpack("<i", data[(3+j+i*rate_window)*SIZE_OF_INT:(4+j+i*rate_window)*SIZE_OF_INT])[0]
            
    return num_bursts_in_recent_trials

def read_training_spread(filename):
    """
    Read spread of training neurons
    """
    with open(filename, "rb") as file:
        data = file.read()
      
        N_TR = struct.unpack("<i", data[:SIZE_OF_INT])[0]
        
        training_spread = np.empty(N_TR, np.float64)
        
        for i in range(N_TR):
            training_spread[i] = struct.unpack("<d", data[(SIZE_OF_INT+i*SIZE_OF_DOUBLE):(SIZE_OF_INT+(i+1)*SIZE_OF_DOUBLE)])[0]
    
        return training_spread
        
def read_training_neurons(filename):
    """
    Read ids of training neurons
    """
    with open(filename, "rb") as file:
        data = file.read()
        file.close()

    N_TR = struct.unpack("<i", data[:SIZE_OF_INT])[0]
    
    training_neurons = np.empty(N_TR, np.int32)
    
    for i in range(N_TR):
        training_neurons[i] = struct.unpack("<i", data[(i+1)*SIZE_OF_INT:(i+2)*SIZE_OF_INT])[0]

    return training_neurons
    
 
def read_fI_HVCRA(filename):
    """
    Function reads data from file for fI curve of HVC(RA) neuron
    """    
    with open(filename, mode = "rb") as file:
        data = file.read()

        num_points = struct.unpack("<i", data[:SIZE_OF_INT])[0]
        ampl_step = struct.unpack("<d", data[SIZE_OF_INT:(SIZE_OF_INT+SIZE_OF_DOUBLE)])[0]
        
        ampl = [float(i) * ampl_step for i in range(num_points)]
        ampl = np.array(ampl)
        
        num_spikes = np.empty(num_points, np.int32)
        num_bursts = np.empty(num_points, np.int32)
        
        start_ind = SIZE_OF_INT + SIZE_OF_DOUBLE
        
        for i in range(num_points):    
            num_spikes[i] = struct.unpack("<i", data[start_ind:(start_ind+SIZE_OF_INT)])[0]
            num_bursts[i] = struct.unpack("<i", data[(start_ind+SIZE_OF_INT):(start_ind+2*SIZE_OF_INT)])[0]
            
            start_ind += 2*SIZE_OF_INT
        
        return ampl, num_spikes, num_bursts

def read_conductance_response(filename):
    """
    Function reads time of response of HVC-RA neuron to excitatory conductance pulse
    """
    with open(filename, mode = "rb") as file:
        data = file.read()

        num_points = struct.unpack("<i", data[:SIZE_OF_INT])[0]
        G_step = struct.unpack("<d", data[SIZE_OF_INT:(SIZE_OF_INT+SIZE_OF_DOUBLE)])[0]
        
        G = [float(i) * G_step for i in range(num_points)]
        G = np.array(G)
        
        burst_onset_times = np.empty(num_points, np.float32)
        spike_times = np.empty(num_points, np.float32)
        
        start_ind = SIZE_OF_INT + SIZE_OF_DOUBLE
        
        for i in range(num_points):    
            burst_onset_times[i] = struct.unpack("<d", data[start_ind:(start_ind+SIZE_OF_DOUBLE)])[0]
            spike_times[i] = struct.unpack("<d", data[(start_ind+SIZE_OF_DOUBLE):(start_ind+2*SIZE_OF_DOUBLE)])[0]
            
            start_ind += 2*SIZE_OF_DOUBLE
        
        return G, burst_onset_times, spike_times

def read_inhAndExc_test(filename):
    """
    Read results of inhibitory and excitatory inputs test
    """
    with open(filename, "rb") as file:
        data = file.read()
        
        N = struct.unpack("<i", data[:SIZE_OF_INT])[0]
        num_trials = struct.unpack("<i", data[SIZE_OF_INT:2*SIZE_OF_INT])[0]
        
        probability_for_relevant_spike = np.empty(N, np.float32)
        probability_for_nonrelevant_spike = np.empty(N, np.float32)
        
        mean_relevant_spike_time = np.empty(N, np.float32)
        std_relevant_spike_time = np.empty(N, np.float32)
        
        mean_spike_time = np.empty(N, np.float32)
        std_spike_time = np.empty(N, np.float32)
        
        ind = 2*SIZE_OF_INT
        
        for i in range(N): 
            probability_for_relevant_spike[i] = struct.unpack("<d", data[ind:(ind+SIZE_OF_DOUBLE)])[0]
            mean_relevant_spike_time[i] = struct.unpack("<d", data[(ind+SIZE_OF_DOUBLE):(ind+2*SIZE_OF_DOUBLE)])[0]
            std_relevant_spike_time[i] = struct.unpack("<d", data[(ind+2*SIZE_OF_DOUBLE):(ind+3*SIZE_OF_DOUBLE)])[0]
            
            probability_for_nonrelevant_spike[i] = struct.unpack("<d", data[(ind+3*SIZE_OF_DOUBLE):(ind+4*SIZE_OF_DOUBLE)])[0]
            mean_spike_time[i] = struct.unpack("<d", data[(ind+4*SIZE_OF_DOUBLE):(ind+5*SIZE_OF_DOUBLE)])[0]
            std_spike_time[i] = struct.unpack("<d", data[(ind+5*SIZE_OF_DOUBLE):(ind+6*SIZE_OF_DOUBLE)])[0]
            
            
            ind += 6*SIZE_OF_DOUBLE
            
    return probability_for_relevant_spike, mean_relevant_spike_time, std_relevant_spike_time, \
            probability_for_nonrelevant_spike, mean_spike_time, std_spike_time
            
def read_clopath_test(filename):
    """
    Read results of clopath STDP model simulation
    """
    with open(filename, "rb") as file:
        data = file.read()
        
        num_datapoints = struct.unpack("<i", data[:SIZE_OF_INT])[0]
        
        print num_datapoints        
        
        time = np.empty(num_datapoints, np.float32)
        vd = np.empty(num_datapoints, np.float32)
        u_minus = np.empty(num_datapoints, np.float32)
        u_plus = np.empty(num_datapoints, np.float32)
        x = np.empty(num_datapoints, np.float32)
        w = np.empty(num_datapoints, np.float32)
        
        
        ind = SIZE_OF_INT
        
        for i in range(num_datapoints): 
            time[i] = struct.unpack("<d", data[ind:(ind+SIZE_OF_DOUBLE)])[0]
            vd[i] = struct.unpack("<d", data[(ind+SIZE_OF_DOUBLE):(ind+2*SIZE_OF_DOUBLE)])[0]
            u_minus[i] = struct.unpack("<d", data[(ind+2*SIZE_OF_DOUBLE):(ind+3*SIZE_OF_DOUBLE)])[0]
            u_plus[i] = struct.unpack("<d", data[(ind+3*SIZE_OF_DOUBLE):(ind+4*SIZE_OF_DOUBLE)])[0]
            x[i] = struct.unpack("<d", data[(ind+4*SIZE_OF_DOUBLE):(ind+5*SIZE_OF_DOUBLE)])[0]
            w[i] = struct.unpack("<d", data[(ind+5*SIZE_OF_DOUBLE):(ind+6*SIZE_OF_DOUBLE)])[0]
            
            ind += 6*SIZE_OF_DOUBLE
            
    return time, vd, u_minus, u_plus, x, w