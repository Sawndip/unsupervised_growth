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
    file.close()
    #print len(data)
    
    N_RA = struct.unpack("<i", data[:SIZE_OF_INT])[0]
    
    #print N_RA     
    
    targets_ID = []
    targets_G = []
    #targets_RA_I_ID = []
    #targets_RA_I_G = []    
    ind = SIZE_OF_INT    
    
    for i in range(N_RA):
        neurons_ID = struct.unpack("<i", data[ind:(ind+SIZE_OF_INT)])[0]
        number_of_targets = struct.unpack("<i", data[(ind+SIZE_OF_INT):(ind+2*SIZE_OF_INT)])[0]
        temp_ID = []
        temp_G = []    
        
        for j in range(number_of_targets):
            ind_targets_ID = ind + 2 * SIZE_OF_INT + j * (SIZE_OF_INT + SIZE_OF_DOUBLE)
            ind_targets_G = ind_targets_ID + SIZE_OF_INT
            
            #print ind_targets_ID
            ID = struct.unpack("<i", data[ind_targets_ID:(ind_targets_ID + SIZE_OF_INT)])[0]
            #print ID
            G = struct.unpack("<d", data[ind_targets_G:(ind_targets_G + SIZE_OF_DOUBLE)])[0]
            #print G            
            temp_ID.append(ID)
            temp_G.append(G)
        
        targets_ID.append(temp_ID)
        targets_G.append(temp_G)
        ind += 2 * SIZE_OF_INT + number_of_targets * (SIZE_OF_INT + SIZE_OF_DOUBLE)

    return (N_RA, targets_ID, targets_G)
    
def read_coordinates(filename):
    """
    Read coordinates of neurons from binary file
    """
    with open(filename, 'rb') as file:
            data = file.read()
    file.close()
    
    N = struct.unpack('<i', data[:SIZE_OF_INT])[0] # get number of neurons
    
    xx = [] # x-coordinates
    yy = [] # y-coordinates

    ind = SIZE_OF_INT # position in file    
    
    for i in range(N):
        x = struct.unpack('<d', data[ind:(ind + SIZE_OF_DOUBLE)])[0]
        y = struct.unpack('<d', data[(ind + SIZE_OF_DOUBLE):(ind + 2*SIZE_OF_DOUBLE)])[0]
        
        xx.append(x)
        yy.append(y)
        ind += 2 * SIZE_OF_DOUBLE
    
    return (xx, yy)


def read_weights(filename):
    """
    Read all synaptic connections of RA to RA neurons
    """
    with open(filename, 'rb') as file:
        data = file.read()
    file.close()        
    N_RA = struct.unpack('<i', data[:SIZE_OF_INT])[0]
    weights = np.zeros((N_RA, N_RA))
    
    pos = SIZE_OF_INT     
    
    for i in range(N_RA):
        for j in range(N_RA):
            weights[i][j] = struct.unpack('<d', data[pos:(pos + SIZE_OF_DOUBLE)])[0]
            pos += SIZE_OF_DOUBLE
    
    return (N_RA, weights)    
    
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
    
    t = np.zeros((dataPointsNumber, 1))
    Vs = np.zeros((dataPointsNumber, 1))
    Is = np.zeros((dataPointsNumber, 1))
    n = np.zeros((dataPointsNumber, 1))
    h = np.zeros((dataPointsNumber, 1))
    Vd = np.zeros((dataPointsNumber, 1))
    Id = np.zeros((dataPointsNumber, 1))
    r = np.zeros((dataPointsNumber, 1))
    c = np.zeros((dataPointsNumber, 1))
    Ca = np.zeros((dataPointsNumber, 1))
    Gexc_d = np.zeros((dataPointsNumber, 1))
    Ginh_d = np.zeros((dataPointsNumber, 1))
    Gexc_s = np.zeros((dataPointsNumber, 1))
    Ginh_s = np.zeros((dataPointsNumber, 1))
    
    flag = np.zeros((dataPointsNumber, 1))
    Ei = np.zeros((dataPointsNumber, 1))

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
    
def read_time_info(filename):
    with open(filename, "rb") as file:
        data = file.read()
        file.close()
    
    trial_number = struct.unpack("<i", data[:SIZE_OF_INT])[0]
    simulation_time = struct.unpack("<d", data[SIZE_OF_INT:(SIZE_OF_INT+SIZE_OF_DOUBLE)])[0]
    N_RA = struct.unpack("<i", data[(SIZE_OF_INT+SIZE_OF_DOUBLE):(2*SIZE_OF_INT+SIZE_OF_DOUBLE)])[0] 
    #print len(data[(2*SIZE_OF_INT+SIZE_OF_DOUBLE):])
    ind = 2*SIZE_OF_INT+SIZE_OF_DOUBLE    
    spike_times = []
    neuron_fired = []
    for i in range(N_RA):
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
    
def read_sim_info(filename):
    with open(filename, "rb") as file:
        data = file.read()
        file.close()
    
    trial_duration = struct.unpack("<d", data[:SIZE_OF_DOUBLE])[0]
    synapses_trials_update = struct.unpack("<i", data[SIZE_OF_DOUBLE:(SIZE_OF_DOUBLE+SIZE_OF_INT)])[0]
    weights_trials_update = struct.unpack("<i", data[(SIZE_OF_DOUBLE+SIZE_OF_INT):(SIZE_OF_DOUBLE+2*SIZE_OF_INT)])[0] 
   
    return (trial_duration, synapses_trials_update, weights_trials_update)
    
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
    
    gaba_potential = []    
    maturation_triggered = []
    mature = []
        
    ind = SIZE_OF_INT
    
    for i in xrange(N_RA):
        gaba_potential.append(struct.unpack("<d", data[ind:(ind + SIZE_OF_DOUBLE)])[0])
        maturation_triggered.append(struct.unpack("<i", data[(ind + SIZE_OF_DOUBLE):(ind + SIZE_OF_DOUBLE + SIZE_OF_INT)])[0])
        mature.append(struct.unpack("<i", data[(ind + SIZE_OF_DOUBLE + SIZE_OF_INT):(ind + SIZE_OF_DOUBLE + 2*SIZE_OF_INT)])[0])
        
        ind += SIZE_OF_DOUBLE + 2*SIZE_OF_INT
    
    return gaba_potential, maturation_triggered, mature

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
    
    num_dend_spikes = []
    mean_burst_time = []
    std_burst_time = []    
    
    ind = 2*SIZE_OF_INT
    
    for i in xrange(N_RA):    
        num_dend_spikes.append(struct.unpack("<i", data[ind:(ind + SIZE_OF_INT)])[0])
        
        mean_burst_time.append(struct.unpack("<d", data[(ind + SIZE_OF_INT):(ind + SIZE_OF_INT + SIZE_OF_DOUBLE)])[0])
        std_burst_time.append(struct.unpack("<d", data[(ind + SIZE_OF_INT + SIZE_OF_DOUBLE):(ind + SIZE_OF_INT + 2 * SIZE_OF_DOUBLE)])[0])
    
        ind = ind + SIZE_OF_INT + 2*SIZE_OF_DOUBLE
    return N_RA, num_trials, num_dend_spikes, mean_burst_time, std_burst_time
    