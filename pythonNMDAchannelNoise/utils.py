# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 18:20:55 2018

@author: jingroup

Scipt contains useful functions
"""
import numpy as np
import bisect
import os
import reading
import matplotlib.pyplot as plt

def plotNormHist(a, nbins, ax, label):
    """
    Plots normalized histogram on axes ax
    """
    hist, bin_edges = np.histogram(a, bins=nbins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
    
    if np.isnan(hist).any():
        print "Histogram has NaN values!"

        if np.isnan(np.array(a)).any():
            print "Array contains NaN values!"
        
        if len(np.where(np.array(a)>0)[0]) > 0:
            print "Array contains negative values!"
        
        #print a
    elif np.sum(hist) == 0:
        print "Histogram has all zero values"
        
        if np.isnan(np.array(a)).any():
            print "Array contains NaN values!"
        
        if len(np.where(np.array(a)>0)[0]) > 0:
            print "Array contains negative values!"
                
        #print a
            
    ax.step(bin_centers, hist / float(np.sum(hist)), label=label, where="pre")

def plotNumSpikes(spike_times, ax, xlabel):
    """
    Plot histogram of spikes
    """
    num_spiked = len(spike_times)    
    
    num_spikes = np.empty(num_spiked, np.int8)
    
    for i, neuron_spikes in enumerate(spike_times):
        num_spikes[i] = len(neuron_spikes)
    
    numSpikesUnique, counts = np.unique(num_spikes, return_counts=True)

    ax.bar(numSpikesUnique, counts, align='center', width=0.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Count')

def getBurstsForNeuron(spike_times, burst_duration):
    """
    Find bursts for a neuron from spike times array. Burst duration are taken to be burst_duration
    
    Input: spike_times - array with spike times
              burst_duration - duration of the burst in ms
    """
    bursts = []
     
    currentBurst = []
        
    for spike in spike_times:
        if len(currentBurst) == 0 :
            currentBurst.append(spike)
        elif spike - currentBurst[0] <= burst_duration:
            currentBurst.append(spike)
        else:
            bursts.append(currentBurst)
            currentBurst = [spike]
    
    bursts.append(currentBurst)
    
    return bursts


def getBursts(spike_times, burst_duration):
    """
    Find bursts for each neuron from spike times array. Burst duration are taken to be burst_duration
    
    Input: spike_times - array with spike times returned by getSpikes function
              burst_duration - duration of the burst in ms
    """
    bursts = []
     
    for spikes in spike_times:
        currentBurst = []
        burstsForNeuron = []
        
        for spike in spikes:
            if len(currentBurst) == 0 :
                currentBurst.append(spike)
            elif spike - currentBurst[0] <= burst_duration:
                currentBurst.append(spike)
            else:
                burstsForNeuron.append(currentBurst)
                currentBurst = [spike]
        
        burstsForNeuron.append(currentBurst)
        
        bursts.append(burstsForNeuron)
        
    return bursts
        

def getSpikes(filename):
    """
    Obtain spikes from the raster data file
    
    Input: filename - absolute path to a file with data
    """
    (trial_number, simulation_time, spike_times, neuron_id) = reading.read_time_info(filename)
    
    
    num_spiked = len(spike_times)    
    
    print "Number of spiked neurons: ",num_spiked
    
    # order neurons based on their first spike time
    ordered_spikes, neuron_ordered_id = zip(*sorted(zip(spike_times, neuron_id)))        
        
    #id_map = {i[0] : j for i,j in zip(neuron_ordered_id, range(num_spiked))}
    
    return spike_times,  neuron_id,  ordered_spikes, neuron_ordered_id


def plotSpikes(spike_times, neuron_id, ax):
    """
    Plot spikes
    """
    for neuron_spikes, id in zip(spike_times, neuron_id):
        for spike in neuron_spikes:
            ax.vlines(spike, id[0] - 0.5, id[0] + 0.5)


def plot_out_weight_distribution(layer, weights, fig):
    """
    Plot output weight distribution for neurons in layer
    """
    fig.clf()
    
    out_weights = []
    
    for i in layer:
        out_weights.extend(weights[i][np.nonzero(weights[i,:])])
    
    ax = fig.add_subplot(111)
    ax.hist(out_weights, bins=20)


def calculate_burst_density_hist(burst_times, bin_width, tmin, tmax):
    """
    Calculate burst density - number of bursts in time bins
    
    Input: burst_times - array with burst times
              bin_width - desired time resolution
              tmin, tmax - range of time
              
    Output:
              time - time bin centers
              burst_density - number of bursts in bins
    """
    size = int((tmax-tmin) / bin_width) + 1
    time = np.array([tmin + float(i)*bin_width for i in range(size)])
    
    #print burst_times
    #print size
    
    num_bursts = np.zeros(size, np.int32)    
    
    for burst_time in burst_times:
        #print burst_time, int(burst_time / bin_width)
        num_bursts[int(round((burst_time - tmin) / bin_width))] += 1
        
    burst_density = num_bursts / bin_width
    
    return time, burst_density

def calculate_burst_density(burst_times, bin_width):
    """
    Calculate burst density - number of bursts in time bins
    
    Input: burst_times - array with burst times
              bin_width - desired time resolution
              
    Output:
              time - time bin centers
              burst_density - number of bursts in bins
    """
    size = int(burst_times[-1] / bin_width) + 1
    time = np.array([float(i)*bin_width + bin_width/2. for i in range(size)])
    
    #print burst_times
    #print size
    
    num_bursts = np.zeros(size, np.int32)    
    
    for burst_time in burst_times:
        #print burst_time, int(burst_time / bin_width)
        num_bursts[int(burst_time / bin_width)] += 1
        
    burst_density = num_bursts / bin_width
    
    return time, burst_density

def calculate_longAndLat(coord):
    """
    Computes longitude and latitude from arrays with 3d coordinates
    
    """
    num_neurons = coord.shape[0]
    
    latitude = np.empty(num_neurons, np.float32)
    longitude = np.empty(num_neurons, np.float32)

    for i in range(num_neurons):
        latitude[i] = np.arcsin(coord[i][2])
        longitude[i] = np.arctan2(coord[i][1], coord[i][0])
    

    return longitude, latitude
    
def get_parallel(latitude, num_points, tolerance, num_iter):
    """
    Computes parallel with constant latitude
    
    Input: latitude;
           num_points - number of points on parallel
           tolerance - value by which function estimated at final value is
           different from zero
           num_iter - max number of iterations to perform
    Output: x,y coordinates of parallel on Mollweide projection
    """
    x = np.empty(num_points, np.float32)
    y = np.empty(num_points, np.float32)
    
    longitude = np.linspace(-np.pi, np.pi, num_points) # array with longitudes   
    
    for i in range(num_points):
        x[i], y[i] = Mollweide_projection(longitude[i], latitude, tolerance, num_iter)
    
    return x,y
 


def get_meridian(longitude, num_points, tolerance, num_iter):
    """
    Computes meridian with constant longitude
    
    Input: longitude;
           num_points - number of points on meridian
           tolerance - value by which function estimated at final value is
           different from zero
           num_iter - max number of iterations to perform
    Output: x,y coordinates of meridian on Mollweide projection
    """
    x = np.empty(num_points, np.float32)
    y = np.empty(num_points, np.float32)
    
    latitude = np.linspace(-np.pi/2.0, np.pi/2.0, num_points) # array with latitudes   
    
    for i in range(num_points):
        x[i], y[i] = Mollweide_projection(longitude, latitude[i], tolerance, num_iter)
    
    return x,y
    
def Mollweide_projection(longitude, latitude, tolerance, num_iter):
    """
    Computes coordinates of Mollweide projection of the point on sphere 
    on 2d plane
    
    Input: longitude and latitude; 
           tolerance - value by which function estimated at final value is
           different from zero
           num_iter - max number of iterations to perform
    Output: x,y coordinates of point on Mollweide projection
    """
    if np.isclose(latitude, np.pi/2.0):
        theta_new = np.pi/2.0
    elif np.isclose(latitude, -np.pi/2.0):
        theta_new = -np.pi/2.0
    else:
        counter = 0 # counter of iterations
        theta_old = latitude
            
        while True:
            theta_new = theta_old - (2*theta_old + np.sin(2*theta_old) - np.pi*np.sin(latitude)) / (2 + 2*np.cos(2*theta_old))
            counter += 1        
            
            if np.abs(2*theta_new + np.sin(2*theta_new) - np.pi*np.sin(latitude)) < tolerance:
                break
            
            if counter > num_iter:
                print "Max number of iterations is reached!"
                break
            
            theta_old = theta_new
    #print np.cos(theta_new)
    x = 2 * np.sqrt(2) * longitude * np.cos(theta_new) / np.pi
    y = np.sqrt(2) * np.sin(theta_new)
    
    return x,y

def find_inhibitory_effect_of_hvcra(source_neurons, RA2I_ids, RA2I_weights, I2RA_ids, I2RA_weights):
    """
    Finds the inhibitory effect source neurons have on other HVC-RA pool neurons
    """
    # find interneurons that receive connections from source neurons
    target_interneurons, num_excitatory_inputs, input_excitatory_weight = \
                find_excited_interneurons(source_neurons, RA2I_ids, RA2I_weights)

    print "Excited interneurons = ",target_interneurons
    #print "Number of excited interneurons = ",len(target_interneurons)
    #print "Excitatory weights to interneurons: ",input_excitatory_weight
    
    nbins = 50

    f = plt.figure()
    
    ax1 = f.add_subplot(211)
    hist, bin_edges = np.histogram(num_excitatory_inputs, bins=nbins)
    #print bin_edges
    #print hist
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
    ax1.step(bin_centers, hist, where="pre")
    ax1.set_xlabel('# of excitatory inputs')
    ax1.set_ylabel('# of interneurons')

    ax2 = f.add_subplot(212)    
    hist, bin_edges = np.histogram(input_excitatory_weight, bins=nbins)
    #print bin_edges
    #print hist
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
    ax2.step(bin_centers, hist, where="pre")
    ax2.set_xlabel('Input excitatory weight')
    ax2.set_ylabel('# of interneurons')
        
    
    # find convergence of inhibition to all pool neurons except for source_neurons:
    N_RA = len(RA2I_ids)
    
    pool_neurons = [i for i in range(N_RA) if i not in source_neurons]
    
    pool_neurons_sorted, num_inhibitory_inputs, input_inhibitory_weight = \
            find_convergence_of_inhibitory_input(target_interneurons, pool_neurons, I2RA_ids, I2RA_weights)
            
    #print "Pool neurons that receive inhibition: ",pool_neurons_sorted
    print "Number of pool neurons that receive inhibition: ",len(pool_neurons_sorted)
    print "Number of inhibitory inputs received: ",num_inhibitory_inputs
    #print "Input inhibitory weights: ",input_inhibitory_weight
    
    f = plt.figure()
    
    ax1 = f.add_subplot(211)
    hist, bin_edges = np.histogram(num_inhibitory_inputs, bins=nbins)
    #print bin_edges
    #print hist
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
    ax1.step(bin_centers, hist, where="pre")
    ax1.set_xlabel('# of inhibitory inputs')
    ax1.set_ylabel('# of neurons')

    ax2 = f.add_subplot(212)    
    hist, bin_edges = np.histogram(input_inhibitory_weight, bins=nbins)
    #print bin_edges
    #print hist
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
    ax2.step(bin_centers, hist, where="pre")
    ax2.set_xlabel('Input inhibitory weight')
    ax2.set_ylabel('# of interneurons') 
    
    plt.show()

def find_convergence_of_inhibitory_input(interneurons, pool_neurons, I2RA_ids, I2RA_weights):
    """
    Counts convergence and total synaptic weight of inhibitory inputs from interneurons to
    pool neurons
    """
    pool_neurons_sorted = sorted(pool_neurons)
    num_inhibitory_inputs = np.zeros(len(pool_neurons), np.int32)
    input_inhibitory_weight = np.zeros(len(pool_neurons), np.float32)
    
    for i in interneurons:
        for target, weight in zip(I2RA_ids[i], I2RA_weights[i]):
            try:
                ind = index(pool_neurons_sorted, target)
                num_inhibitory_inputs[ind] += 1
                input_inhibitory_weight[ind] += weight
            except:
                continue
            
    return pool_neurons_sorted, num_inhibitory_inputs, input_inhibitory_weight
            
def find_excited_interneurons(source_neurons, RA2I_ids, RA2I_weights):
    """
    Finds all interneurons that receive connections from source neurons
    """
    target_neurons = []
    num_excitatory_inputs = []
    input_excitatory_weight = []
    
    for source in source_neurons:
        for target, weight in zip(RA2I_ids[source], RA2I_weights[source]):
            try:
                ind = index(target_neurons, target)
                
                input_excitatory_weight[ind] += weight
                num_excitatory_inputs[ind] += 1
            except:
                ind = bisect.bisect_left(target_neurons, target)
                target_neurons.insert(ind, target)
                input_excitatory_weight.insert(ind, weight)
                num_excitatory_inputs.insert(ind, 1)
                
    return target_neurons, num_excitatory_inputs, input_excitatory_weight

def get_num_active_and_super_synapses(dirname, end_trial, trialStep):
    """
    Calculate time sequence of total number of active and super connections in the network
    """
    num_timepoints = end_trial / trialStep + 1
    
    num_active = np.zeros(num_timepoints, np.int32)
    num_super = np.zeros(num_timepoints, np.int32)
    

    current_trial = 0
    timepoint = 0

    while current_trial <= end_trial:   
        fileActive = os.path.join(dirname, "RA_RA_active_connections_" + str(current_trial) + ".bin")
        fileSuper = os.path.join(dirname, "RA_RA_super_connections_" + str(current_trial) + ".bin")
        
        (_, _, active_synapses) = reading.read_synapses(fileActive)
        (_, _, super_synapses) = reading.read_synapses(fileSuper)        
        
        for i in range(len(active_synapses)):
            num_active[timepoint] += len(active_synapses[i])
            num_super[timepoint] += len(super_synapses[i])
                
        timepoint += 1
        current_trial += trialStep
        
    return num_active, num_super


def get_total_input_weight_time_sequence(dirname, end_trial, trialStep, target_neurons):
    """
    Outputs total excitatory conductance input to target_neurons at different
        points during simulation
    """
    num_timepoints = end_trial / trialStep + 1
    
    input_weights = np.zeros((len(target_neurons), num_timepoints), np.float32)

    current_trial = 0
    timepoint = 0

    while current_trial <= end_trial:   
        fileWeights = os.path.join(dirname, "weights_" + str(current_trial) + ".bin")
        fileActive = os.path.join(dirname, "RA_RA_active_connections_" + str(current_trial) + ".bin")
        
        (N, _, weights) = reading.read_weights(fileWeights)
        (_, _, active_synapses) = reading.read_synapses(fileActive)
        
        for source_id in range(N):
            for i, target_id in enumerate(target_neurons):
                if target_id in active_synapses[source_id]:
                    input_weights[i][timepoint] += weights[source_id][target_id]
                
        timepoint += 1
        current_trial += trialStep
        
    return input_weights

def get_source_input_weight_time_sequence(dirname, end_trial, trialStep, source_neurons, target_neurons):
    """
    Outputs total excitatory conductance input from source_neurons to target_neurons at different
        points during simulation
    """
    num_timepoints = end_trial / trialStep + 1
    
    input_weights = np.zeros((len(target_neurons), num_timepoints), np.float32)

    current_trial = 0
    timepoint = 0

    while current_trial <= end_trial:   
        fileWeights = os.path.join(dirname, "weights_" + str(current_trial) + ".bin")
        (_, _, weights) = reading.read_weights(fileWeights)
        
        for source_id in source_neurons:
            for i, target_id in enumerate(target_neurons):
                input_weights[i][timepoint] += weights[source_id][target_id]
                
        timepoint += 1
        current_trial += trialStep
        
    return input_weights

def get_maturation_time_sequence(dirname, end_trial, trialStep, neurons):
    """
    Outputs excitatory conductance input from source_neurons to target_neuron at different
        points during simulation
    """
    
    num_timepoints = end_trial / trialStep + 1
    
    maturation_indicators = np.zeros((len(neurons), num_timepoints), np.float32)

    current_trial = 0
    timepoint = 0

    while current_trial <= end_trial:   
        fileMature = os.path.join(dirname, "mature_" + str(current_trial) + ".bin")
        (_, _, mat_indicators) = reading.read_mature_indicators(fileMature)
        
        for i, neuron_id in enumerate(neurons):
            maturation_indicators[i][timepoint] = mat_indicators[neuron_id]
                
        timepoint += 1
        current_trial += trialStep
        
    return maturation_indicators         

def get_weight_time_sequence_from_source(dirname, end_trial, trialStep, source_neuron, target_neurons):
    """
    Outputs excitatory conductance input from source_neuron to target_neurons at different
        points during simulation
    """
    num_timepoints = end_trial / trialStep + 1
    
    synaptic_weights = np.zeros((len(target_neurons), num_timepoints), np.float32)

    current_trial = 0
    timepoint = 0

    while current_trial <= end_trial:   
        fileWeights = os.path.join(dirname, "weights_" + str(current_trial) + ".bin")
        (_, _, weights) = reading.read_weights(fileWeights)
        
        for i, target_id in enumerate(target_neurons):
            synaptic_weights[i][timepoint] = weights[source_neuron][target_id]
                
        timepoint += 1
        current_trial += trialStep
        
    return synaptic_weights               

    
def get_weight_time_sequence(dirname, end_trial, trialStep, source_neurons, target_neuron):
    """
    Outputs excitatory conductance input from source_neurons to target_neuron at different
        points during simulation
    """
    num_timepoints = end_trial / trialStep + 1
    
    synaptic_weights = np.zeros((len(source_neurons), num_timepoints), np.float32)

    current_trial = 0
    timepoint = 0

    while current_trial <= end_trial:   
        fileWeights = os.path.join(dirname, "weights_" + str(current_trial) + ".bin")
        (_, _, weights) = reading.read_weights(fileWeights)
        
        for i, source_id in enumerate(source_neurons):
            synaptic_weights[i][timepoint] = weights[source_id][target_neuron]
                
        timepoint += 1
        current_trial += trialStep
        
    return synaptic_weights               


def get_stdp_timeDif_sequence(dirname, end_trial, trialStep, N_RA, target_neurons, source_neuron):
    """
    Outputs time difference between first spike of source neuron and dendritic spike times of target neurons
        - quantity used in STDP rules
    """
    num_timepoints = end_trial / trialStep + 1
    
    timeDif_sequence = np.zeros((len(target_neurons), num_timepoints), np.float32)

    current_trial = 0
    timepoint = 0

    while current_trial <= end_trial:   
        fileSomaSpikeTimes = os.path.join(dirname, "spike_times_soma_" + str(current_trial) + ".bin")
        fileDendSpikeTimes = os.path.join(dirname, "spike_times_dend_" + str(current_trial) + ".bin")
        
        (_, _, soma_spike_times_raw, neuron_spiked) = reading.read_time_info(fileSomaSpikeTimes)
    
        first_soma_spike_times = np.empty(N_RA, np.float32)
        first_soma_spike_times.fill(-1.0)

        for n, time in zip(neuron_spiked, soma_spike_times_raw):
            first_soma_spike_times[n[0]] = time[0]
        
        (_, _, dend_spike_times_raw, neuron_spiked) = reading.read_time_info(fileDendSpikeTimes)
    
        first_dend_spike_times = np.empty(N_RA, np.float32)
        first_dend_spike_times.fill(-1.0)

        for n, time in zip(neuron_spiked, dend_spike_times_raw):
            first_dend_spike_times[n[0]] = time[0]
                
        
        for i, target in enumerate(target_neurons):
            timeDif_sequence[i][timepoint] = first_dend_spike_times[target] - first_soma_spike_times[source_neuron]
                
        timepoint += 1
        current_trial += trialStep
        
    return timeDif_sequence
    

def get_first_spike_time_sequence_relative_to_training(dirname, end_trial, trialStep, N, neurons, training_neurons):
    """
    Outputs first spike times of neurons relative to first spike time of training neurons
        at different points during simulation
    """
    num_timepoints = end_trial / trialStep + 1
    
    spike_times_sequence = np.zeros((len(neurons), num_timepoints), np.float32)

    current_trial = 0
    timepoint = 0

    while current_trial <= end_trial:   
        fileSpikeTimes = os.path.join(dirname, "spike_times_soma_" + str(current_trial) + ".bin")
        
        spike_times = neuron_first_spike_time_relative_to_training(N, neurons, training_neurons, fileSpikeTimes)
        
        for i, spike in enumerate(spike_times):
            spike_times_sequence[i][timepoint] = spike
                
        timepoint += 1
        current_trial += trialStep
        
    return spike_times_sequence
    

def neuron_first_spike_time_relative_to_training(N, neurons, training_neurons, filename):
    """
    Outputs first spike times of neurons relative to first spike time of training neurons
    """
    (_, _, spike_times_raw, neuron_spiked) = reading.read_time_info(filename)
    
    
    #print neuron_spiked
    #print spike_times_raw    
    
    first_spike_times = np.empty(N, np.float32)
    first_spike_times.fill(-1.0)

    for n, time in zip(neuron_spiked, spike_times_raw):
        first_spike_times[n[0]] = time[0]

    # calculate mean spike time of training neurons
    mean_training_first_spike = np.mean(first_spike_times[training_neurons])
    
    return first_spike_times[neurons] - mean_training_first_spike    


def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect.bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError

def get_event_times_neurons(filename, neurons):
    """
    Reads event times (spikes in somatic or dendritic compartments )of neurons in trial
    !!! Array neurons has to be sorted!!!
    """
    (_, _, event_times, neuron_fired) = reading.read_time_info(filename)
    
    event_times_neurons = [[] for i in range(len(neurons))]
    
    for neuron_event_times, neuron_fired_id in zip(event_times, neuron_fired):
        try: 
            ind = index(neurons, neuron_fired_id[0])
            event_times_neurons[ind] = neuron_event_times
            
        except ValueError:
            continue
            
    return event_times_neurons
