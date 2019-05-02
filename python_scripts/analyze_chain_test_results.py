# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 19:37:06 2018

@author: jingroup

Script plots network test results
"""
import reading
import matplotlib.pyplot as plt
import numpy as np
import utils
import os

TRAINING_KICK_TIME = 100.0

simname = "matTrans62"
trial = 72000
testTrial = 0
testDir = "/mnt/hodgkin/eugene/results/immature/clusters/test/" + simname + "/trial" + str(trial)
dataDir =  "/mnt/hodgkin/eugene/results/immature/clusters/" + simname

filenameJitter = os.path.join(testDir, "jitter.bin")
filenameTrial = os.path.join(testDir, "test_spike_times_soma_" + str(testTrial) + ".bin")
filenameTrialInterneurons = os.path.join(testDir, "test_spike_times_interneuron_" + str(testTrial) + ".bin")
filenameMature = os.path.join(dataDir, "mature_"+str(trial)+".bin")

N, num_test_trials, \
    probability_soma_spike, average_num_soma_spikes_in_trial, mean_first_soma_spike_time, std_first_soma_spike_time,\
    probability_dend_spike, average_num_dend_spikes_in_trial, mean_first_dend_spike_time, std_first_dend_spike_time = reading.read_jitter(filenameJitter) 

#filenameReplacement = "/mnt/hodgkin/eugene/results/immature/clusters/matTrans62/replacement_history_32400.bin"

#time_from_previous_replacement = reading.read_replacement_history(filenameReplacement)
(_, _, mature_indicators) = reading.read_mature_indicators(filenameMature)


plt.figure()
plt.hist(probability_soma_spike)
plt.xlabel('Probability to fire a soma spike')
plt.ylabel('Count')

print len(np.where(np.array(probability_soma_spike) > 0.5)[0])


#robust_neurons = set(np.where(np.array(probability_soma_spike) > 0.5)[0])
robust_neurons = np.where(np.array(probability_soma_spike) > 0.5)[0]

mature_neurons = np.where(mature_indicators == 1)[0]

print mean_first_soma_spike_time

mean_first_soma_spike_time_mature = np.array(mean_first_soma_spike_time)[mature_neurons]

smallest_first_spike_time_mature = np.min(mean_first_soma_spike_time_mature)
mean_first_soma_spike_time_mature -= smallest_first_spike_time_mature
std_first_soma_spike_time_mature = np.array(std_first_soma_spike_time)[mature_neurons]



if len(np.where(std_first_soma_spike_time_mature < 0)[0]) > 0:
    print "Negative jitter of neuron that is mature!"

else:
    
    for n in np.where(std_first_soma_spike_time_mature > 10*np.median(std_first_soma_spike_time_mature))[0]:
        print "Neuron {0} has jitter {1}, burst time {2}, firing robustness {3}, <# spikes> {4} <# dend spikes> {5}".format(mature_neurons[n], \
            std_first_soma_spike_time_mature[n], mean_first_soma_spike_time_mature[n], probability_soma_spike[mature_neurons[n]], \
            average_num_soma_spikes_in_trial[mature_neurons[n]], average_num_dend_spikes_in_trial[mature_neurons[n]]) 

    # plot jitter
    f = plt.figure()
    
    ax1 = f.add_subplot(121)
    ax1.scatter(mean_first_soma_spike_time_mature, std_first_soma_spike_time_mature)
    ax1.set_xlabel('Burst time (ms)')
    ax1.set_ylabel('Std on first spike time (ms)')
    
    ax2 = f.add_subplot(122)
    ax2.hist(std_first_soma_spike_time_mature)
    ax2.set_ylabel('Counts')
    ax2.set_xlabel('Std on first spike time (ms)')

    mean_first_soma_spike_time_mature_sorted = sorted(mean_first_soma_spike_time_mature)

    # plot average first spikes and burst density    
    f = plt.figure()   
    plt.suptitle('Average')

    ax1 = f.add_subplot(121)
    for i, t in enumerate(mean_first_soma_spike_time_mature_sorted):
        ax1.vlines(t, i-0.5, i+0.5)
    ax1.set_ylabel('id')
    ax1.set_xlabel('Time (ms)')
    ax1.set_title('Ordered 1st somatic spikes')
    ax1.set_xlim([0, 150])
    ax1.set_ylim([-10, 310])

    BIN_WIDTH = 1.0
    
    time, burst_density = utils.calculate_burst_density_hist(mean_first_soma_spike_time_mature_sorted, BIN_WIDTH, -10.0, 600.0)
    
    ax2 = f.add_subplot(122)
    ax2.plot(time, burst_density)    
    ax2.set_ylabel('Burst density (1/ms)')
    ax2.set_xlabel('Time (ms)')
    ax2.set_xlim([0, 350])
    
    # plot somatic spikes and burst density for a single trial
    mature_neurons_set = set(mature_neurons)    
    spike_times_s,  neuron_id_s,  ordered_spikes_s, neuron_ordered_id_s = utils.getSpikes(os.path.join(testDir, filenameTrial))



    # find earliest firing time of a starter neuron
    training_neurons = reading.read_training_neurons(os.path.join(dataDir, "training_neurons.bin"))

    BURST_DURATION = 30.0
    bursts = utils.getBursts(spike_times_s, BURST_DURATION)


    min_first_spike_training = 1e6
    
    for spike_times, neuron_ids in zip(spike_times_s, neuron_id_s):
        if neuron_ids[0] in training_neurons:
            if spike_times[0] < min_first_spike_training:
                min_first_spike_training = spike_times[0]
    
    fsi = []
    first_spike_time_for_fsi = []
    
    for burstsNeuron,  id in zip(bursts, neuron_id_s):
        if len(burstsNeuron) == 1:
            if burstsNeuron[0][0] >= min_first_spike_training and len(burstsNeuron[0]) > 1:
                first_spike_time_for_fsi.append(burstsNeuron[0][0] - min_first_spike_training)
                fsi.append(burstsNeuron[0][1] - burstsNeuron[0][0])
        else:
            print "Neuron {0} produced {1} bursts\n".format(id[0], len(burstsNeuron))
    
    plt.figure()
    plt.scatter(first_spike_time_for_fsi, fsi)
    plt.xlabel('First spike time (ms)')
    plt.ylabel('First spike time interval (ms)')       
    plt.ylim([0,8])
    plt.xlim([0,350])
    
    ordered_spikes_s_mature = []
    neuron_ordered_id_s_mature = []
    
    for spikes, neuron_id in zip(spike_times_s, neuron_id_s):
        if neuron_id[0] in mature_neurons_set:
            ordered_spikes_s_mature.append(spikes)
            neuron_ordered_id_s_mature.append(neuron_id)
 
    ordered_spikes_s_mature = sorted(ordered_spikes_s_mature)

    smallest_spike_time = ordered_spikes_s_mature[0][0]
    for spikes in ordered_spikes_s_mature:
        for i in range(len(spikes)):
            spikes[i] -= smallest_spike_time
    
    first_spikes_s = [s[0] for s in ordered_spikes_s_mature]
    
    time, burst_density = utils.calculate_burst_density_hist(first_spikes_s, BIN_WIDTH, -10.0, 600.0)

    f = plt.figure()
    plt.suptitle('Single trial')

    ax1 = f.add_subplot(121)
    utils.plotSpikes(ordered_spikes_s_mature, [[i] for i in range(len(ordered_spikes_s_mature))], ax1)
    ax1.set_ylabel('id')
    ax1.set_xlabel('Time (ms)')
    ax1.set_title('Ordered somatic spikes')
    ax1.set_xlim([0, 350])
    ax1.set_ylim([-10, 650])

    ax2 = f.add_subplot(122)
    ax2.plot(time, burst_density)    
    ax2.set_ylabel('Burst density (1/ms)')
    ax2.set_xlabel('Time (ms)')
    ax2.set_xlim([0, 350])
    ax2.set_ylim([0, 11])
    
    
    f = plt.figure()
    plt.title('Single trial interneurons')
    
    spike_times_in,  neuron_id_in,  _, _ = utils.getSpikes(os.path.join(testDir, filenameTrialInterneurons))

    all_spikes_in = [spike for spikes in spike_times_in for spike in spikes]
    
    time, burst_density_interneurons = utils.calculate_burst_density(sorted(all_spikes_in), BIN_WIDTH)
    plt.plot(time, burst_density_interneurons)
    plt.ylabel('Spike density (1/ms)')
    plt.xlabel('Time (ms)')
    plt.xlim([0, 350])
    
   

#==============================================================================
# print len(robust_neurons)
# 
# spike_times_s,  neuron_id_s,  ordered_spikes_s, neuron_ordered_id_s = utils.getSpikes(filenameTrial)
# 
# 
# #print spike_times_s
# #print neuron_id_s
# 
# ordered_spikes_s_robust = []
# neuron_ordered_id_s_robust = []
# 
# for spikes, neuron_id in zip(spike_times_s, neuron_id_s):
#     if neuron_id[0] in robust_neurons:
#         ordered_spikes_s_robust.append(spikes)
#         neuron_ordered_id_s_robust.append(neuron_id)
# 
# #print ordered_spikes_s_robust
# #print neuron_ordered_id_s_robust
#         
# ordered_spikes_s_robust = sorted(ordered_spikes_s_robust)
# 
# 
# #print ordered_spikes_s_robust
# #==============================================================================
# # smallest_spike_time = 1e6
# # for spikes in ordered_spikes_s_robust:
# #     for spike in spikes:
# #         if spike < smallest_spike_time:
# #             smallest_spike_time = spike
# # 
# #==============================================================================
# 
# #ordered_spikes_s_robust = [map(int, x) for x in values]
# 
# smallest_spike_time = ordered_spikes_s_robust[0][0]
# for spikes in ordered_spikes_s_robust:
#     for i in range(len(spikes)):
#         spikes[i] -= smallest_spike_time
#     
#     
# f = plt.figure()
# 
# ax = f.add_subplot(111)
# utils.plotSpikes(ordered_spikes_s_robust, [[i] for i in range(len(ordered_spikes_s_robust))], ax)
# ax.set_ylabel('id')
# ax.set_xlabel('Time (ms)')
# ax.set_title('Ordered somatic spikes')
# ax.set_xlim([0, 250])
# ax.set_ylim([-10, 650])
# 
#==============================================================================

# plot burst density
#

plt.show()