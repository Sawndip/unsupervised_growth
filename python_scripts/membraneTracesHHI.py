# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 17:43:01 2017

@author: jingroup

Script plots membrane traces for same HVC(I) neuron
"""
import reading
import matplotlib.pyplot as plt

dataDir = "/home/eugene/hodgkinData/gabaMaturation310117/membraneTraces/I/"

CURRENT_INJECTION = 100 # time of current injection

def get_spike_times(filename):
    """
    Extract neuronal spike times
    """    
    (t, _, _, _, _, _, _, _, _, flag, _) = reading.read_hhi(filename)

    flag_previous = flag[0]    

    spike_times = []
    
    for i in range(1, len(flag)):
        if flag[i] - flag_previous < -0.5:
            spike_times.append(t[i] - CURRENT_INJECTION)
        flag_previous = flag[i]
        

    return spike_times


def get_spikes_of_neuron(num_trials, neuron_id, dataDir):
    """
    Extract spike times of neuron with id neuron_id from different trials
    from files in directory dataDir
    """

    spike_times = []

    for i in range(num_trials):
        filename = dataDir + "I" + str(neuron_id) + "_trial" + str(i+1) + ".bin"
        spike_times.append(get_spike_times(filename))

    return spike_times

def get_spikes_of_neurons(num_trials, neurons, dataDir):
    """
    Extract spike times of neurons with id in list neurons from different trials
    from files in directory dataDir
    """
    spike_times = []
     
    for n in neurons:
        spike_times.append(get_spikes_of_neuron(num_trials, n, dataDir))
         
    return spike_times

neurons = range(20)
num_trials = 20

spike_times = get_spikes_of_neurons(num_trials, neurons, dataDir)
print spike_times


f = plt.figure()

colors = ['b', 'r', 'g', 'm', 'c', 'k']

ax = f.add_subplot(111)
for i in range(len(neurons)):
    for j in range(num_trials):
        for k in range(len(spike_times[i][j])):
            ax.vlines(spike_times[i][j][k], i*(num_trials+1)+j-0.25, i*(num_trials+1)+j+0.25, colors=colors[i%len(colors)], linewidth=5.0)

#ax.scatter(spike_times, random_ID)
#ax2.set_yticks(random_ID_stronglyConnected)

plt.yticks([])
xmin, xmax = ax.get_xlim()
ax.set_xlim([0, xmax+5])
ax.set_ylabel("trial #")
ax.set_xlabel("spike time (ms)")
ax.set_title("Spike times of HVC(I) neurons")


plt.show()
