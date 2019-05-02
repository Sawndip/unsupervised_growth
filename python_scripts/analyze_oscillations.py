# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 15:23:48 2018

@author: jingroup

Script analyzes oscillations for many simulations
"""
import reading
import matplotlib.pyplot as plt
import numpy as np
import utils
import os

BIN_WIDTH = 1.0


def estimate_oscillations(filenameTrial, filenameMature):
    """
    Estimate oscillation prominence for mature neurons in trial
    """
    
    (_, _, mature_indicators) = reading.read_mature_indicators(filenameMature)
    
    
    mature_neurons = np.where(mature_indicators == 1)[0]
    
     # plot somatic spikes and burst density for a single trial
    mature_neurons_set = set(mature_neurons)    
    spike_times_s,  neuron_id_s,  ordered_spikes_s, neuron_ordered_id_s = utils.getSpikes(os.path.join(testDir, filenameTrial))
    
    
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
    
    ax = f.add_subplot(111)
    ax.plot(time, burst_density)    
    ax.set_ylabel('Burst density (1/ms)')
    ax.set_xlabel('Time (ms)')
    ax.set_xlim([0, 250])
    ax.set_ylim([0, 11])
    
    first_nonzero = np.where(burst_density > 0)[0][0]
    last_nonzero = np.where(burst_density > 0)[0][-1]
    
    mean_bd = np.mean(burst_density[first_nonzero:last_nonzero])
    std_bd = np.std(burst_density[first_nonzero:last_nonzero])
    
    print "Interval for non-zero burst density: ", time[first_nonzero], time[last_nonzero]
    print "Mean burst density: ", mean_bd
    print "Std burst density: ", std_bd
    print "Std / mean: ",std_bd / mean_bd

    return std_bd / mean_bd


simnames = ["matTrans68", "matTrans62", "matTrans69", "matTrans66", "matTrans64"]
trials = [37400, 58400, 27800, 71400, 57200]

velocities = [0.5, 1.0, 1.33, 2.0, 10.0]

oscillations = []

for simname, trial in zip(simnames, trials):
    testDir = "/mnt/hodgkin/eugene/results/immature/clusters/test/" + simname + "/trial" + str(trial)
    dataDir =  "/mnt/hodgkin/eugene/results/immature/clusters/" + simname
    
    filenameTrial = os.path.join(testDir, "test_spike_times_soma_5.bin")
    filenameMature = os.path.join(dataDir, "mature_"+str(trial)+".bin")

    oscillations.append(estimate_oscillations(filenameTrial, filenameMature))


plt.figure()
plt.scatter(velocities, oscillations)
plt.xlabel('Conduction velocity (units of normal)' )
plt.ylabel('Coefficient of variation of burst density')

def fsigmoid(x, a, b, c, d):
    return d + c / (1.0 + np.exp(-a*(x-b)))

from scipy.optimize import curve_fit

popt, pcov = curve_fit( fsigmoid, velocities, oscillations, bounds=([0.0, 0.1, 1.0, 0.25], [10.0,  2.0, 3.0, 0.75]) )

print popt

x = np.linspace(0, 15, 1000)
y = fsigmoid(x, *popt)

plt.plot(x, y, 'k')

plt.show()