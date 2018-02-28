# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 11:28:55 2017

@author: jingroup

Script analyzes deliveries from and to mature neurons relative to the first spike of postsynaptic cell
"""

import reading
import numpy as np
import os
import matplotlib.pyplot as plt
import numpy as np

trial_number = 5300

network_dir = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays3/"
#network_name = "out50_active30_1ms"
#sim_name = "ee0.015_ie0.03_newModel_burstOnsets_0_"

fileMature = os.path.join(network_dir, "mature_" + str(trial_number) +".bin")
fileActiveSynapses = os.path.join(network_dir, "RA_RA_active_connections_" + str(trial_number) +".bin")
fileAxonalDelaysRA2RA = os.path.join(network_dir, "axonal_delays_RA2RA_"+ str(trial_number) + ".bin")
fileSpikeTimes = "/home/eugene/Output/networks/chainGrowth/matureTest/test_spike_times_soma_2.bin"


#MARGIN_LATE = 0.5 # margin for the burst coming late
MARGIN_LATE = 0.0 # margin for the burst coming late

(N_RA, _, mature_indicators) = reading.read_mature_indicators(fileMature)
(_, _, active_synapses) = reading.read_synapses(fileActiveSynapses)
(_, _, axonal_delays_RA2RA) = reading.read_axonal_delays(fileAxonalDelaysRA2RA)

print "Number of HVC(RA) neurons = ",N_RA

(_, _, spike_times_soma, neuron_fired_soma) = reading.read_time_info(fileSpikeTimes)

#print "Dedritic spikes: ", neuron_fired_dend
#print "Dendritic spike times: ", spike_times_dend

#print "Somatic spikes: ", neuron_fired_soma
#print "Somatic spike times: ", spike_times_soma


#print ordered_soma_spikes_raw
#print ordered_soma_raw


first_spike_times = np.empty(N_RA, np.float32)
first_spike_times.fill(-1.0)

for n, time in zip(neuron_fired_soma, spike_times_soma):
    first_spike_times[n] = time[0]


delivered_times_and_id = []
arrivals = []

mature_neurons = np.where(mature_indicators == 1)[0]

for i in mature_neurons:
    if first_spike_times[i] > 0:
        for target in active_synapses[i]:
            if first_spike_times[target] > 0:
                time_difference = first_spike_times[i] + axonal_delays_RA2RA[i][target] - first_spike_times[target]
                #if time_difference > MARGIN_LATE:
                    #print "time difference = {0} source: {1} target: {2}".format(time_difference, i, target)
                arrivals.append(time_difference)


num_arrivals = len(arrivals)
num_late_arrivals = sum(arrival > MARGIN_LATE for arrival in arrivals)

print "Total number of deliveries: ",num_arrivals
print "Number of late deliveries: ",num_late_arrivals
print "Fraction of late deliveries: ",float(num_late_arrivals) / float(num_arrivals)


plt.figure()

plt.hist(arrivals, bins = 100)
#plt.title('Statistics of deliveries')
plt.xlabel('1st spike input time - target burst onset time (ms)')
plt.ylabel('# of inputs')
plt.xlim([-15,15])
plt.show()

#(_, simulation_time, burst_times_raw, neuron_fired) = reading.read_time_info(fileSpikePattern)
#burst_labels = reading.read_burst_labels(fileBurstLabels)

#print delivered_times_and_id

#print burst_times
#print neuron_fired
#N_RA = 6000



