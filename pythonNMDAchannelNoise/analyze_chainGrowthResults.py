# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 12:54:34 2018

@author: jingroup

Script analyzes chain growth results
"""
import reading
import utils
#import matplotlib
import numpy as np
#matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import os

TRIAL_DURATION = 500.0
BURST_DURATION = 30.0
A_D = 10000.0
BIN_WIDTH = 1.0

#dataDir = "/mnt/hodgkin/eugene/Output/networks/chainGrowth/matTrans83"
dataDir = "/mnt/hodgkin/eugene/results/immature/clusters/matTrans63"
#trial = 6800
trial = 25200


N_RA, N_I = reading.read_num_neurons(os.path.join(dataDir, "num_neurons.bin"))

print "Number of HVC-RA neurons: ",N_RA
print "Number of HVC-I neurons: ",N_I

##############################################
# get somatic and dendritic spike times
##############################################
spike_times_s,  neuron_id_s,  ordered_spikes_s, neuron_ordered_id_s = utils.getSpikes(os.path.join(dataDir, "spike_times_soma_"+str(trial)+".bin"))
spike_times_d,  neuron_id_d,  ordered_spikes_d, neuron_ordered_id_d = utils.getSpikes(os.path.join(dataDir, "spike_times_dend_"+str(trial)+".bin"))

##########################################
# plot histograms of the number of spikes
##########################################
f = plt.figure()
ax1 = f.add_subplot(121)
utils.plotNumSpikes(spike_times_s, ax1, '# somatic spikes')

ax2 = f.add_subplot(122)
utils.plotNumSpikes(spike_times_d, ax2, '# dendritic spikes')

##############################################
# plot unordered somatic and dendritic spikes
##############################################
f = plt.figure()

ax1 = f.add_subplot(121)
utils.plotSpikes(spike_times_s, neuron_id_s, ax1)
ax1.set_ylabel('id')
ax1.set_xlabel('Time (ms)')
ax1.set_xlim([-5, TRIAL_DURATION])
ax1.set_ylim([-25, N_RA+25])
ax1.set_title('Somatic spikes')

ax2 = f.add_subplot(122)
utils.plotSpikes(spike_times_d, neuron_id_d, ax2)
ax2.set_ylabel('id')
ax2.set_xlabel('Time (ms)')
ax2.set_xlim([-5, TRIAL_DURATION])
ax2.set_ylim([-25, N_RA+25])
ax2.set_title('Dendritic spikes')

##############################################
# plot burst density
##############################################
training_spread = reading.read_training_spread(os.path.join(dataDir, "training_spread.bin"))
training_neurons = reading.read_training_neurons(os.path.join(dataDir, "training_neurons.bin"))

N_TR = len(training_neurons)

print "Number of training neurons: ",N_TR

bursts = utils.getBursts(spike_times_s, BURST_DURATION)
        
first_spikes_in_bursts = [[] for i in range(N_RA)]

for burstsNeuron,  id in zip(bursts, neuron_id_s):
    for burst in burstsNeuron:
        first_spikes_in_bursts[id[0]].append(burst[0])
        
        
all_first_spikes_in_burst = []

for first_spikes in first_spikes_in_bursts:
    if len(first_spikes) > 0:
        all_first_spikes_in_burst.extend(first_spikes)
        
min_first_spike = np.min(all_first_spikes_in_burst)

all_first_spikes_in_burst = [s - min_first_spike for s in all_first_spikes_in_burst]

#ordered_dend_spikes = [d - min_dend_spike_time for d in ordered_dend_spikes]

time, burst_density = utils.calculate_burst_density(sorted(all_first_spikes_in_burst), BIN_WIDTH)

#start_time = 0.0
#end_time = 120.0
        
#########################################################################
# calculate continuity for the burst pattern produced by mature neurons #
#########################################################################
(_, _, mature_indicators) = reading.read_mature_indicators(os.path.join(dataDir, "mature_"+str(trial)+".bin"))

print "Number of mature HVC-RA: ",len(np.where(mature_indicators == 1)[0])

min_first_spike_mature = 1e6
max_first_spike_mature = -1e6

for id in np.where(mature_indicators == 1)[0]:
    for first_spike in first_spikes_in_bursts[id]:
        if first_spike < min_first_spike_mature:
            min_first_spike_mature = first_spike
        if first_spike > max_first_spike_mature:
            max_first_spike_mature = first_spike

min_first_spike_mature -= min_first_spike
max_first_spike_mature -= min_first_spike

   
print "(Min first spike mature, max first spike mature) = ",min_first_spike_mature, max_first_spike_mature
print "Mean burst density: ",np.mean(burst_density[(time > min_first_spike_mature) & (time < max_first_spike_mature)])
print "Std burst density: ",np.std(burst_density[(time > min_first_spike_mature) & (time < max_first_spike_mature)])
print "std / mean = ",np.std(burst_density[(time > min_first_spike_mature) & (time < max_first_spike_mature)]) / np.mean(burst_density[(time > min_first_spike_mature) & (time < max_first_spike_mature)])

plt.figure()
plt.plot(time, burst_density)
plt.xlabel('Time (ms)')
plt.ylabel('# of bursts / ms')

#########################################################################
# calculate first spike interval for neurons firing after starters      #
#########################################################################

# find earliest firing time of a starter neuron
min_first_spike_training = 1e6

for n in training_neurons:
    if len(first_spikes_in_bursts[n]) > 0 and first_spikes_in_bursts[n][0] < min_first_spike_training:
        min_first_spike_training = first_spikes_in_bursts[n][0]

fsi = []
first_spike_time_for_fsi = []

for burstsNeuron,  id in zip(bursts, neuron_id_s):
    if len(burstsNeuron) == 1:
        if burstsNeuron[0][0] >= min_first_spike_training and len(burstsNeuron[0]) > 1:
            first_spike_time_for_fsi.append(burstsNeuron[0][0])
            fsi.append(burstsNeuron[0][1] - burstsNeuron[0][0])
    else:
        print "Neuron {0} produced {1} bursts\n".format(id[0], len(burstsNeuron))

plt.figure()
plt.scatter(first_spike_time_for_fsi, fsi)
plt.xlabel('First spike time (ms)')
plt.ylabel('First spike time interval (ms)')       

###########################################################
# calculate delivery time of the first spike in the burst #
###########################################################
(_, _, active_synapses) = reading.read_synapses(os.path.join(dataDir, "RA_RA_active_connections_"+str(trial)+".bin"))
(_, _, super_synapses) = reading.read_synapses(os.path.join(dataDir, "RA_RA_super_connections_"+str(trial)+".bin"))

(_, _, axonal_delays_RA2I) = reading.read_axonal_delays(os.path.join(dataDir, "axonal_delays_RA2I_"+str(trial)+".bin"))
(_, _, axonal_delays_RA2RA) = reading.read_axonal_delays(os.path.join(dataDir, "axonal_delays_RA2RA_"+str(trial)+".bin"))
(_, _, axonal_delays_I2RA) = reading.read_axonal_delays(os.path.join(dataDir, "axonal_delays_I2RA_"+str(trial)+".bin"))


nbins = 50

f = plt.figure()

ax = f.add_subplot(111)

utils.plotNormHist([delay for delays in axonal_delays_RA2I for delay in delays], nbins, ax, "HVC(RA) -> HVC(I)")
utils.plotNormHist([delay for delays in axonal_delays_I2RA for delay in delays], nbins, ax, "HVC(I) -> HVC(RA)")
utils.plotNormHist([delay for delays in axonal_delays_RA2RA for delay in delays], nbins, ax, "random HVC(RA) -> HVC(RA)")
utils.plotNormHist([axonal_delays_RA2RA[i][target] for i, targets in enumerate(active_synapses) for target in targets], nbins, ax, "active HVC(RA) -> HVC(RA)")
utils.plotNormHist([axonal_delays_RA2RA[i][target] for i, targets in enumerate(super_synapses) for target in targets ], nbins, ax, "super HVC(RA) -> HVC(RA)")

ax.set_xlabel('Axonal time delay (ms)')
ax.set_ylabel('norm # of connections')
#ax.set_xlim([0, 0.85])
#ax.set_ylim([0, 0.035])

ax.legend()


###########################################################
# calculate delivery time of the first spike in the burst #
###########################################################

ncol = 3
nrow = 3

num_examples = ncol * nrow

np.random.seed(1991)

mature_spiked_not_starters = []


for n in np.where(mature_indicators == 1)[0]:
    if n not in training_neurons and len(first_spikes_in_bursts) > 0:
        mature_spiked_not_starters.append(n)

if len(mature_spiked_not_starters) >= ncol * nrow:
    example_neurons = list(np.random.choice(mature_spiked_not_starters, size=num_examples, replace=False))
else:
    example_neurons = []
    
example_arrivals = [[] for i in range(num_examples)]
example_presyn = [[] for i in range(num_examples)]

input_arrivals_all = []
input_arrivals_matureToMature = []
presynSpikeDiff = []

for i in range(N_RA):
    if len(first_spikes_in_bursts[i]) > 0:
        for j, target in enumerate(active_synapses[i]):
            if len(first_spikes_in_bursts[target]) > 0:
                for source_first_spike in first_spikes_in_bursts[i]:
                    for target_first_spike in first_spikes_in_bursts[target]:
                        time_difference_arrival = source_first_spike + axonal_delays_RA2RA[i][target] - target_first_spike
                        time_difference_presyn = source_first_spike - target_first_spike
                        
                        #if time_difference_arrival > 20:
                        #    print "time difference = {0} source: {1} target: {2} delay = {3}".format(time_difference, i, target, axonal_delays_RA2RA[i][target])
                        
                        input_arrivals_all.append(time_difference_arrival)
                        
                        
                        if mature_indicators[i] == 1 and mature_indicators[target] == 1:
                            input_arrivals_matureToMature.append(time_difference_arrival)
                            presynSpikeDiff.append(time_difference_presyn)

                            
                        if target in example_neurons:
                            ind = example_neurons.index(target)
                            example_arrivals[ind].append(time_difference_arrival)
                            example_presyn[ind].append(time_difference_presyn)



print "Fraction of late input arrivals: ",float(sum(input > 0.0 for input in input_arrivals_all)) / len(input_arrivals_all)

if len(input_arrivals_matureToMature) > 0:
    print "Fraction of late input arrivals between mature neurons: ",float(sum(input > 0.0 for input in input_arrivals_matureToMature)) / len(input_arrivals_matureToMature)

f = plt.figure()

ax1 = f.add_subplot(131)
ax1.hist(input_arrivals_all, bins = 200)
ax1.set_title('Between all HVC-RA neurons')
ax1.set_xlabel('first spike arrival time - target first spike time (ms)')
ax1.set_ylabel('# of deliveries')
ax1.set_xlim([-15,15])

ax2 = f.add_subplot(132)
ax2.hist(input_arrivals_matureToMature, bins = 50)
ax2.set_title('Between mature HVC-RA neurons')
ax2.set_xlabel('first spike arrival time - target first spike time (ms)')
ax2.set_ylabel('# of deliveries')
ax2.set_xlim([-15,15])

ax3 = f.add_subplot(133)
ax3.hist(presynSpikeDiff, bins = 50)
ax3.set_title('Between mature HVC-RA neurons')
ax3.set_xlabel('first source spike time - target first spike time (ms)')
ax3.set_ylabel('# of spikes')
ax3.set_xlim([-15,15])

# now plot input times for example neurons
if len(example_neurons) > 0:
    f, axarr = plt.subplots(nrows=nrow, ncols=ncol)
    plt.suptitle('Examples of input times distributions')
    for i in range(nrow):
        for j in range(ncol):
            axarr[i,j].hist(example_arrivals[i*ncol+j], bins = 20)
            axarr[i,j].set_xlim([-15,15])
            axarr[i,j].set_title('neuron {0}'.format(example_neurons[i*ncol+j]))

# now plot source spike times for example neurons
if len(example_neurons) > 0:
    f, axarr = plt.subplots(nrows=nrow, ncols=ncol)
    plt.suptitle('Examples of presyn times distributions')
    
    for i in range(nrow):
        for j in range(ncol):
            axarr[i,j].hist(example_presyn[i*ncol+j], bins = 20)
            axarr[i,j].set_xlim([-15,15])
            axarr[i,j].set_title('neuron {0}'.format(example_neurons[i*ncol+j]))

##############################################
# plot ordered somatic and dendritic spikes
##############################################
f = plt.figure()

ax1 = f.add_subplot(121)
utils.plotSpikes(ordered_spikes_s, [[i] for i in range(len(ordered_spikes_s))], ax1)
ax1.set_ylabel('id')
ax1.set_xlabel('Time (ms)')
ax1.set_title('Ordered somatic spikes')

ax2 = f.add_subplot(122)
utils.plotSpikes(ordered_spikes_d,[[i] for i in range(len(ordered_spikes_d))], ax2)
ax2.set_ylabel('id')
ax2.set_xlabel('Time (ms)')
ax2.set_title('Ordered dendritic spikes')

##############################################
# get interneuron spike times
##############################################
spike_times,  neuron_id,  _, _ = utils.getSpikes(os.path.join(dataDir, "spike_times_interneuron_"+str(trial)+".bin"))

##########################################
# plot histograms of the number of spikes
##########################################
f = plt.figure()
ax = f.add_subplot(111)
utils.plotNumSpikes(spike_times, ax, '# spikes')

##############################################
# plot unordered interneuron spikes
##############################################
f = plt.figure()

ax1 = f.add_subplot(111)
utils.plotSpikes(spike_times, neuron_id, ax1)
ax1.set_ylabel('id')
ax1.set_xlabel('Time (ms)')
ax1.set_xlim([-5, TRIAL_DURATION])
ax1.set_ylim([-25, N_I+25])
ax1.set_title('Interneuron spikes')

##############################################
# plot number of active and supersynapses
##############################################
trial_number, num_active_synapses, num_supersynapses = reading.read_num_synapses(os.path.join(dataDir, "num_synapses.bin"))
        
f = plt.figure()

ax1 = f.add_subplot(211)
ax1.plot(trial_number, num_active_synapses)
ax1.set_ylabel("# active synapses")

ax2 = f.add_subplot(212)
ax2.plot(trial_number, num_supersynapses)
ax2.set_xlabel("Time (# trials)")
ax2.set_ylabel("# supersynapses")

##############################################
# plot weight distributions
##############################################
(_, _, weights_RA2RA) = reading.read_weights(os.path.join(dataDir, "weights_"+str(trial)+".bin"))  
(_, targets_ID, weights_RA2I, syn_lengths, axonal_delays) = reading.read_connections(os.path.join(dataDir, "RA_I_connections_"+str(trial)+".bin"))  
(_, targets_ID, weights_I2RA, syn_lengths, axonal_delays) = reading.read_connections(os.path.join(dataDir, "I_RA_connections_"+str(trial)+".bin"))  

f = plt.figure()

ax1 = f.add_subplot(131)
ax1.hist([w/A_D  for weights in weights_RA2RA for w in weights], bins=50)
ax1.set_ylabel("Counts")
ax1.set_xlabel("Synaptic weight HVC-RA -> HVC-RA")
ax1.set_yscale('log')

ax2 = f.add_subplot(132)
ax2.hist([w  for weights in weights_RA2I for w in weights], bins=50)
ax2.set_ylabel("Counts")
ax2.set_xlabel("Synaptic weight HVC-RA -> HVC-I")
ax2.set_yscale('log')

ax3 = f.add_subplot(133)
ax3.hist([w  for weights in weights_I2RA for w in weights], bins=50)
ax3.set_ylabel("Counts")
ax3.set_xlabel("Synaptic weight HVC-I -> HVC-RA")
ax3.set_yscale('log')

plt.show()
