# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 18:21:03 2018

@author: jingroup

Script analyzes statistics of all weights, weights of active synapses
and weights of supersynapses
"""
import matplotlib.pyplot as plt
import reading
import os
import numpy as np

trial_number = 4000

CONVERTION_CONSTANT = 10 # 1 weight in the model equals 10 pS
MAX_WEIGHT = 50

#dirname = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/maturationTransition1/"
dirname = "/home/eugene/results/immature/clusters/matTrans10/"

fileWeights = os.path.join(dirname, "weights_" + str(trial_number) + ".bin")
fileActiveSynapses = os.path.join(dirname, "RA_RA_active_connections_" + str(trial_number) + ".bin")
fileSuperSynapses = os.path.join(dirname, "RA_RA_super_connections_" + str(trial_number) + ".bin")
fileTraining = os.path.join(dirname, "training_neurons.bin")
fileMature = os.path.join(dirname, "mature_" + str(trial_number) + ".bin")
fileAxonalDelaysRA2RA =  os.path.join(dirname, "axonal_delays_RA2RA_" + str(trial_number) + ".bin")

#fileActive = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays3/RA_RA_active_connections_5300.bin"

#(_, _, active_synapses) = reading.read_synapses(fileActive)

(_, _, axonal_delays_RA2RA) = reading.read_axonal_delays(fileAxonalDelaysRA2RA)
(_, _, active_synapses) = reading.read_synapses(fileActiveSynapses)
(_, _, super_synapses) = reading.read_synapses(fileSuperSynapses)
(_, _, weights) = reading.read_weights(fileWeights) 
(_, _, mature_indicators) = reading.read_remodeled_indicators(fileMature)

training_neurons = reading.read_training_neurons(fileTraining)


print "Mature neurons: ",[i for i in np.where(mature_indicators == 1)[0] if i not in training_neurons]
print "Training neurons: ",training_neurons
for i in training_neurons:
    print "Training neuron {0} has {1} supersynapses : {2}".format(i, len(super_synapses[i]), super_synapses[i])

print axonal_delays_RA2RA[228][231]
print weights[228][231]

print axonal_delays_RA2RA[726][231]
print weights[726][231]


#print weights[156]

weights_active_synapses = []

for source_id, active_targets in enumerate(active_synapses):
    for active_target_id in active_targets:
        if (weights[source_id][active_target_id] > MAX_WEIGHT):
            print "Scaled synapse {0} -> {1} with weight {2}".format(source_id, active_target_id, CONVERTION_CONSTANT * weights[source_id][active_target_id])
        weights_active_synapses.append(CONVERTION_CONSTANT * weights[source_id][active_target_id])

weights_super_synapses = []

for source_id, super_targets in enumerate(super_synapses):
    for super_target_id in super_targets:
        weights_super_synapses.append(CONVERTION_CONSTANT * weights[source_id][super_target_id])     


nbins = 50

f = plt.figure()

ax1 = f.add_subplot(111)

hist, bin_edges = np.histogram(weights_active_synapses, bins=nbins)
#print bin_edges
#print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
ax1.step(bin_centers, hist, label="active", where="pre")

hist, bin_edges = np.histogram(weights_super_synapses, bins=nbins)
#print bin_edges
#print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
ax1.step(bin_centers, hist, label="super", where="pre")

hist, bin_edges = np.histogram(weights[:][:] * CONVERTION_CONSTANT, bins=nbins)
#print bin_edges
#print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
ax1.step(bin_centers, hist, label="all", where="pre")


ax1.set_xlabel('weights [pS]')
ax1.set_ylabel('# of synapses')

ax1.legend()

plt.show()
        

