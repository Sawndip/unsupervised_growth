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

trial_number = 100

CONVERTION_CONSTANT = 10 # 1 weight in the model equals 10 pS

dirname = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays/"

fileWeights = os.path.join(dirname, "weights_" + str(trial_number) + ".bin")
fileActiveSynapses = os.path.join(dirname, "RA_RA_active_connections_" + str(trial_number) + ".bin")
fileSuperSynapses = os.path.join(dirname, "RA_RA_super_connections_" + str(trial_number) + ".bin")


(_, _, active_synapses) = reading.read_synapses(fileActiveSynapses)
(_, _, super_synapses) = reading.read_synapses(fileSuperSynapses)
(_, _, weights) = reading.read_weights(fileWeights) 

print weights[156]

weights_active_synapses = []

for source_id, active_targets in enumerate(active_synapses):
    for active_target_id in active_targets:
        weights_active_synapses.append(CONVERTION_CONSTANT * weights[source_id][active_target_id])

weights_super_synapses = []

for source_id, super_targets in enumerate(super_synapses):
    for super_target_id in super_targets:
        weights_super_synapses.append(CONVERTION_CONSTANT * weights[source_id][super_target_id])     


nbins = 50

f = plt.figure()

ax1 = f.add_subplot(111)

hist, bin_edges = np.histogram(weights_active_synapses, bins=nbins)
print bin_edges
print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
ax1.step(bin_centers, hist, label="active", where="pre")

hist, bin_edges = np.histogram(weights_super_synapses, bins=nbins)
print bin_edges
print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
ax1.step(bin_centers, hist, label="super", where="pre")

hist, bin_edges = np.histogram(weights[:][:], bins=nbins)
print bin_edges
print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
ax1.step(bin_centers, hist, label="all", where="pre")


ax1.set_xlabel('weights [pS]')
ax1.set_ylabel('# of synapses')

ax1.legend()

plt.show()
        

