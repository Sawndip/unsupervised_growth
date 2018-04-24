# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 18:58:35 2018

@author: jingroup

Script analyzes axonal delays between neurons
"""
import reading
import matplotlib.pyplot as plt
import numpy as np
import os

dirname = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/maturationTransition3/"

trial_number = 0

fileAxonalDelaysRA2I = os.path.join(dirname, "axonal_delays_RA2I_" + str(trial_number) + ".bin")
fileAxonalDelaysRA2RA = os.path.join(dirname, "axonal_delays_RA2RA_" + str(trial_number) + ".bin")
fileAxonalDelaysI2RA = os.path.join(dirname, "axonal_delays_I2RA_" + str(trial_number) + ".bin")

fileActive = os.path.join(dirname, "RA_RA_active_connections_" + str(trial_number) + ".bin")
fileSuper = os.path.join(dirname, "RA_RA_super_connections_" + str(trial_number) + ".bin")

(_, _, active_synapses) = reading.read_synapses(fileActive)
(_, _, super_synapses) = reading.read_synapses(fileSuper)


(_, _, axonal_delays_RA2I) = reading.read_axonal_delays(fileAxonalDelaysRA2I)
(_, _, axonal_delays_RA2RA) = reading.read_axonal_delays(fileAxonalDelaysRA2RA)
(_, _, axonal_delays_I2RA) = reading.read_axonal_delays(fileAxonalDelaysI2RA)

super_axonal_delays_RA2RA = []

for i, targets in enumerate(super_synapses):
    for target in targets:
        super_axonal_delays_RA2RA.append(axonal_delays_RA2RA[i][target])
        


#print axonal_delays_RA2I 

all_axonal_delays_RA2I = [delay for delays in axonal_delays_RA2I for delay in delays]
all_axonal_delays_RA2RA = [delay for delays in axonal_delays_RA2RA for delay in delays]
all_axonal_delays_I2RA = [delay for delays in axonal_delays_I2RA for delay in delays]

nbins = 50

f = plt.figure()

ax1 = f.add_subplot(111)

hist, bin_edges = np.histogram(all_axonal_delays_RA2I, bins=nbins)
print bin_edges
print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
ax1.step(bin_centers, hist / float(np.sum(hist)), label="HVC(RA) -> HVC(I)", where="pre")

hist, bin_edges = np.histogram(all_axonal_delays_RA2RA, bins=nbins)
print bin_edges
print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
ax1.step(bin_centers, hist / float(np.sum(hist)), label="HVC(RA) -> HVC(RA)", where="pre")

#hist, bin_edges = np.histogram(super_axonal_delays_RA2RA, bins=nbins)
#print bin_edges
#print hist
#bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
#ax1.step(bin_centers, hist / float(np.sum(hist)), label="super HVC(RA) -> HVC(RA)", where="pre")


hist, bin_edges = np.histogram(all_axonal_delays_I2RA, bins=nbins)
print bin_edges
print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
ax1.step(bin_centers, hist / float(np.sum(hist)), label="HVC(I) -> HVC(RA)", where="pre")

ax1.set_xlabel('Axonal time delay (ms)')
ax1.set_ylabel('norm # of connections')

ax1.legend()

plt.show()

