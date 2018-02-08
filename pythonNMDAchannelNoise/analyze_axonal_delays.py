# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 18:58:35 2018

@author: jingroup

Script analyzes axonal delays between neurons
"""
import reading
import matplotlib.pyplot as plt
import numpy as np

fileAxonalDelaysRA2I = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays/axonal_delays_RA2I_0.bin"
fileAxonalDelaysRA2RA = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays/axonal_delays_RA2RA_0.bin"
fileAxonalDelaysI2RA = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays/axonal_delays_I2RA_0.bin"


(_, _, axonal_delays_RA2I) = reading.read_axonal_delays(fileAxonalDelaysRA2I)
(_, _, axonal_delays_RA2RA) = reading.read_axonal_delays(fileAxonalDelaysRA2RA)
(_, _, axonal_delays_I2RA) = reading.read_axonal_delays(fileAxonalDelaysI2RA)

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
ax1.step(bin_centers, hist / float(np.max(hist)), label="HVC(RA) -> HVC(I)", where="pre")

hist, bin_edges = np.histogram(all_axonal_delays_RA2RA, bins=nbins)
print bin_edges
print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
ax1.step(bin_centers, hist / float(np.max(hist)), label="HVC(RA) -> HVC(RA)", where="pre")

hist, bin_edges = np.histogram(all_axonal_delays_I2RA, bins=nbins)
print bin_edges
print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
ax1.step(bin_centers, hist / float(np.max(hist)), label="HVC(I) -> HVC(RA)", where="pre")

ax1.set_xlabel('Axonal time delay (ms)')
ax1.set_ylabel('norm # of connections')

ax1.legend()

plt.show()

