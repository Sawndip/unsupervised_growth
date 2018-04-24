# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 16:17:42 2018

@author: jingroup

Script analyzes spatial distribution of connections in the network
"""
import matplotlib.pyplot as plt
import reading
import numpy as np

file_RA2I = "/home/eugene/Output/networks/chainGrowth/network2000/RA_I_connections.bin"
file_I2RA = "/home/eugene/Output/networks/chainGrowth/network2000/I_RA_connections.bin"

#file_RA2I = "/home/eugene/Output/networks/chainGrowth/test/RA_I_connections.bin"
#file_I2RA = "/home/eugene/Output/networks/chainGrowth/test/I_RA_connections.bin"


file_RA2RA = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/noImmatureOut2/RA_RA_super_connections.bin"

(N_RA, targets_id_RA2I, _, syn_lengths_RA2I, _) = reading.read_connections(file_RA2I)
#(_, targets_id_RA2RA, _, syn_lengths_RA2RA, _) = reading.read_connections(file_RA2RA)
(N_I, targets_id_I2RA, _, syn_lengths_I2RA, _) = reading.read_connections(file_I2RA)


#############################################
### Calculate fraction of connected targets
#############################################
fractions = []

for targets in targets_id_RA2I:
    fractions.append(float(len(targets)) / float(N_I))

print "Fraction of connected HVC-I by single HVC-RA: ", np.mean(np.array(fractions))

fractions = []

for targets in targets_id_I2RA:
    fractions.append(float(len(targets)) / float(N_RA))

print "Fraction of connected HVC-RA by single HVC-I: ", np.mean(np.array(fractions))
    
    
all_syn_lengths_RA2I = np.array([l for lengths in syn_lengths_RA2I for l in lengths])
#all_syn_lengths_RA2RA = np.array([l for lengths in syn_lengths_RA2RA for l in lengths])
all_syn_lengths_I2RA = np.array([l for lengths in syn_lengths_I2RA for l in lengths])

print all_syn_lengths_I2RA

nbins = 40

f = plt.figure()

ax1 = f.add_subplot(111)

hist, bin_edges = np.histogram(all_syn_lengths_RA2I, bins=nbins)
print bin_edges
print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
ax1.step(bin_centers, hist / float(N_RA), label="HVC(RA) -> HVC(I)", where="pre")

#hist, bin_edges = np.histogram(all_syn_lengths_RA2RA, bins=nbins)
#print bin_edges
#print hist
#bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
#x1.step(bin_centers, hist / float(N_RA), label="HVC(RA) -> HVC(RA)", where="pre")


hist, bin_edges = np.histogram(all_syn_lengths_I2RA, bins=nbins)
print bin_edges
print hist
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
ax1.step(bin_centers, hist / float(N_I), label="HVC(I) -> HVC(RA)", where="pre")


ax1.set_xlabel('Synaptic length ($\mu m$)')
ax1.set_ylabel('# of connections')

ax1.legend()

plt.show()

