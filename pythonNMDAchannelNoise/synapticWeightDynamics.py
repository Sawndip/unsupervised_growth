# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 10:36:42 2017

@author: jingroup

Script plots synaptic weight dynamics of a synapse
"""
import reading
import matplotlib.pyplot as plt

filenameWeights = "/home/eugene/Output/networks/test230217/weightsTimeSequence.bin"
filenameStates = "/home/eugene/Output/networks/test230217/maturation_time_sequence.bin"

(target, t, remodeled, mature, gaba_potential, firing_rate) = reading.read_maturation_time_sequence(filenameStates)

(source, target, t, weights) = reading.read_synaptic_weights_time_sequence(filenameWeights)

# plot synamptic weight dynamics
f = plt.figure()

ax = f.add_subplot(111)

source_neuron = 0 # real id of source neuron
target_neuron = 32 # real id of target neuron

ind_source = source.index(source_neuron) # index of source neuron in the list
ind_target = target.index(target_neuron) # index of target neuron in the list


ax.set_title("Synaptic weight {0} -> {1} vs. time".format(source_neuron, target_neuron))
ax.plot(t, weights[ind_source][ind_target])
ax.set_xlabel("t (# trial)")
ax.set_ylabel("w")

# plot firing rate and gaba potential
neuron_id = 32
xmax = 3000 # max trial to show

f = plt.figure()

ax1 = f.add_subplot(211)

ax1.set_title("Firing rate of neuron {0}".format(neuron_id))
ax1.plot(t, firing_rate[neuron_id])
ax1.set_ylabel("r")
ax1.set_xlim([0, xmax])

ax2 = f.add_subplot(212)

ax2.set_title("GABA reverse potential of {0}".format(neuron_id))
ax2.plot(t, gaba_potential[neuron_id])
ax2.set_xlabel("t (# trial)")
ax2.set_ylabel("$E_{GABA}$ mV")
ax2.set_xlim([0, xmax])


plt.show()   