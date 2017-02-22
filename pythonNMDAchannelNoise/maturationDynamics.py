# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 11:17:52 2017

@author: jingroup

Script contasts firing rate and GABA reverse potential
"""
import reading
import matplotlib.pyplot as plt

filenameStates = "/home/eugene/Output/networks/test170217/maturation_time_sequence.bin"

(target, t, remodeled, mature, gaba_potential, firing_rate) = reading.read_maturation_time_sequence(filenameStates)

print target

neuron_id = 19 # neuron id which state to display

f1 = plt.figure()

plt.suptitle("Maturation for neuron {0}".format(neuron_id))

ax1 = f1.add_subplot(211)

ax1.plot(t, firing_rate[neuron_id])
ax1.set_ylabel("r")
#ax1.set_xlim([0, 2000])

ax2 = f1.add_subplot(212)
ax2.plot(t, gaba_potential[neuron_id])

ax2.set_ylabel("$E_{GABA}$ (mV)")
ax2.set_xlabel("time (# trial)")
#ax2.set_xlim([0, 2000])

plt.show()