# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 14:45:53 2016

@author: jingroup
"""

from reading import *
import matplotlib.pyplot as plt


fileSimInfo = "/home/eugene/Output/sim_info.bin"
fileSynapticInfo = "/home/eugene/Output/synaptic_info.bin"

(trial_duration, synapses_trials_update, weights_trials_update) = read_sim_info(fileSimInfo)

print "trial_duration", trial_duration
print "synapses_trials_update", synapses_trials_update
print "weights_trials_update", weights_trials_update

(num_active, num_super) = read_synaptic_info(fileSynapticInfo)
t_synapses = [i*trial_duration*synapses_trials_update/1000 for i in  xrange(len(num_active))]

print "num_active", num_active

f = plt.figure()

ax1 = f.add_subplot(211)
ax1.plot(t_synapses, num_active)
ax1.set_ylabel("# of active synapses")

ax2 = f.add_subplot(212)
ax2.plot(t_synapses, num_super)
ax2.set_ylabel("# of super synapses")
ax2.set_xlabel("time (s)")

plt.show()