# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 10:40:25 2017

@author: jingroup

Script analyses neuron states including remodeling, maturation and gaba reverse 
potential
"""

import reading
import matplotlib.pyplot as plt

import numpy as np

filenameStates = "/home/eugene/Output/growthTest/neuron_states.bin"
filenameWeights = "/home/eugene/Output/growthTest/weights_time.bin"

(t, remodeled, gaba_potential, firing_rate) = reading.read_maturation_time_sequence(filenameStates)
(t, weights) = reading.read_synaptic_weights_time_sequence(filenameWeights)

#print t
#print remodeled
#print mature
#print gaba_potential

num_neurons = 20

x = list(t)
x.append(t[-1] + t[1] - t[0])
y = range(num_neurons + 1)

print y

xticklabels = []
xticklabels.extend(t)
#xticklabels.append(t[-1] + (t[1] - t[0]))

X, Y = np.meshgrid(x, y)

plt.figure()

plt.pcolor(X, Y, gaba_potential, cmap='plasma')   
plt.title("Reverse Gaba potential")
    
plt.axis([X.min(), X.max(), Y.min(), Y.max()])
cbar = plt.colorbar()

#plt.yticks(range(len(num_neurons)), target)
plt.xlabel("time (# trial)")
plt.ylabel("neuron id")

plt.figure()

plt.pcolor(X, Y, firing_rate, cmap='plasma', vmin=0, vmax=0.2)   
plt.title("Firing rate")
    
plt.axis([X.min(), X.max(), Y.min(), Y.max()])
cbar = plt.colorbar()

#plt.yticks(range(len(target)), target)
plt.xlabel("time (# trial)")
plt.ylabel("neuron id")

 
f2 = plt.figure()

ax1 = f2.add_subplot(111)
ax1.spy(remodeled)
ax1.set_title("Remodeled")
ax1.set_aspect('auto', 'box')
ax1.invert_yaxis()
ticks =  ax1.get_xticks()
#plt.xticks(ticks[1:-2], xticklabels)
ax1.xaxis.tick_bottom()


#plt.xticks(ticks[1:-2], xticklabels)

#ax2.set_xticklabels(xticklabels)

# plot how specific synaptic weight changes

print weights.shape
print len(t)

f3 = plt.figure()

ax = f3.add_subplot(111)

source_id = 0 # real id of source neuron
target_id = 16 # real id of target neuron

ax.set_title("Synaptic weight {0} -> {1} vs. time".format(source_id, target_id))
ax.plot(t, weights[:,source_id,target_id])
ax.set_xlabel("t (# trial)")
ax.set_ylabel("w")

plt.show()   


#==============================================================================
# f = plt.figure()
# 
# ax1 = f.add_subplot(211)
# ax1.plot(t, num_active)
# ax1.set_ylabel("# active synapses")
# 
# ax2 = f.add_subplot(212)
# ax2.plot(t, num_super)
# ax2.set_ylabel("# super synapses")
# ax2.set_xlabel("t (# trial)")
# 
# plt.show()
# 
# 
#==============================================================================

