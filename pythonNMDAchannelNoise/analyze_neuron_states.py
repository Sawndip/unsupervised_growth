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

filenameStates = "/home/eugene/Output/networks/test230217/maturation_time_sequence.bin"
filenameWeights = "/home/eugene/Output/networks/test230217/weightsTimeSequence.bin"

(target, t, remodeled, mature, gaba_potential, firing_rate) = reading.read_maturation_time_sequence(filenameStates)
(source, target, t, weights) = reading.read_synaptic_weights_time_sequence(filenameWeights)

print source
print target
#print weights
#print t
#print remodeled
#print mature
#print gaba_potential

num_target = len(target)
num_source = len(source)

x = list(t)
x.append(t[-1] + t[1] - t[0])
y = range(num_target + 1)

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

plt.yticks(range(len(target)), target)
plt.xlabel("time (# trial)")
plt.ylabel("neuron id")

plt.figure()

plt.pcolor(X, Y, firing_rate, cmap='plasma', vmin=0, vmax=1)   
plt.title("Firing rate")
    
plt.axis([X.min(), X.max(), Y.min(), Y.max()])
cbar = plt.colorbar()

plt.yticks(range(len(target)), target)
plt.xlabel("time (# trial)")
plt.ylabel("neuron id")

 
f2 = plt.figure()

ax1 = f2.add_subplot(121)
ax1.spy(remodeled)
ax1.set_title("Remodeled")
ax1.set_aspect('auto', 'box')
ax1.invert_yaxis()
ticks =  ax1.get_xticks()
#plt.xticks(ticks[1:-2], xticklabels)
ax1.xaxis.tick_bottom()

ax2 = f2.add_subplot(122)
ax2.spy(mature)
ax2.set_title("Mature")
ax2.set_aspect('auto', 'box')
ax2.invert_yaxis()
ax2.xaxis.tick_bottom()

#plt.xticks(ticks[1:-2], xticklabels)

#ax2.set_xticklabels(xticklabels)

# plot how specific synaptic weight changes

f3 = plt.figure()

ax = f3.add_subplot(111)

source_neuron = 0 # real id of source neuron
target_neuron = 25 # real id of target neuron

ind_source = source.index(source_neuron) # index of source neuron in the list
ind_target = target.index(target_neuron) # index of target neuron in the list


ax.set_title("Synaptic weight {0} -> {1} vs. time".format(source_neuron, target_neuron))
ax.plot(t, weights[ind_source][ind_target])
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
