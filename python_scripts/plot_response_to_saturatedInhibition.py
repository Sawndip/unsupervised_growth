# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 13:47:47 2018

@author: jingroup

Script plots response to saturated inhibition
"""
import numpy as np
import matplotlib.pyplot as plt

# model with E_GABA = -51.5 mV
num_spikes_to_saturated = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 6, 7, 8, 9, 10, 10, 11, 11, 11, 11, 11, 11, 10, 10, 9, 9, 8, 6, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
Gsaturated = [float(i) * 0.25 for i in range(len(num_spikes_to_saturated))]
mean_spike_delay_saturated = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 106.82, 44.22, 33, 27.74, 24.6, 22.54, 21.1, 20.1, 19.38, 18.9, 18.6, 18.5, 18.56, 18.8, 19.26, 20, 21.1, 22.78, 25.52, 30.5, 43.28, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
#std_spike_delay_saturated = [-1, -1, 37.5049, 19.9595, 26.2335, 10.385, 8.04711, 9.02335, 6.952, 6.75338, 6.94311, 7.65506, 9.06835, 7.81932, 7.26411, 6.37052, 7.51844, 6.3687, 7.44579, 8.83594, 6.81979, 9.53978, 7.76698, 8.28635, 8.27966, 9.78356, 8.08715, 6.28005, 7.59593, 6.38946, 7.12088, 6.70433, 5.59488, 12.9413, 7.75818, 9.36325, 8.89761, 7.69177, 8.70106, 6.58262, 7.94616, 9.76032]

#TRIAL = 100
#mean_spike_delay_saturated = [mean - TRIAL for mean in mean_spike_delay_saturated]

f = plt.figure()

ax1 = f.add_subplot(211)
#ax1.plot(Gsaturated, num_spikes_to_saturated, label='saturated inhibition')
ax1.plot(Gsaturated, num_spikes_to_saturated)

ax1.set_ylabel('# spikes')
ax1.set_xlim([0,10])

#ax1.legend(loc=4)

ax2 = f.add_subplot(212)
#ax2.plot(Gsaturated, mean_spike_delay_saturated, label='saturated inhibition')
ax2.plot(Gsaturated, mean_spike_delay_saturated)
ax2.set_ylabel('< first spike delay (ms) >')
ax2.set_xlabel('G_inh (mS/cm^2)')
ax2.set_xlim([0,10])
ax2.set_ylim([0, 110])

#ax2.set_xlabel('G_inh (mS/cm^2)')
#ax2.set_ylim([0, max(max(mean_spike_delay_saturated) + 5, max(mean_spike_delay_kicks) + 5)])
#==============================================================================
# ax2.legend(loc=4)
# 
# ax3 = f.add_subplot(313)
# #ax3.plot(Gsaturated, std_spike_delay_saturated, label='saturated inhibition')
# ax3.plot(Gkick_1, std_spike_delay_kicks_1, label=label1)
# ax3.plot(Gkick_2, std_spike_delay_kicks_2, label=label2)
# ax3.plot(Gkick_3, std_spike_delay_kicks_3, label=label3)
# ax3.set_ylabel('<std of spike delays (ms)>')
# ax3.set_xlabel('G_inh (mS/cm^2)')
# ax3.set_xlim([0,120])
# ax3.set_ylim([0,30])
# #ax3.set_ylim([0, max(max(std_spike_delay_saturated) + 5, max(std_spike_delay_kicks) + 5)])
# ax3.legend()
# 
#==============================================================================
#plt.legend()
plt.show()