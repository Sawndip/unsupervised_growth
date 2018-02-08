# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 10:01:57 2018

@author: jingroup

Script plots responses of neurons to different inhibitory inputs
"""
import matplotlib.pyplot as plt

# model with E_L = -80 mV
Gkick = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5]
num_spikes_to_kicks = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.002, 0.006, 0.004, 0.01, 0.036, 0.028, 0.056, 0.066, 0.088]
mean_spike_delay_kicks = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 9.96, 17.2533, 14.93, 14.948, 14.0411, 13.5014, 13.3414, 13.737, 13.3477]
std_spike_delay_kicks = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1.50098, 3.6911, 4.69041, 3.02865, 1.57419, 2.83242, 2.79134, 2.63214]

# model with E_L = -65 mV
#Gkick = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7]
#num_spikes_to_kicks = [0, 0, 0, 0, 0, 0.004, 0.008, 0.02, 0.022, 0.056, 0.058, 0.048, 0.072, 0.064, 0.08]
#mean_spike_delay_kicks = [-1, -1, -1, -1, -1, 43.84, 34.02, 30.742, 37.44, 32.4043, 31.1545, 33.4367, 32.1522, 32.0481, 30.663]
#std_spike_delay_kicks = [-1, -1, -1, -1, -1, 25.1164, 8.68167, 8.01545, 11.065, 7.14262, 7.15878, 8.60652, 8.90736, 7.50936, 7.09437]

# model with E_L = -80 mV
num_spikes_to_saturated = [0, 0, 0, 0, 0, 0, 0, 0.006, 0.04, 0.054, 0.074, 0.09, 0.054, 0.056, 0.062, 0.058, 0.074, 0.068, 0.094, 0.082, 0.06, 0.078, 0.058, 0.064, 0.07, 0.06, 0.064, 0.054, 0.086, 0.07, 0.052, 0.068, 0.06, 0.072, 0.08, 0.08, 0.058, 0.078, 0.068, 0.084, 0.05, 0.076]
Gsaturated = [float(i) * 0.25 for i in range(len(num_spikes_to_saturated))]
mean_spike_delay_saturated = [-1, -1, -1, -1, -1, -1, -1, 90.2467, 60.035, 69.9778, 78.0649, 81.7689, 92.3459, 99.3743, 107.737, 109.11, 110.194, 109.662, 110.4, 110.641, 109.879, 110.328, 111.024, 110.879, 111.23, 111.904, 113.106, 111.175, 111.776, 112.486, 112.109, 112.533, 112.548, 112.712, 113.124, 112.817, 112.523, 112.751, 113.728, 113.272, 113.541, 113.975]
std_spike_delay_saturated = [-1, -1, -1, -1, -1, -1, -1, 27.0691, 27.8375, 30.0645, 29.4894, 29.9642, 24.4028, 22.5498, 5.27811, 3.70505, 3.36845, 3.11654, 2.85402, 3.2277, 2.46281, 2.87244, 2.68974, 2.55255, 2.85803, 2.95038, 3.02457, 2.40917, 3.29368, 2.79898, 2.85102, 2.4733, 3.0345, 3.28651, 2.46024, 2.62871, 2.81148, 3.07703, 3.16802, 2.51587, 2.41487, 2.80641]

# model with E_L = -65 mV
#num_spikes_to_saturated = [0, 0, 0.00666667, 0.0266667, 0.0333333, 0.04, 0.03, 0.0566667, 0.0533333, 0.0533333, 0.06, 0.0766667, 0.0533333, 0.0866667, 0.0766667, 0.0533333, 0.0766667, 0.0733333, 0.0866667, 0.06, 0.09, 0.06, 0.07, 0.0733333, 0.0533333, 0.07, 0.06, 0.0633333, 0.0566667, 0.0466667, 0.0766667, 0.0633333, 0.09, 0.0433333, 0.0633333, 0.0633333, 0.0633333, 0.0466667, 0.0733333, 0.0466667, 0.0566667, 0.0966667]
#Gsaturated = [float(i) * 0.25 for i in range(len(num_spikes_to_saturated))]
#mean_spike_delay_saturated = [-1, -1, 98.46, 101.475, 115.678, 126.913, 121.353, 125.145, 125.421, 124.465, 125.361, 127.91, 127.351, 126.584, 127.276, 127.451, 127.621, 127.677, 131.506, 130.407, 129.633, 130.171, 131.564, 131.469, 129.926, 130.458, 132.851, 130.854, 130.079, 130.033, 132.461, 132.537, 131.266, 133.372, 133.34, 132.578, 133.254, 133.644, 134.405, 135.503, 134.635, 135.37]
#std_spike_delay_saturated = [-1, -1, 37.5049, 19.9595, 26.2335, 10.385, 8.04711, 9.02335, 6.952, 6.75338, 6.94311, 7.65506, 9.06835, 7.81932, 7.26411, 6.37052, 7.51844, 6.3687, 7.44579, 8.83594, 6.81979, 9.53978, 7.76698, 8.28635, 8.27966, 9.78356, 8.08715, 6.28005, 7.59593, 6.38946, 7.12088, 6.70433, 5.59488, 12.9413, 7.75818, 9.36325, 8.89761, 7.69177, 8.70106, 6.58262, 7.94616, 9.76032]

TRIAL = 100
mean_spike_delay_saturated = [mean - TRIAL for mean in mean_spike_delay_saturated]

f = plt.figure()

ax1 = f.add_subplot(311)
ax1.plot(Gsaturated, num_spikes_to_saturated, label='saturated inhibition')
ax1.plot(Gkick, num_spikes_to_kicks, label='inhibitory kick')
ax1.set_ylabel('<# of somatic spikes>')
ax1.set_xlim([0,7])

ax1.legend()

ax2 = f.add_subplot(312)
ax2.plot(Gsaturated, mean_spike_delay_saturated, label='saturated inhibition')
ax2.plot(Gkick, mean_spike_delay_kicks, label='inhibitory kick')
ax2.set_ylabel('<spike delay (ms)>')
ax2.set_xlim([0,7])
ax2.set_ylim([0, max(max(mean_spike_delay_saturated) + 5, max(mean_spike_delay_kicks) + 5)])
ax2.legend()

ax3 = f.add_subplot(313)
ax3.plot(Gsaturated, std_spike_delay_saturated, label='saturated inhibition')
ax3.plot(Gkick, std_spike_delay_kicks, label='inhibitory kick')
ax3.set_ylabel('<std of spike delays (ms)>')
ax3.set_xlabel('G_inh (mS/cm^2)')
ax3.set_xlim([0,10])
ax3.set_ylim([0, max(max(std_spike_delay_saturated) + 5, max(std_spike_delay_kicks) + 5)])
ax3.legend()

plt.legend()
plt.show()