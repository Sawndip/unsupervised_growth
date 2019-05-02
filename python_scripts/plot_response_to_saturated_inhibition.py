# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 13:21:35 2018

@author: jingroup

Script plots burst delay responses to inhibitory kicks
"""
import matplotlib.pyplot as plt

Gstep = 0.25
num_steps = 20

Ginh = [float(i) * Gstep for i in range(num_steps)]

num_spikes_soma = [0, 1.052, 1.21, 0.858, 1.068, 0.898, 1.118, 1.09, 1.062, 1.058, 1.142, 1.186, 1.066, 1.154, 1.206, 1.198, 1.246, 1.076, 1.272, 1.16]

mean_soma_spike_times = [-1, 65.5646, 62.8022, 75.4006, 74.365, 66.9402, 72.7536, 70.6491, 68.3749, 69.4, 69.298, 66.7131, 68.3536, 70.8533, 71.1148, 71.6013, 66.2936, 66.0223, 72.5328, 67.3528]

std_soma_spike_times = [-1, 25.6948, 32.9747, 38.0621, 38.1627, 35.3033, 36.8947, 36.7315, 34.6871, 34.693, 35.2005, 36.3129, 33.9388, 36.4378, 37.1195, 36.0699, 34.5159, 33.8833, 36.2534, 35.9664]

num_spikes_dend = [0, 0.168, 0.14, 0.05, 0.058, 0.024, 0.046, 0.038, 0.028, 0.026, 0.028, 0.028, 0.02, 0.034, 0.034, 0.034, 0.02, 0.016, 0.036, 0.022]

mean_dend_spike_times = [-1, 63.4883, 65.0883, 110.681, 114.333, 117.705, 116.231, 118.088, 117.657, 118.2, 119.467, 108.376, 119.066, 119.467, 123.834, 121.572, 120.888, 121.142, 121.892, 121.591]

std_dend_spike_times = [-1, 25.5779, 34.4416, 18.8241, 5.55172, 4.60074, 4.43395, 5.53474, 4.3888, 3.07292, 5.79858, 38.8119, 5.32339, 4.98347, 5.44741, 2.79926, 3.44856, 3.47444, 4.77657, 2.98761]



plt.figure()
plt.plot(Ginh, num_spikes_soma)
plt.xlabel('Gi (mS/cm^2)')
plt.ylabel('<# of somatic spikes>')

plt.figure()
plt.plot(Ginh, mean_soma_spike_times)
plt.xlabel('Gi (mS/cm^2)')
plt.ylabel('<soma spike time> (ms)')

plt.figure()
plt.plot(Ginh, std_soma_spike_times)
plt.xlabel('Gi (mS/cm^2)')
plt.ylabel('<std of soma spike time (ms)')

plt.figure()
plt.plot(Ginh, num_spikes_dend)
plt.xlabel('Gi (mS/cm^2)')
plt.ylabel('<# of dendritic spikes>')

plt.figure()
plt.plot(Ginh, mean_dend_spike_times)
plt.xlabel('Gi (mS/cm^2)')
plt.ylabel('<dend spike time> (ms)')

plt.figure()
plt.plot(Ginh, std_dend_spike_times)
plt.xlabel('Gi (mS/cm^2)')
plt.ylabel('<std of dend spike time (ms)')

plt.show()
