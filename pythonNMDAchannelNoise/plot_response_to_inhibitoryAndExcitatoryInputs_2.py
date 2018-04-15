# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 15:23:25 2018

@author: jingroup

Script plots probability to spike and spike delays of HVC-RA neurons
"""
import matplotlib.pyplot as plt
import numpy as np

CONVERTION_CONSTANT = 10.0
num_targets = 10
#==============================================================================
# 
# ##### Comparison between mature and immature models #####
# #synaptic_weight_immature = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5]) * CONVERTION_CONSTANT
# #probability_to_spike_immature = [0, 0.01, 0.364, 0.938, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# #num_spikes_immature = np.array([0, 5, 182, 469, 501, 538, 670, 911, 993, 999, 1000, 1001, 1008, 1048, 1224, 1403, 1484, 1500, 1500, 1500]) / float(num_targets)
# #first_spike_delay_immature = [-1, 34.956, 26.0764, 18.484, 13.3945, 11.1331, 10.0219, 9.19688, 8.61476, 8.15608, 7.80652, 7.5028, 7.23292, 7.02664, 6.79252, 6.5942, 6.42956, 6.27428, 6.1336, 6.016]
# 
# synaptic_weight_immature = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5]) * CONVERTION_CONSTANT
# probability_to_spike_immature = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_spikes_immature = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 10, 20, 20, 20, 20, 20, 30, 30]) / float(num_targets)
# first_spike_delay_immature = [-1, -1, -1, -1, -1, -1, -1, -1, -1, 13.78, 9.94, 8.52, 7.66, 7.08, 6.62, 6.28, 5.98, 5.74, 5.54, 5.36]
# 
# 
# #synaptic_weight_mature = np.array([0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950])* CONVERTION_CONSTANT
# #probability_to_spike_mature = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# #num_spikes_mature = np.array([0, 0, 0, 0, 0, 529, 3193, 2542, 2500, 2500, 2011, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000]) / float(num_targets)
# #first_spike_delay_mature = [-1, -1, -1, -1, -1, 11.7548, 9.15028, 7.71468, 6.79788, 6.40596, 6.14052, 5.97544, 5.82764, 5.73208, 5.6288, 5.54756, 5.47104, 5.4022, 5.34828, 5.30048]
# 
# synaptic_weight_mature = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975])* CONVERTION_CONSTANT
# probability_to_spike_mature = [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_spikes_mature = np.array([0, 0, 0, 20, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40]) / float(num_targets)
# first_spike_delay_mature = [-1, -1, -1, 12.44, 8.72, 7.48, 6.72, 6.18, 5.78, 5.46, 5.18, 4.94, 4.74, 4.56, 4.42, 4.28, 4.16, 4.06, 3.96, 3.88, 3.8, 3.72, 3.64, 3.58, 3.5, 3.46, 3.4, 3.34, 3.3, 3.26, 3.22, 3.18, 3.14, 3.1, 3.08, 3.04, 3.02, 2.98, 2.96, 2.94]
# 
# 
# f = plt.figure()
# 
# ax1 = f.add_subplot(311)
# ax1.plot(synaptic_weight_immature, probability_to_spike_immature, label='immature')
# #ax1.plot(synaptic_weight_mature, probability_to_spike_mature, label='mature')
# 
# ax1.set_ylabel('probability to spike')
# ax1.set_ylim([0, 1.1])
# 
# ax2 = f.add_subplot(312)
# ax2.plot(synaptic_weight_immature, num_spikes_immature, label='immature')
# #ax2.plot(synaptic_weight_mature, num_spikes_mature, label='mature')
# ax2.set_ylabel('# of spikes')
# ax2.set_ylim([0, max(num_spikes_immature) + 0.5])
# 
# ax3 = f.add_subplot(313)
# ax3.plot(synaptic_weight_immature, first_spike_delay_immature, label='immature')
# ax3.set_ylim([0, max(first_spike_delay_immature) + 0.5])
# #ax3.plot(synaptic_weight_mature, first_spike_delay_mature, label='mature')
# ax3.set_ylabel('delay to first spike (ms)')
# ax3.set_xlabel('Single input weight (pS)')
# 
# 
# f = plt.figure()
# 
# ax1 = f.add_subplot(311)
# #ax1.plot(synaptic_weight_immature, probability_to_spike_immature, label='immature')
# ax1.plot(synaptic_weight_mature, probability_to_spike_mature, label='mature')
# 
# ax1.set_ylabel('probability to spike')
# ax1.set_ylim([0, 1.1])
# 
# ax2 = f.add_subplot(312)
# #ax2.plot(synaptic_weight_immature, num_spikes_immature, label='immature')
# ax2.plot(synaptic_weight_mature, num_spikes_mature, label='mature')
# ax2.set_ylabel('# of spikes')
# ax2.set_ylim([0, max(num_spikes_mature) + 0.5])
# 
# ax3 = f.add_subplot(313)
# #ax3.plot(synaptic_weight_immature, first_spike_delay_immature, label='immature')
# ax3.plot(synaptic_weight_mature, first_spike_delay_mature, label='mature')
# ax3.set_ylim([0, max(first_spike_delay_immature) + 0.5])
# 
# ax3.set_ylabel('delay to first spike (ms)')
# ax3.set_xlabel('Single input weight (pS)')
# 
#==============================================================================

#==============================================================================
# 
# ##### Comparison between different inhibition strength for immature model #####
# num_targets = 1000
# 
# synaptic_weight_kick0 = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5]) * CONVERTION_CONSTANT
# probability_to_spike_kick0 = [0, 0, 0.001, 0.004, 0.023, 0.048, 0.132, 0.277, 0.456, 0.666, 0.81, 0.924, 0.977, 0.991, 0.998, 0.999, 1, 1, 1, 1]
# num_spikes_kick0 = np.array([0, 0, 1, 4, 23, 50, 137, 290, 506, 781, 1012, 1311, 1577, 1761, 1949, 2210, 2393, 2568, 2773, 2946]) / float(num_targets)
# mean_first_spike_delay_kick0 = [-1, -1, 21.36, 13.9, 13.8609, 12.9025, 12.55, 12.089, 11.5159, 10.504, 9.74862, 8.74273, 8.00154, 7.30521, 6.85637, 6.28142, 6.00324, 5.90964, 5.65382, 5.49476]
# std_first_spike_delay_kick0 = [0, 0, -1, 3.99276, 2.93644, 3.34343, 3.76677, 3.61812, 3.23429, 3.02505, 2.81422, 2.48489, 1.8374, 1.4787, 1.25725, 2.38737, 2.21516, 0.71554, 0.661401, 0.61508]
# 
# synaptic_weight_kick5 = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5]) * CONVERTION_CONSTANT
# probability_to_spike_kick5 = [0.035, 0.08, 0.108, 0.197, 0.314, 0.424, 0.603, 0.718, 0.799, 0.912, 0.952, 0.982, 0.994, 0.998, 0.999, 1, 1, 1, 1, 1]
# num_spikes_kick5 = np.array([36, 82, 111, 199, 324, 441, 637, 792, 927, 1140, 1276, 1479, 1660, 1770, 1922, 2118, 2240, 2402, 2543, 2679]) / float(num_targets)
# mean_first_spike_delay_kick5 = [15.6451, 16.5952, 15.7928, 13.8654, 13.4664, 12.5637, 11.7178, 11.0731, 10.6471, 9.44662, 8.7412, 8.10287, 7.69394, 7.28918, 6.91161, 6.46382, 6.28844, 6.13906, 5.95692, 5.797]
# std_first_spike_delay_kick5 = [3.96409, 4.73036, 4.3252, 4.63187, 3.88098, 3.74575, 3.41557, 2.97771, 3.20348, 2.36521, 1.80009, 1.68938, 1.26554, 1.25101, 0.984933, 2.3405, 2.19598, 0.631381, 0.635124, 0.56638]
# 
# synaptic_weight_kick10 = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5]) * CONVERTION_CONSTANT
# probability_to_spike_kick10 = [0.09, 0.146, 0.169, 0.244, 0.345, 0.441, 0.555, 0.642, 0.724, 0.853, 0.907, 0.951, 0.98, 0.987, 0.991, 0.997, 0.999, 1, 1, 1]
# num_spikes_kick10 = np.array([92, 150, 173, 247, 356, 457, 587, 709, 831, 1027, 1133, 1319, 1472, 1564, 1694, 1894, 1986, 2173, 2272, 2393]) / float(num_targets)
# mean_first_spike_delay_kick10 = [17.0424, 17.0956, 16.1182, 14.5943, 14.0525, 13.7825, 12.7806, 12.3448, 12.0167, 10.8732, 10.2099, 9.36412, 8.90937, 8.48219, 7.85045, 7.30299, 7.10539, 6.8163, 6.64284, 6.39374]
# std_first_spike_delay_kick10 = [5.29658, 4.99667, 4.90885, 4.65403, 3.85821, 4.17085, 3.85236, 3.46809, 3.65119, 3.45149, 2.90326, 2.44927, 2.05641, 2.07459, 1.40964, 2.56961, 2.40508, 0.887006, 0.886433, 0.73157]
# 
# synaptic_weight_kick15 = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5]) * CONVERTION_CONSTANT
# probability_to_spike_kick15 = [0.111, 0.152, 0.161, 0.213, 0.286, 0.325, 0.427, 0.479, 0.557, 0.686, 0.762, 0.816, 0.898, 0.924, 0.945, 0.973, 0.982, 0.997, 0.998, 0.999]
# num_spikes_kick15 = np.array([113, 156, 166, 214, 290, 332, 447, 513, 611, 774, 882, 1002, 1164, 1244, 1345, 1546, 1628, 1815, 1898, 2035]) / float(num_targets)
# mean_first_spike_delay_kick15 = [19.3514, 18.717, 18.1748, 16.6409, 16.0567, 15.7346, 15.0814, 14.8162, 14.351, 13.0669, 12.6803, 11.723, 11.2309, 10.8823, 9.89306, 9.24687, 8.93646, 8.28269, 8.09894, 7.55479]
# std_first_spike_delay_kick15 = [5.61841, 4.83567, 5.2619, 5.01571, 4.33129, 4.5107, 4.24223, 4.20581, 4.24715, 3.89294, 3.77323, 3.52318, 3.06115, 3.38339, 2.57247, 3.39704, 3.13746, 1.83473, 1.70325, 1.31319]
# 
# synaptic_weight_kick20 = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5]) * CONVERTION_CONSTANT
# probability_to_spike_kick20 = [0.109, 0.16, 0.138, 0.174, 0.244, 0.25, 0.308, 0.338, 0.428, 0.524, 0.558, 0.634, 0.736, 0.746, 0.803, 0.882, 0.886, 0.948, 0.963, 0.972]
# num_spikes_kick20 = np.array([110, 165, 144, 176, 247, 252, 318, 352, 455, 574, 616, 714, 851, 880, 971, 1158, 1181, 1408, 1461, 1601]) / float(num_targets)
# mean_first_spike_delay_kick20 = [20.2809, 19.8833, 19.5212, 18.7182, 17.9441, 17.6401, 17.3367, 17.0398, 16.7198, 15.6511, 15.0089, 14.2818, 14.2542, 13.6915, 12.8528, 12.2533, 12.0409, 10.9039, 10.8458, 9.9201]
# std_first_spike_delay_kick20 = [4.85372, 4.84793, 4.70903, 5.04526, 5.12962, 4.50171, 4.36873, 4.39459, 4.36824, 4.42888, 4.1105, 4.36406, 3.9416, 3.90286, 3.75719, 4.46942, 4.26588, 3.06411, 3.11527, 2.73964]
# 
# synaptic_weight_kick25 = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5]) * CONVERTION_CONSTANT
# probability_to_spike_kick25 = [ 0.112, 0.164, 0.139, 0.157, 0.209, 0.204, 0.242, 0.276, 0.331, 0.4, 0.438, 0.494, 0.554, 0.602, 0.631, 0.727, 0.74, 0.825, 0.84, 0.871]
# num_spikes_kick25 = np.array([112, 169, 142, 159, 212, 205, 251, 289, 346, 429, 460, 541, 610, 665, 699, 835, 864, 1015, 1065, 1186]) / float(num_targets)
# mean_first_spike_delay_kick25 = [21.5659, 20.8666, 20.5814, 20.0243, 19.3953, 18.9109, 18.9307, 18.672, 18.3742, 17.8393, 16.9702, 16.4642, 16.4319, 16.2167, 15.5473, 15.0042, 14.8901, 14.0278, 13.799, 12.7961]
# std_first_spike_delay_kick25 = [5.49766, 4.69547, 4.35992, 4.86284, 5.18471, 4.6375, 4.32503, 5.7807, 4.53998, 4.81248, 4.30459, 4.63778, 4.01599, 4.32283, 4.38144, 5.00093, 4.97303, 3.81253, 3.9233, 3.83485]
# 
# synaptic_weight_kick30 = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5]) * CONVERTION_CONSTANT
# probability_to_spike_kick30 = [0.1, 0.154, 0.129, 0.144, 0.184, 0.173, 0.207, 0.232, 0.289, 0.332, 0.365, 0.395, 0.449, 0.474, 0.503, 0.586, 0.596, 0.658, 0.698, 0.75]
# num_spikes_kick30 = np.array([100, 158, 132, 147, 189, 174, 214, 241, 301, 347, 377, 430, 482, 506, 536, 647, 649, 757, 802, 907]) / float(num_targets)
# mean_first_spike_delay_kick30 = [22.4788, 21.5929, 21.2998, 21.4625, 20.5282, 20.1632, 20.1042, 19.9489, 19.8777, 18.9143, 18.7546, 18.1044, 17.8837, 18.0343, 17.583, 16.9617, 17.2068, 16.3179, 16.0732, 15.1621]
# std_first_spike_delay_kick30 = [5.10509, 4.55224, 4.16508, 5.73504, 5.21771, 4.56061, 4.11033, 4.50235, 4.83798, 4.58248, 4.84206, 4.40437, 4.10876, 4.4218, 4.44448, 5.4314, 5.53317, 4.09034, 4.19795, 4.26035]
# 
# synaptic_weight_kick50 = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5]) * CONVERTION_CONSTANT
# probability_to_spike_kick50 = [0.093, 0.139, 0.126, 0.117, 0.146, 0.12, 0.157, 0.158, 0.191, 0.228, 0.223, 0.235, 0.247, 0.26, 0.28, 0.323, 0.327, 0.352, 0.381, 0.387]
# num_spikes_kick50 = np.array([94, 141, 130, 120, 147, 121, 158, 166, 198, 232, 228, 246, 258, 270, 289, 334, 336, 377, 399, 418]) / float(num_targets)
# mean_first_spike_delay_kick50 = [24.6587, 24.1666, 24.1179, 24.5932, 24.0507, 23.5205, 23.5918, 23.9299, 23.5266, 23.1304, 22.9808, 22.4196, 22.0092, 22.3553, 21.8514, 21.6154, 21.6097, 21.1893, 21.1946, 21.0812]
# std_first_spike_delay_kick50 = [4.88623, 4.88954, 5.03754, 5.79506, 5.62079, 4.83928, 5.34716, 4.35618, 4.75031, 4.64649, 4.41691, 4.43706, 4.66359, 4.79874, 4.50396, 6.45694, 6.72855, 4.13652, 4.7921, 5.22634]
# 
# 
# 
# f = plt.figure()
# 
# ax1 = f.add_subplot(311)
# ax1.plot(synaptic_weight_kick0, probability_to_spike_kick0, label='Ginh = 0 ms/cm^2')
# ax1.plot(synaptic_weight_kick5, probability_to_spike_kick5, label='Ginh = 5 ms/cm^2')
# ax1.plot(synaptic_weight_kick10, probability_to_spike_kick10, label='Ginh = 10 ms/cm^2')
# ax1.plot(synaptic_weight_kick15, probability_to_spike_kick15, label='Ginh = 15 ms/cm^2')
# ax1.plot(synaptic_weight_kick20, probability_to_spike_kick20, label='Ginh = 20 ms/cm^2')
# ax1.plot(synaptic_weight_kick25, probability_to_spike_kick25, label='Ginh = 25 ms/cm^2')
# ax1.plot(synaptic_weight_kick30, probability_to_spike_kick30, label='Ginh = 30 ms/cm^2')
# ax1.plot(synaptic_weight_kick50, probability_to_spike_kick50, label='Ginh = 50 ms/cm^2')
# ax1.set_ylabel('probability to spike')
# ax1.set_ylim([0, 1.1])
# ax1.legend(loc=4)
# 
# ax2 = f.add_subplot(312)
# ax2.plot(synaptic_weight_kick0, num_spikes_kick0, label='Ginh = 0 ms/cm^2')
# ax2.plot(synaptic_weight_kick5, num_spikes_kick5, label='Ginh = 5 ms/cm^2')
# ax2.plot(synaptic_weight_kick10, num_spikes_kick10, label='Ginh = 10 ms/cm^2')
# ax2.plot(synaptic_weight_kick15, num_spikes_kick15, label='Ginh = 15 ms/cm^2')
# ax2.plot(synaptic_weight_kick20, num_spikes_kick20, label='Ginh = 20 ms/cm^2')
# ax2.plot(synaptic_weight_kick25, num_spikes_kick25, label='Ginh = 25 ms/cm^2')
# ax2.plot(synaptic_weight_kick30, num_spikes_kick30, label='Ginh = 30 ms/cm^2')
# ax2.plot(synaptic_weight_kick50, num_spikes_kick50, label='Ginh = 50 ms/cm^2')
# ax2.set_ylabel('# of spikes')
# ax2.set_ylim([0, 4])
# ax2.legend(loc=2)
# 
# ax3 = f.add_subplot(313)
# ax3.plot(synaptic_weight_kick0, mean_first_spike_delay_kick0, label='Ginh = 0 ms/cm^2')
# ax3.plot(synaptic_weight_kick5, mean_first_spike_delay_kick5, label='Ginh = 5 ms/cm^2')
# ax3.plot(synaptic_weight_kick10, mean_first_spike_delay_kick10, label='Ginh = 10 ms/cm^2')
# ax3.plot(synaptic_weight_kick15, mean_first_spike_delay_kick15, label='Ginh = 15 ms/cm^2')
# ax3.plot(synaptic_weight_kick20, mean_first_spike_delay_kick20, label='Ginh = 20 ms/cm^2')
# ax3.plot(synaptic_weight_kick25, mean_first_spike_delay_kick25, label='Ginh = 25 ms/cm^2')
# ax3.plot(synaptic_weight_kick30, mean_first_spike_delay_kick30, label='Ginh = 30 ms/cm^2')
# ax3.plot(synaptic_weight_kick50, mean_first_spike_delay_kick50, label='Ginh = 50 ms/cm^2')
# ax3.set_ylim([0, 25])
# #ax3.plot(synaptic_weight_mature, first_spike_delay_mature, label='mature')
# ax3.set_ylabel('delay to first spike (ms)')
# ax3.set_xlabel('Single input weight (pS)')
# ax3.legend()
#==============================================================================



#==============================================================================
# ##### Comparison between different inhibition strength for mature model #####
# num_targets = 100
# 
# synaptic_weight_kick0 = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975]) * CONVERTION_CONSTANT
# probability_to_spike_kick0 = [0, 0, 0, 0.88, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_spikes_kick0 = np.array([0, 0, 0, 324, 519, 501, 500, 498, 497, 493, 488, 484, 486, 474, 478, 468, 469, 458, 459, 467, 457, 447, 447, 439, 448, 446, 437, 432, 440, 429, 438, 427, 425, 428, 424, 423, 424, 424, 420, 412]) / float(num_targets)
# mean_first_spike_delay_kick0 = [-1, -1, -1, 11.9861, 8.9696, 7.6454, 6.8626, 6.3368, 5.9144, 5.5404, 5.322, 5.0694, 4.851, 4.661, 4.4882, 4.4224, 4.2748, 4.2466, 4.1256, 4.0464, 3.8732, 3.8676, 3.742, 3.7448, 3.594, 3.6046, 3.6178, 3.4962, 3.4898, 3.388, 3.2974, 3.2494, 3.2634, 3.3676, 3.219, 3.184, 3.1078, 3.1084, 3.0806, 3.0418]
# std_first_spike_delay_kick0 = [0, 0, 0, 1.69543, 0.771096, 0.445048, 0.376777, 0.264028, 0.233205, 0.20207, 0.168739, 0.177194, 0.144875, 0.157849, 0.138625, 0.139444, 0.134047, 0.143831, 0.132972, 0.121867, 0.121613, 0.114443, 0.112115, 0.121974, 0.118799, 0.103469, 0.0971054, 0.116643, 0.107327, 0.0920694, 0.108643, 0.101959, 0.100787, 0.0870982, 0.0939106, 0.0930732, 0.104983, 0.10031, 0.0838688, 0.094499]
# 
# synaptic_weight_kick0_5 = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975]) * CONVERTION_CONSTANT
# probability_to_spike_kick0_5 = [ 0, 0, 0, 0, 0, 0.13, 0.89, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_spikes_kick0_5 = np.array([0, 0, 0, 0, 0, 13, 137, 569, 589, 542, 507, 500, 500, 499, 499, 481, 485, 459, 453, 468, 458, 445, 443, 430, 439, 437, 429, 429, 429, 418, 431, 425, 422, 420, 418, 418, 416, 419, 417, 409]) / float(num_targets)
# mean_first_spike_delay_kick0_5 = [-1, -1, -1, -1, -1, 12.7369, 11.3939, 9.3648, 7.9058, 7.0604, 6.5968, 5.8166, 5.525, 5.1874, 4.8282, 4.514, 4.3846, 4.268, 4.1354, 4.063, 3.8934, 3.8776, 3.7468, 3.7536, 3.5954, 3.6054, 3.6184, 3.496, 3.4894, 3.3878, 3.2972, 3.2492, 3.2628, 3.3676, 3.2186, 3.184, 3.1076, 3.1084, 3.0806, 3.0418]
# std_first_spike_delay_kick0_5 = [0, 0, 0, 0, 0, 2.45074, 2.15294, 1.25933, 0.821432, 0.574755, 0.450613, 0.431601, 0.343322, 0.34037, 0.254083, 0.182297, 0.179986, 0.155076, 0.140529, 0.131268, 0.132606, 0.12011, 0.115022, 0.126615, 0.120398, 0.105289, 0.098706, 0.116792, 0.106986, 0.0921963, 0.108768, 0.101957, 0.100786, 0.0870982, 0.0937762, 0.0930732, 0.104786, 0.10031, 0.0838688, 0.094499]
# 
# synaptic_weight_kick1_0 = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975]) * CONVERTION_CONSTANT
# probability_to_spike_kick1_0 = [0, 0, 0, 0, 0, 0, 0, 0.02, 0.36, 0.74, 0.98, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_spikes_kick1_0 = np.array([0, 0, 0, 0, 0, 0, 0, 2, 37, 81, 155, 239, 646, 614, 572, 495, 499, 462, 452, 470, 460, 443, 440, 424, 433, 431, 420, 419, 423, 412, 420, 418, 411, 410, 416, 410, 413, 413, 410, 407]) / float(num_targets)
# mean_first_spike_delay_kick1_0 = [-1, -1, -1, -1, -1, -1, -1, 17.48, 12.4933, 11.2541, 9.8151, 7.555, 6.8848, 6.3896, 5.6248, 4.7044, 4.6126, 4.2954, 4.1462, 4.0816, 3.9204, 3.8888, 3.7518, 3.761, 3.5974, 3.6058, 3.6188, 3.496, 3.4894, 3.3874, 3.297, 3.2492, 3.2628, 3.3676, 3.2182, 3.184, 3.1076, 3.1084, 3.0806, 3.0418]
# std_first_spike_delay_kick1_0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7.34792, 2.19686, 2.97049, 1.71367, 0.071485, 0.080838, 0.0895339, 0.0981032, 0.0765087, 0.0698732, 0.0584998, 0.062085, 0.0571693, 0.0564055, 0.049353, 0.0432661, 0.0395219, 0.0354429, 0.0313024, 0.0315082, 0.0308656, 0.0307081, 0.0309813, 0.0264071, 0.0327457, 0.0320826, 0.0315242, 0.0328449, 0.027634]
# 
# synaptic_weight_kick1_5 = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975]) * CONVERTION_CONSTANT
# probability_to_spike_kick1_5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0.07, 0.38, 0.79, 0.94, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_spikes_kick1_5 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 7, 38, 79, 127, 170, 466, 364, 460, 454, 473, 468, 440, 438, 421, 428, 423, 416, 414, 417, 409, 415, 410, 407, 408, 414, 407, 412, 410, 407, 406]) / float(num_targets)
# mean_first_spike_delay_kick1_5 = [-1, -1, -1, -1, -1, -1, -1, -1, -1, 14.78, 13.34, 10.14, 11.2691, 8.96319, 7.3932, 4.9198, 5.1402, 4.2996, 4.158, 4.106, 3.9602, 3.9002, 3.7566, 3.7694, 3.5992, 3.6062, 3.6188, 3.496, 3.4892, 3.387, 3.2962, 3.249, 3.2628, 3.3676, 3.2182, 3.184, 3.1076, 3.1082, 3.0806, 3.0416]
# std_first_spike_delay_kick1_5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 2.91255, 2.66104, 3.09248, 2.50976, 1.74585, 0.592243, 0.85929, 0.191179, 0.156515, 0.161445, 0.174252, 0.131303, 0.121408, 0.137259, 0.124557, 0.108122, 0.101058, 0.116792, 0.107061, 0.0922174, 0.108867, 0.101916, 0.100786, 0.0870982, 0.0937262, 0.0930732, 0.104786, 0.100266, 0.0838688, 0.0946095]
# 
# synaptic_weight_kick2_0 = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975]) * CONVERTION_CONSTANT
# probability_to_spike_kick2_0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.02, 0, 0.25, 0.53, 0.94, 0.93, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_spikes_kick2_0 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 25, 54, 165, 146, 459, 456, 479, 482, 442, 433, 421, 427, 420, 413, 409, 412, 405, 411, 407, 405, 407, 409, 406, 408, 407, 402, 401]) / float(num_targets)
# mean_first_spike_delay_kick2_0 = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 6.1, -1, 10.0256, 10.2958, 5.7634, 5.92882, 4.3106, 4.1738, 4.1424, 4.025, 3.9138, 3.7618, 3.7798, 3.6002, 3.606, 3.619, 3.4958, 3.489, 3.3868, 3.2958, 3.2488, 3.2628, 3.3674, 3.2178, 3.184, 3.1072, 3.1082, 3.0806, 3.0416]
# std_first_spike_delay_kick2_0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.79196, 0, 2.78013, 3.48548, 1.97782, 1.65937, 0.214023, 0.168456, 0.197648, 0.239383, 0.140467, 0.126046, 0.142715, 0.126762, 0.109894, 0.102035, 0.116802, 0.107059, 0.0922544, 0.108704, 0.101715, 0.100786, 0.0866739, 0.0941477, 0.0930732, 0.104544, 0.100266, 0.0838688, 0.0946095]
# 
# synaptic_weight_kick2_5 = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975]) * CONVERTION_CONSTANT
# probability_to_spike_kick2_5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0.07, 0.62, 0.62, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_spikes_kick2_5 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 7, 75, 69, 450, 457, 471, 450, 445, 432, 423, 423, 415, 412, 408, 407, 405, 410, 406, 403, 404, 406, 404, 402, 405, 401, 400]) / float(num_targets)
# mean_first_spike_delay_kick2_5 = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 6.98, 7.13429, 5.86968, 6.60548, 4.3594, 4.1892, 4.1994, 4.1466, 3.9238, 3.7664, 3.7868, 3.6012, 3.6062, 3.6186, 3.4958, 3.4888, 3.3864, 3.2954, 3.2484, 3.2626, 3.3674, 3.2178, 3.184, 3.107, 3.108, 3.0806, 3.0416]
# std_first_spike_delay_kick2_5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 2.30543, 1.76575, 2.48737, 0.263657, 0.187259, 0.273672, 0.414798, 0.152837, 0.130171, 0.146369, 0.127883, 0.111253, 0.10207, 0.116802, 0.106641, 0.0921957, 0.108687, 0.101272, 0.100691, 0.0866739, 0.0941477, 0.0930732, 0.104345, 0.100383, 0.0838688, 0.0946095]
# 
# synaptic_weight_kick3_0 = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975]) * CONVERTION_CONSTANT
# probability_to_spike_kick3_0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.02, 0.31, 0.32, 1, 1, 1, 0.99, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_spikes_kick3_0 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 39, 32, 392, 449, 409, 322, 443, 431, 423, 418, 413, 409, 408, 407, 403, 404, 403, 403, 404, 405, 403, 402, 404, 401, 400]) / float(num_targets)
# mean_first_spike_delay_kick3_0 = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 4.79, 5.26129, 6.30375, 4.4722, 4.2048, 4.3082, 4.44848, 3.926, 3.772, 3.7988, 3.6032, 3.6072, 3.6186, 3.4956, 3.4888, 3.3862, 3.2952, 3.2482, 3.2626, 3.3674, 3.2174, 3.184, 3.107, 3.108, 3.0806, 3.0416]
# std_first_spike_delay_kick3_0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.155563, 1.14083, 2.35807, 0.428159, 0.220726, 0.494504, 1.35103, 0.178535, 0.133817, 0.15468, 0.130785, 0.112941, 0.102976, 0.116604, 0.106641, 0.0921437, 0.10825, 0.101149, 0.100691, 0.0866739, 0.0938796, 0.0930732, 0.104345, 0.100383, 0.0838688, 0.0946095]
# 
# synaptic_weight_kick5_0 = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975]) * CONVERTION_CONSTANT
# probability_to_spike_kick5_0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.07, 0.03, 0.62, 0.97, 0.77, 0.61, 0.99, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_spikes_kick5_0 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 3, 139, 355, 147, 110, 396, 439, 437, 411, 406, 402, 401, 401, 400, 402, 400, 401, 400, 401, 400, 401, 400, 400, 400]) / float(num_targets)
# mean_first_spike_delay_kick5_0 = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 4.53429, 4.26667, 4.53161, 4.26289, 4.39792, 4.2682, 3.99051, 3.7914, 3.8774, 3.6086, 3.61, 3.6188, 3.4948, 3.4882, 3.3852, 3.2934, 3.2478, 3.2616, 3.3672, 3.2158, 3.184, 3.1064, 3.1078, 3.0806, 3.0416]
# std_first_spike_delay_kick5_0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.331806, 0.141892, 0.578722, 0.356113, 0.539147, 0.682558, 0.261773, 0.162537, 0.234442, 0.139211, 0.11743, 0.106014, 0.11671, 0.106481, 0.0914803, 0.107461, 0.101221, 0.10037, 0.0864809, 0.0937349, 0.0930732, 0.103938, 0.100097, 0.0838688, 0.0946095]
# 
# 
# f = plt.figure()
# 
# ax1 = f.add_subplot(311)
# ax1.plot(synaptic_weight_kick0, probability_to_spike_kick0, label='Ginh = 0 ms/cm^2')
# ax1.plot(synaptic_weight_kick0_5, probability_to_spike_kick0_5, label='Ginh = 0.5 ms/cm^2')
# ax1.plot(synaptic_weight_kick1_0, probability_to_spike_kick1_0, label='Ginh = 1.0 ms/cm^2')
# ax1.plot(synaptic_weight_kick1_5, probability_to_spike_kick1_5, label='Ginh = 1.5 ms/cm^2')
# ax1.plot(synaptic_weight_kick2_0, probability_to_spike_kick2_0, label='Ginh = 2.0 ms/cm^2')
# ax1.plot(synaptic_weight_kick2_5, probability_to_spike_kick2_5, label='Ginh = 2.5 ms/cm^2')
# ax1.plot(synaptic_weight_kick3_0, probability_to_spike_kick3_0, label='Ginh = 3.0 ms/cm^2')
# ax1.plot(synaptic_weight_kick5_0, probability_to_spike_kick5_0, label='Ginh = 5.0 ms/cm^2')
# ax1.set_ylabel('probability to spike')
# ax1.set_ylim([0, 1.1])
# ax1.legend(loc=4)
# 
# ax2 = f.add_subplot(312)
# ax2.plot(synaptic_weight_kick0, num_spikes_kick0, label='Ginh = 0 ms/cm^2')
# ax2.plot(synaptic_weight_kick0_5, num_spikes_kick0_5, label='Ginh = 0.5 ms/cm^2')
# ax2.plot(synaptic_weight_kick1_0, num_spikes_kick1_0, label='Ginh = 1.0 ms/cm^2')
# ax2.plot(synaptic_weight_kick1_5, num_spikes_kick1_5, label='Ginh = 1.5 ms/cm^2')
# ax2.plot(synaptic_weight_kick2_0, num_spikes_kick2_0, label='Ginh = 2.0 ms/cm^2')
# ax2.plot(synaptic_weight_kick2_5, num_spikes_kick2_5, label='Ginh = 2.5 ms/cm^2')
# ax2.plot(synaptic_weight_kick3_0, num_spikes_kick3_0, label='Ginh = 3.0 ms/cm^2')
# ax2.plot(synaptic_weight_kick5_0, num_spikes_kick5_0, label='Ginh = 5.0 ms/cm^2')
# ax2.set_ylabel('# of spikes')
# ax2.set_ylim([0, 7.1])
# ax2.legend(loc=2)
# 
# ax3 = f.add_subplot(313)
# ax3.plot(synaptic_weight_kick0, mean_first_spike_delay_kick0, label='Ginh = 0 ms/cm^2')
# ax3.plot(synaptic_weight_kick0_5, mean_first_spike_delay_kick0_5, label='Ginh = 0.5 ms/cm^2')
# ax3.plot(synaptic_weight_kick1_0, mean_first_spike_delay_kick1_0, label='Ginh = 1.0 ms/cm^2')
# ax3.plot(synaptic_weight_kick1_5, mean_first_spike_delay_kick1_5, label='Ginh = 1.5 ms/cm^2')
# ax3.plot(synaptic_weight_kick2_0, mean_first_spike_delay_kick2_0, label='Ginh = 2.0 ms/cm^2')
# ax3.plot(synaptic_weight_kick2_5, mean_first_spike_delay_kick2_5, label='Ginh = 2.5 ms/cm^2')
# ax3.plot(synaptic_weight_kick3_0, mean_first_spike_delay_kick3_0, label='Ginh = 3.0 ms/cm^2')
# ax3.plot(synaptic_weight_kick5_0, mean_first_spike_delay_kick5_0, label='Ginh = 5.0 ms/cm^2')
# ax3.set_ylim([0, 25])
# #ax3.plot(synaptic_weight_mature, first_spike_delay_mature, label='mature')
# ax3.set_ylabel('delay to first spike (ms)')
# ax3.set_xlabel('Single input weight (pS)')
# ax3.legend()
# 
#==============================================================================


#==============================================================================
# ##### Comparison between different noise level for immature model #####
# num_targets = 100
# 
# synaptic_weight_sd0_15 = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5]) * CONVERTION_CONSTANT
# probability_to_spike_sd0_15 = [0.02, 0.01, 0.02, 0.02, 0.12, 0.13, 0.19, 0.37, 0.59, 0.63, 0.77, 0.92, 0.95, 0.94, 1, 1, 1, 1, 1, 1]
# num_spikes_sd0_15 = np.array([2, 1, 2, 2, 14, 14, 26, 41, 79, 82, 114, 145, 152, 169, 210, 238, 265, 266, 313, 309]) / float(num_targets)
# mean_first_spike_delay_sd0_15 = [7.97, 81.64, 41.02, 18.27, 24.9267, 16.5708, 13.8579, 10.2589, 9.94508, 7.99238, 8.51169, 7.89217, 7.89432, 7.46596, 6.784, 6.4436, 5.9748, 5.2708, 4.3872, 4.2902]
# std_first_spike_delay_sd0_15 = [49.2571, -1, 49.95, 18.8515, 19.9306, 6.98254, 12.7389, 19.2945, 4.19392, 12.945, 12.3792, 5.82072, 1.61988, 1.72409, 1.53928, 1.17856, 1.03214, 6.24021, 9.90695, 9.58268]
# 
# synaptic_weight_sd0_10 = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5]) * CONVERTION_CONSTANT
# probability_to_spike_sd0_10 = [0, 0, 0, 0, 0.04, 0.02, 0.11, 0.2, 0.52, 0.65, 0.8, 0.93, 0.98, 0.97, 1, 1, 1, 1, 1, 1]
# num_spikes_sd0_10 = np.array([0, 0, 0, 0, 4, 2, 14, 22, 58, 73, 103, 133, 144, 165, 200, 230, 248, 255, 295, 305]) / float(num_targets)
# mean_first_spike_delay_sd0_10 = [-1, -1, -1, -1, 12.185, 11.15, 12.4782, 11.504, 10.5542, 10.8742, 9.83125, 8.69269, 7.97122, 7.47918, 6.705, 6.4792, 6.0194, 5.9, 5.6588, 5.4744]
# std_first_spike_delay_sd0_10 = [0, 0, 0, 0, 4.40247, 0.947523, 5.78491, 3.25964, 2.53278, 2.9247, 2.50339, 2.61691, 1.61323, 1.60589, 1.02756, 0.939572, 0.759444, 0.749373, 0.682567, 0.560729]
# 
# synaptic_weight_sd0_05 = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5]) * CONVERTION_CONSTANT
# probability_to_spike_sd0_05 = [0, 0, 0, 0, 0.02, 0.01, 0.03, 0.1, 0.52, 0.66, 0.84, 0.96, 0.98, 0.99, 1, 1, 1, 1, 1, 1]
# num_spikes_sd0_05 = np.array([0, 0, 0, 0, 2, 1, 3, 11, 56, 67, 94, 127, 143, 164, 199, 220, 242, 250, 287, 296]) / float(num_targets)
# mean_first_spike_delay_sd0_05 = [-1, -1, -1, -1, 16.15, 10.9, 10.3933, 12.142, 11.7985, 11.2261, 10.2064, 8.83312, 7.86388, 7.35535, 6.708, 6.4904, 6.0662, 5.8908, 5.6582, 5.4518]
# std_first_spike_delay_sd0_05 = [0, 0, 0, 0, 1.56978, -1, 0.91309, 2.87647, 3.28202, 2.57909, 2.66209, 2.00515, 1.13478, 1.16516, 0.759793, 0.77701, 0.625354, 0.576401, 0.561341, 0.465063]
# 
# 
# f = plt.figure()
# 
# ax1 = f.add_subplot(311)
# ax1.plot(synaptic_weight_sd0_05, probability_to_spike_sd0_05, label='$\sigma_d = 0.05$')
# ax1.plot(synaptic_weight_sd0_10, probability_to_spike_sd0_10, label='$\sigma_d = 0.10$')
# ax1.plot(synaptic_weight_sd0_15, probability_to_spike_sd0_15, label='$\sigma_d = 0.15$')
# ax1.set_ylabel('probability to spike')
# ax1.set_ylim([0, 1.1])
# ax1.legend(loc=4)
# 
# ax2 = f.add_subplot(312)
# ax2.plot(synaptic_weight_sd0_05, num_spikes_sd0_05, label='$\sigma_d = 0.05$')
# ax2.plot(synaptic_weight_sd0_10, num_spikes_sd0_10, label='$\sigma_d = 0.10$')
# ax2.plot(synaptic_weight_sd0_15, num_spikes_sd0_15, label='$\sigma_d = 0.15$')
# ax2.set_ylabel('# of spikes')
# ax2.set_ylim([0, 7.1])
# ax2.legend(loc=2)
# 
# ax3 = f.add_subplot(313)
# ax3.plot(synaptic_weight_sd0_05, mean_first_spike_delay_sd0_05, label='$\sigma_d = 0.05$')
# ax3.plot(synaptic_weight_sd0_10, mean_first_spike_delay_sd0_10, label='$\sigma_d = 0.10$')
# ax3.plot(synaptic_weight_sd0_15, mean_first_spike_delay_sd0_15, label='$\sigma_d = 0.15$')
# 
# ax3.set_ylim([0, 25])
# #ax3.plot(synaptic_weight_mature, first_spike_delay_mature, label='mature')
# ax3.set_ylabel('delay to first spike (ms)')
# ax3.set_xlabel('Single input weight (pS)')
# ax3.legend()
# 
# ##### Comparison between different inhibition strength for immature model E_GABA = -80 mV; sigma_s = sigma_d = 0.15#####
# num_targets = 100
# 
# synaptic_weight_kick0 = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5]) * CONVERTION_CONSTANT
# probability_to_spike_kick0 = [0.02, 0.01, 0.02, 0.02, 0.12, 0.13, 0.19, 0.37, 0.59, 0.63, 0.77, 0.92, 0.95, 0.94, 1, 1, 1, 1, 1, 1]
# num_spikes_kick0 = np.array([2, 1, 2, 2, 14, 14, 26, 41, 79, 82, 114, 145, 152, 169, 210, 238, 265, 266, 313, 309]) / float(num_targets)
# mean_first_spike_delay_kick0 = [7.97, 81.64, 41.02, 18.27, 24.9267, 16.5708, 13.8579, 10.2589, 9.94508, 7.99238, 8.51169, 7.89217, 7.89432, 7.46596, 6.784, 6.4436, 5.9748, 5.2708, 4.3872, 4.2902]
# std_first_spike_delay_kick0 = [49.2571, -1, 49.95, 18.8515, 19.9306, 6.98254, 12.7389, 19.2945, 4.19392, 12.945, 12.3792, 5.82072, 1.61988, 1.72409, 1.53928, 1.17856, 1.03214, 6.24021, 9.90695, 9.58268]
# 
# synaptic_weight_kick1_0 = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5]) * CONVERTION_CONSTANT
# probability_to_spike_kick1_0 = [0.01, 0.01, 0.01, 0.01, 0.03, 0.02, 0.08, 0.05, 0.18, 0.2, 0.36, 0.48, 0.58, 0.64, 0.88, 0.93, 0.94, 0.99, 1, 1]
# num_spikes_kick1_0 = np.array([1, 1, 1, 1, 3, 2, 10, 6, 21, 22, 42, 56, 67, 81, 124, 146, 171, 172, 221, 227]) / float(num_targets)
# mean_first_spike_delay_kick1_0 = [-26.86, 81.6, 76, 33.66, 52.12, 23.7, 21.3325, -4.692, 9.96556, 6.557, 9.64056, 11.0042, 11.54, 9.8125, 8.94136, 8.32946, 7.45638, 6.48303, 5.553, 5.3206]
# std_first_spike_delay_kick1_0 = [-1, -1, -1, -1, 23.4408, 17.3383, 18.1859, 51.6647, 5.57892, 23.2702, 18.1854, 9.57102, 3.89605, 3.07487, 3.16, 2.34188, 1.77348, 6.53378, 10.1676, 9.81256]
# 
# synaptic_weight_kick2_0 = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5]) * CONVERTION_CONSTANT
# probability_to_spike_kick2_0 = [0.01, 0.01, 0.01, 0.01, 0.02, 0.01, 0.02, 0.03, 0.04, 0.04, 0.1, 0.14, 0.1, 0.24, 0.41, 0.51, 0.7, 0.73, 0.82, 0.88]
# num_spikes_kick2_0 = np.array([1, 1, 1, 1, 2, 1, 3, 3, 5, 4, 10, 14, 10, 24, 47, 59, 85, 90, 125, 127]) / float(num_targets)
# mean_first_spike_delay_kick2_0 = [-26.86, 81.58, 75.74, 34.86, 63.7, 36.3, 36.79, -17.44, 15.4, -20.835, 6.944, 11.98, 19.33, 12.6975, 11.5702, 12.0063, 11.334, 10.1162, 7.67195, 7.51159]
# std_first_spike_delay_kick2_0 = [-1, -1, -1, -1, 17.4231, -1, 32.4562, 68.388, 19.9164, 46.1073, 35.5873, 18.4381, 13.0669, 4.35803, 4.44122, 4.43241, 5.021, 9.80106, 11.9218, 10.8334]
# 
# 
# 
# f = plt.figure()
# 
# ax1 = f.add_subplot(311)
# ax1.plot(synaptic_weight_kick0, probability_to_spike_kick0, label='Ginh = 0 ms/cm^2')
# ax1.plot(synaptic_weight_kick1_0, probability_to_spike_kick1_0, label='Ginh = 1.0 ms/cm^2')
# ax1.plot(synaptic_weight_kick2_0, probability_to_spike_kick2_0, label='Ginh = 2.0 ms/cm^2')
# ax1.set_ylabel('probability to spike')
# ax1.set_ylim([0, 1.1])
# ax1.legend(loc=4)
# 
# ax2 = f.add_subplot(312)
# ax2.plot(synaptic_weight_kick0, num_spikes_kick0, label='Ginh = 0 ms/cm^2')
# ax2.plot(synaptic_weight_kick1_0, num_spikes_kick1_0, label='Ginh = 1.0 ms/cm^2')
# ax2.plot(synaptic_weight_kick2_0, num_spikes_kick2_0, label='Ginh = 2.0 ms/cm^2')
# ax2.set_ylabel('# of spikes')
# ax2.set_ylim([0, 7.1])
# ax2.legend(loc=2)
# 
# ax3 = f.add_subplot(313)
# ax3.plot(synaptic_weight_kick0, mean_first_spike_delay_kick0, label='Ginh = 0 ms/cm^2')
# ax3.plot(synaptic_weight_kick1_0, mean_first_spike_delay_kick1_0, label='Ginh = 1.0 ms/cm^2')
# ax3.plot(synaptic_weight_kick2_0, mean_first_spike_delay_kick2_0, label='Ginh = 2.0 ms/cm^2')
# ax3.set_ylim([0, 25])
# #ax3.plot(synaptic_weight_mature, first_spike_delay_mature, label='mature')
# ax3.set_ylabel('delay to first spike (ms)')
# ax3.set_xlabel('Single input weight (pS)')
# ax3.legend()
#==============================================================================


###############################################################
### KLT current, same dendritic size and coupling constant, ###
### wider calcium burst and rest potential -55 mV           ###
###############################################################

#==============================================================================
# ################### Kick 0 Different neuron ages, no noise #########################
# 
# num_neurons = 10
# 
# weights_kick0_age0 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390]
# p_fire_soma_kick0_age0 = [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_soma_spikes_kick0_age0 = np.array([0, 0, 0, 0, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 30]) / np.float(num_neurons)
# mean_spike_delay_kick0_age0 = [-1, -1, -1, -1, 9.84, 7.98, 7.08, 6.46, 6, 5.64, 5.34, 5.12, 4.92, 4.74, 4.6, 4.46, 4.34, 4.24, 4.14, 4.04, 3.96, 3.88, 3.8, 3.74, 3.68, 3.62, 3.56, 3.52, 3.46, 3.42, 3.36, 3.32, 3.28, 3.24, 3.22, 3.18, 3.14, 3.12, 3.08, 3.06]
# std_spike_delay_kick0_age0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# 
# weights_kick0_age01 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390]
# p_fire_soma_kick0_age01 = [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_soma_spikes_kick0_age01 = np.array([0, 0, 0, 0, 10, 10, 10, 10, 10, 10, 10, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20]) / np.float(num_neurons)
# mean_spike_delay_kick0_age01 = [-1, -1, -1, -1, 11.22, 8.46, 7.46, 6.78, 6.28, 5.9, 5.58, 5.32, 5.1, 4.92, 4.78, 4.64, 4.5, 4.4, 4.28, 4.2, 4.1, 4.02, 3.94, 3.88, 3.82, 3.76, 3.7, 3.64, 3.58, 3.54, 3.48, 3.44, 3.4, 3.36, 3.32, 3.28, 3.26, 3.22, 3.18, 3.16]
# std_spike_delay_kick0_age01 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# 
# weights_kick0_age03 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390]
# p_fire_soma_kick0_age03 = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_soma_spikes_kick0_age03 = np.array([0, 0, 0, 0, 0, 10, 10, 10, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20]) / np.float(num_neurons)
# mean_spike_delay_kick0_age03 = [-1, -1, -1, -1, -1, 9.72, 8.28, 7.46, 6.88, 6.5, 6.14, 5.82, 5.56, 5.34, 5.16, 5, 4.86, 4.72, 4.6, 4.5, 4.4, 4.32, 4.24, 4.16, 4.08, 4.02, 3.96, 3.9, 3.84, 3.78, 3.74, 3.68, 3.64, 3.6, 3.54, 3.5, 3.46, 3.44, 3.4, 3.36]
# std_spike_delay_kick0_age03 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# 
# weights_kick0_age05 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390]
# p_fire_soma_kick0_age05 = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_soma_spikes_kick0_age05 = np.array([0, 0, 0, 0, 0, 10, 10, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20]) / np.float(num_neurons)
# mean_spike_delay_kick0_age05 = [-1, -1, -1, -1, -1, 12.56, 9.26, 8.18, 7.6, 7.08, 6.64, 6.3, 6, 5.76, 5.54, 5.34, 5.18, 5.04, 4.9, 4.78, 4.68, 4.58, 4.5, 4.42, 4.34, 4.26, 4.18, 4.12, 4.06, 4, 3.96, 3.9, 3.84, 3.8, 3.76, 3.72, 3.66, 3.62, 3.6, 3.56]
# std_spike_delay_kick0_age05 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# 
# weights_kick0_age100 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390]
# p_fire_soma_kick0_age100 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_soma_spikes_kick0_age100 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40]) / np.float(num_neurons)
# mean_spike_delay_kick0_age100 = [-1, -1, -1, -1, -1, -1, -1, -1, -1, 10.98, 9.86, 8.98, 8.38, 7.92, 7.54, 7.22, 6.92, 6.68, 6.44, 6.24, 6.04, 5.86, 5.7, 5.56, 5.44, 5.32, 5.2, 5.1, 5.02, 4.92, 4.84, 4.78, 4.7, 4.64, 4.58, 4.5, 4.46, 4.4, 4.34, 4.3]
# std_spike_delay_kick0_age100 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# 
# 
# 
# label1 = "age 0"
# label2 = "age 0.1"
# label3 = "age 0.3"
# label4 = "age 0.5"
# label5 = "age 100"
# 
# 
# f = plt.figure()
# 
# ax1 = f.add_subplot(411)
# ax1.plot(weights_kick0_age0, p_fire_soma_kick0_age0, label=label1)
# ax1.plot(weights_kick0_age01, p_fire_soma_kick0_age01, label=label2)
# ax1.plot(weights_kick0_age03, p_fire_soma_kick0_age03, label=label3)
# ax1.plot(weights_kick0_age05, p_fire_soma_kick0_age05, label=label4)
# ax1.plot(weights_kick0_age100, p_fire_soma_kick0_age100, label=label5)
# ax1.set_ylabel('p')
# 
# plt.legend(loc=9, bbox_to_anchor=(0.5, 1.3), ncol=5)
# 
# ax2 = f.add_subplot(412)
# ax2.plot(weights_kick0_age0, num_soma_spikes_kick0_age0, label=label1)
# ax2.plot(weights_kick0_age01, num_soma_spikes_kick0_age01, label=label2)
# ax2.plot(weights_kick0_age03, num_soma_spikes_kick0_age03, label=label3)
# ax2.plot(weights_kick0_age05, num_soma_spikes_kick0_age05, label=label4)
# ax2.plot(weights_kick0_age100, num_soma_spikes_kick0_age100, label=label5)
# ax2.set_ylabel('# spikes')
# 
# ax3 = f.add_subplot(413)
# ax3.plot(weights_kick0_age0, mean_spike_delay_kick0_age0, label=label1)
# ax3.plot(weights_kick0_age01, mean_spike_delay_kick0_age01, label=label2)
# ax3.plot(weights_kick0_age03, mean_spike_delay_kick0_age03, label=label3)
# ax3.plot(weights_kick0_age05, mean_spike_delay_kick0_age05, label=label4)
# ax3.plot(weights_kick0_age100, mean_spike_delay_kick0_age100, label=label5)
# ax3.set_ylabel('spike delay (ms)')
# 
# ax4 = f.add_subplot(414)
# ax4.plot(weights_kick0_age0, std_spike_delay_kick0_age0, label=label1)
# ax4.plot(weights_kick0_age01, std_spike_delay_kick0_age01, label=label2)
# ax4.plot(weights_kick0_age03, std_spike_delay_kick0_age03, label=label3)
# ax4.plot(weights_kick0_age05, std_spike_delay_kick0_age05, label=label4)
# ax4.plot(weights_kick0_age100, std_spike_delay_kick0_age100, label=label5)
# ax4.set_ylabel('std spike delay (ms)')
# ax4.set_xlabel('Single weight')
# 
#==============================================================================

################### Age100 Different inhibitory kicks, no noise #########################
#==============================================================================
# 
# num_neurons = 10
# 
# weights_kick0_age100 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390]
# p_fire_soma_kick0_age100 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_soma_spikes_kick0_age100 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40]) / np.float(num_neurons)
# mean_spike_delay_kick0_age100 = [-1, -1, -1, -1, -1, -1, -1, -1, -1, 10.98, 9.86, 8.98, 8.38, 7.92, 7.54, 7.22, 6.92, 6.68, 6.44, 6.24, 6.04, 5.86, 5.7, 5.56, 5.44, 5.32, 5.2, 5.1, 5.02, 4.92, 4.84, 4.78, 4.7, 4.64, 4.58, 4.5, 4.46, 4.4, 4.34, 4.3]
# std_spike_delay_kick0_age100 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# 
# weights_kick05_age100 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390]
# p_fire_soma_kick05_age100 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_soma_spikes_kick05_age100 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40]) / np.float(num_neurons)
# mean_spike_delay_kick05_age100 = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 10.76, 10.08, 9.28, 8.7, 8.22, 7.82, 7.46, 7.14, 6.82, 6.52, 6.24, 5.98, 5.74, 5.52, 5.32, 5.16, 5, 4.86, 4.74, 4.64, 4.54, 4.44]
# std_spike_delay_kick05_age100 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# 
# weights_kick1_age100 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390]
# p_fire_soma_kick1_age100 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_soma_spikes_kick1_age100 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 40, 10, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40]) / np.float(num_neurons)
# mean_spike_delay_kick1_age100 = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 11.64, 11.16, 9.08, 8.7, 7.98, 7.32, 6.68, 6.08, 5.6, 5.26, 5, 4.78]
# std_spike_delay_kick1_age100 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# 
# weights_kick15_age100 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490]
# p_fire_soma_kick15_age100 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_soma_spikes_kick15_age100 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 50, 40, 40, 40, 40, 40, 40, 40, 40, 40]) / np.float(num_neurons)
# mean_spike_delay_kick15_age100 = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 6.7, 5.48, 4.98, 4.86, 4.6, 4.4, 4.26, 4.16, 4.06, 4, 3.94, 3.9]
# std_spike_delay_kick15_age100 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# 
# 
# 
# label1 = "G_inh = 0.0"
# label2 = "G_inh = 0.5"
# label3 = "G_inh = 1.0"
# label4 = "G_inh = 1.5"
# 
# 
# f = plt.figure()
# 
# ax1 = f.add_subplot(411)
# ax1.plot(weights_kick0_age100, p_fire_soma_kick0_age100, label=label1)
# ax1.plot(weights_kick05_age100, p_fire_soma_kick05_age100, label=label2)
# ax1.plot(weights_kick1_age100, p_fire_soma_kick1_age100, label=label3)
# ax1.plot(weights_kick15_age100, p_fire_soma_kick15_age100, label=label4)
# ax1.set_ylabel('p')
# 
# plt.legend(loc=9, bbox_to_anchor=(0.5, 1.3), ncol=5)
# 
# ax2 = f.add_subplot(412)
# ax2.plot(weights_kick0_age100, num_soma_spikes_kick0_age100, label=label1)
# ax2.plot(weights_kick05_age100, num_soma_spikes_kick05_age100, label=label2)
# ax2.plot(weights_kick1_age100, num_soma_spikes_kick1_age100, label=label3)
# ax2.plot(weights_kick15_age100, num_soma_spikes_kick15_age100, label=label4)
# ax2.set_ylabel('# spikes')
# 
# ax3 = f.add_subplot(413)
# ax3.plot(weights_kick0_age100, mean_spike_delay_kick0_age100, label=label1)
# ax3.plot(weights_kick05_age100, mean_spike_delay_kick05_age100, label=label2)
# ax3.plot(weights_kick1_age100, mean_spike_delay_kick1_age100, label=label3)
# ax3.plot(weights_kick15_age100, mean_spike_delay_kick15_age100, label=label4)
# ax3.set_ylabel('spike delay (ms)')
# 
# ax4 = f.add_subplot(414)
# ax4.plot(weights_kick0_age100, std_spike_delay_kick0_age100, label=label1)
# ax4.plot(weights_kick05_age100, std_spike_delay_kick05_age100, label=label2)
# ax4.plot(weights_kick1_age100, std_spike_delay_kick1_age100, label=label3)
# ax4.plot(weights_kick15_age100, std_spike_delay_kick15_age100, label=label4)
# ax4.set_ylabel('std spike delay (ms)')
# ax4.set_xlabel('Single weight')
#==============================================================================

#==============================================================================
# ################### Age0 Different inhibitory kicks, no noise #########################
# 
# num_neurons = 10
# 
# weights_kick0_age0 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390]
# p_fire_soma_kick0_age0 = [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_soma_spikes_kick0_age0 = np.array([0, 0, 0, 0, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 30]) / np.float(num_neurons)
# mean_spike_delay_kick0_age0 = [-1, -1, -1, -1, 9.84, 7.98, 7.08, 6.46, 6, 5.64, 5.34, 5.12, 4.92, 4.74, 4.6, 4.46, 4.34, 4.24, 4.14, 4.04, 3.96, 3.88, 3.8, 3.74, 3.68, 3.62, 3.56, 3.52, 3.46, 3.42, 3.36, 3.32, 3.28, 3.24, 3.22, 3.18, 3.14, 3.12, 3.08, 3.06]
# std_spike_delay_kick0_age0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# 
# weights_kick05_age0 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490]
# p_fire_soma_kick05_age0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_soma_spikes_kick05_age0 = np.array([00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20]) / np.float(num_neurons)
# mean_spike_delay_kick05_age0 = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 6.24, 5.26, 4.84, 4.58, 4.38, 4.22, 4.1, 4, 3.9, 3.82, 3.74, 3.68, 3.62, 3.56, 3.5, 3.46, 3.42, 3.36, 3.32, 3.28, 3.24, 3.22, 3.18, 3.14, 3.12, 3.08, 3.06, 3.02, 3, 2.98, 2.94, 2.92, 2.9, 2.88, 2.86, 2.84, 2.82]
# std_spike_delay_kick05_age0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# 
# 
# weights_kick1_age0 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490]
# p_fire_soma_kick1_age0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_soma_spikes_kick1_age0 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 20, 20, 20]) / np.float(num_neurons)
# mean_spike_delay_kick1_age0 = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 4.72, 4.4, 4.18, 4.04, 3.92, 3.84, 3.76, 3.68, 3.62, 3.56, 3.5, 3.46, 3.42, 3.36, 3.32, 3.28, 3.24, 3.22, 3.18, 3.14, 3.12, 3.08, 3.06, 3.02, 3, 2.98, 2.94, 2.92, 2.9, 2.88, 2.86, 2.84, 2.82]
# std_spike_delay_kick1_age0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# 
# 
# weights_kick15_age0 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490]
# p_fire_soma_kick15_age0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# num_soma_spikes_kick15_age0 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]) / np.float(num_neurons)
# mean_spike_delay_kick15_age0 = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 4.82, 4.32, 4.1, 3.96, 3.84, 3.76, 3.68, 3.62, 3.56, 3.5, 3.46, 3.42, 3.36, 3.32, 3.28, 3.24, 3.22, 3.18, 3.14, 3.12, 3.08, 3.06, 3.02, 3, 2.98, 2.94, 2.92, 2.9, 2.88, 2.86, 2.84, 2.82]
# std_spike_delay_kick15_age0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# 
# 
# 
# label1 = "G_inh = 0.0"
# label2 = "G_inh = 0.5"
# label3 = "G_inh = 1.0"
# label4 = "G_inh = 1.5"
# 
# 
# f = plt.figure()
# 
# ax1 = f.add_subplot(411)
# ax1.plot(weights_kick0_age0, p_fire_soma_kick0_age0, label=label1)
# ax1.plot(weights_kick05_age0, p_fire_soma_kick05_age0, label=label2)
# ax1.plot(weights_kick1_age0, p_fire_soma_kick1_age0, label=label3)
# ax1.plot(weights_kick15_age0, p_fire_soma_kick15_age0, label=label4)
# ax1.set_ylabel('p')
# 
# plt.legend(loc=9, bbox_to_anchor=(0.5, 1.3), ncol=5)
# 
# ax2 = f.add_subplot(412)
# ax2.plot(weights_kick0_age0, num_soma_spikes_kick0_age0, label=label1)
# ax2.plot(weights_kick05_age0, num_soma_spikes_kick05_age0, label=label2)
# ax2.plot(weights_kick1_age0, num_soma_spikes_kick1_age0, label=label3)
# ax2.plot(weights_kick15_age0, num_soma_spikes_kick15_age0, label=label4)
# ax2.set_ylabel('# spikes')
# 
# ax3 = f.add_subplot(413)
# ax3.plot(weights_kick0_age0, mean_spike_delay_kick0_age0, label=label1)
# ax3.plot(weights_kick05_age0, mean_spike_delay_kick05_age0, label=label2)
# ax3.plot(weights_kick1_age0, mean_spike_delay_kick1_age0, label=label3)
# ax3.plot(weights_kick15_age0, mean_spike_delay_kick15_age0, label=label4)
# ax3.set_ylabel('spike delay (ms)')
# 
# ax4 = f.add_subplot(414)
# ax4.plot(weights_kick0_age0, std_spike_delay_kick0_age0, label=label1)
# ax4.plot(weights_kick05_age0, std_spike_delay_kick05_age0, label=label2)
# ax4.plot(weights_kick1_age0, std_spike_delay_kick1_age0, label=label3)
# ax4.plot(weights_kick15_age0, std_spike_delay_kick15_age0, label=label4)
# ax4.set_ylabel('std spike delay (ms)')
# ax4.set_xlabel('Single weight')
#==============================================================================

################### Kick = 0 Different age, noise sigma_s = 0.1 sigma_d = 0.2 #########################
num_neurons = 200

weights_kick0_age0 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]
p_fire_soma_kick0_age0 = [0.165, 0.145, 0.275, 0.565, 0.75, 0.91, 0.985, 0.975, 1, 1, 1, 1, 1, 1, 1]
num_soma_spikes_kick0_age0 = np.array([35, 29, 57, 123, 160, 206, 218, 226, 224, 231, 240, 264, 303, 329, 351]) / np.float(num_neurons)
mean_spike_delay_kick0_age0 = [61.1175, 41.421, 18.3746, 11.6283, 10.8445, 8.45598, 7.50235, 6.98155, 6.2569, 5.8366, 5.6851, 5.36281, 5.1509, 4.9106, 4.8432]
std_spike_delay_kick0_age0 = [25.3568, 29.7596, 22.4709, 12.1712, 11.1198, 3.41769, 1.81887, 3.1421, 1.17583, 1.13616, 1.05611, 0.900157, 0.803417, 0.715089, 0.852478]

weights_kick0_age01 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]
p_fire_soma_kick0_age01 = [0.06, 0.045, 0.15, 0.425, 0.665, 0.88, 0.98, 0.985, 0.995, 1, 1, 1, 1, 1, 1]
num_soma_spikes_kick0_age01 = np.array([13, 9, 30, 87, 134, 183, 203, 208, 225, 249, 298, 356, 403, 416, 429]) / np.float(num_neurons)
mean_spike_delay_kick0_age01 = [58.6633, 18.0667, 10.5048, 10.9163, 10.9394, 8.68846, 7.8499, 7.31959, 6.55045, 6.02465, 5.8637, 5.5229, 5.3291, 5.0984, 4.9609]
std_spike_delay_kick0_age01 = [27.2543, 15.0546, 3.36924, 8.56766, 10.1911, 1.83028, 1.80584, 3.16687, 1.17971, 0.99648, 0.934328, 0.851494, 0.764882, 0.673431, 0.827935]

weights_kick0_age03 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]
p_fire_soma_kick0_age03 = [0.005, 0.005, 0.06, 0.225, 0.48, 0.815, 0.95, 0.985, 0.995, 1, 1, 1, 1, 1, 1]
num_soma_spikes_kick0_age03 = np.array([1, 1, 12, 45, 96, 165, 222, 277, 361, 415, 427, 417, 420, 407, 403]) / np.float(num_neurons)
mean_spike_delay_kick0_age03 = [-1, 17.46, 10.565, 11.0509, 10.1215, 9.54626, 8.58063, 7.79178, 7.14342, 6.6255, 6.3813, 6.0158, 5.772, 5.5299, 5.3218]
std_spike_delay_kick0_age03 = [0, -1, 2.87785, 2.5328, 1.93418, 1.81011, 1.75983, 1.34743, 1.13395, 0.998463, 0.937642, 0.765892, 0.754809, 0.650331, 0.650159]

weights_kick0_age05 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]
p_fire_soma_kick0_age05 = [0.005, 0, 0.03, 0.075, 0.28, 0.69, 0.9, 0.955, 0.995, 1, 1, 1, 1, 1, 1]
num_soma_spikes_kick0_age05 = np.array([1, 0, 6, 15, 58, 150, 273, 383, 504, 569, 584, 583, 583, 579, 573]) / np.float(num_neurons)
mean_spike_delay_kick0_age05 = [-1, -1, 10.32, 11.4493, 10.685, 10.297, 9.35133, 8.44607, 7.78794, 7.1765, 6.8906, 6.4818, 6.1867, 5.9126, 5.6845]
std_spike_delay_kick0_age05 = [0, 0, 1.82761, 2.68739, 1.87911, 1.84457, 1.7635, 1.33139, 1.21048, 1.05991, 0.945859, 0.745692, 0.752447, 0.649107, 0.62312]

weights_kick0_age100 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]
p_fire_soma_kick0_age100 = [0, 0, 0, 0, 0, 0.005, 0.04, 0.155, 0.485, 0.87, 0.975, 1, 1, 1, 1]
num_soma_spikes_kick0_age100 = np.array([0, 0, 0, 0, 0, 1, 24, 96, 359, 689, 790, 810, 805, 803, 801]) / np.float(num_neurons)
mean_spike_delay_kick0_age100 = [-1, -1, -1, -1, -1, 11.76, 12.5575, 12.6658, 12.2596, 11.0883, 10.0031, 9.1218, 8.5331, 8.0472, 7.6382]
std_spike_delay_kick0_age100 = [0, 0, 0, 0, 0, -1, 1.60069, 2.74318, 2.6639, 1.7793, 1.4152, 0.857581, 0.758424, 0.618473, 0.519212]



label1 = "age 0"
label2 = "age 0.1"
label3 = "age 0.3"
label4 = "age 0.5"
label5 = "age 100"


f = plt.figure()

ax1 = f.add_subplot(411)
ax1.plot(weights_kick0_age0, p_fire_soma_kick0_age0, label=label1)
ax1.plot(weights_kick0_age01, p_fire_soma_kick0_age01, label=label2)
ax1.plot(weights_kick0_age03, p_fire_soma_kick0_age03, label=label3)
ax1.plot(weights_kick0_age05, p_fire_soma_kick0_age05, label=label4)
ax1.plot(weights_kick0_age100, p_fire_soma_kick0_age100, label=label5)
ax1.set_ylabel('p')

plt.legend(loc=9, bbox_to_anchor=(0.5, 1.3), ncol=5)

ax2 = f.add_subplot(412)
ax2.plot(weights_kick0_age0, num_soma_spikes_kick0_age0, label=label1)
ax2.plot(weights_kick0_age01, num_soma_spikes_kick0_age01, label=label2)
ax2.plot(weights_kick0_age03, num_soma_spikes_kick0_age03, label=label3)
ax2.plot(weights_kick0_age05, num_soma_spikes_kick0_age05, label=label4)
ax2.plot(weights_kick0_age100, num_soma_spikes_kick0_age100, label=label5)
ax2.set_ylabel('# spikes')

ax3 = f.add_subplot(413)
ax3.plot(weights_kick0_age0, mean_spike_delay_kick0_age0, label=label1)
ax3.plot(weights_kick0_age01, mean_spike_delay_kick0_age01, label=label2)
ax3.plot(weights_kick0_age03, mean_spike_delay_kick0_age03, label=label3)
ax3.plot(weights_kick0_age05, mean_spike_delay_kick0_age05, label=label4)
ax3.plot(weights_kick0_age100, mean_spike_delay_kick0_age100, label=label5)
ax3.set_ylabel('spike delay (ms)')

ax4 = f.add_subplot(414)
ax4.plot(weights_kick0_age0, std_spike_delay_kick0_age0, label=label1)
ax4.plot(weights_kick0_age01, std_spike_delay_kick0_age01, label=label2)
ax4.plot(weights_kick0_age03, std_spike_delay_kick0_age03, label=label3)
ax4.plot(weights_kick0_age05, std_spike_delay_kick0_age05, label=label4)
ax4.plot(weights_kick0_age100, std_spike_delay_kick0_age100, label=label5)
ax4.set_ylabel('std spike delay (ms)')
ax4.set_xlabel('Single weight')


plt.show()
