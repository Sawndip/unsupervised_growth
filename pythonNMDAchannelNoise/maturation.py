# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 13:58:48 2017

@author: jingroup

Script analyzes file with maturation info
"""

import reading
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

filename = "/home/eugene/Output/networks/gabaMaturation100217/mature.bin"

trial_number, gaba_potential, firing_rate, remodeled, mature = reading.read_maturation_info(filename)

print firing_rate
#colors = [str(x) for x in firing_rate]
plt.figure()

colors = cm.coolwarm(np.array(firing_rate))
plot = plt.scatter(firing_rate, firing_rate, c = firing_rate, cmap = 'coolwarm')
plt.clf()
plt.colorbar(plot)
plt.bar(range(len(gaba_potential)), gaba_potential, edgecolor='none', color=colors)
plt.ylim([min(gaba_potential), max(gaba_potential) + 1])

plt.gca().invert_yaxis()
plt.xlabel("neuron id")
plt.ylabel("GABA potential")

plt.show()