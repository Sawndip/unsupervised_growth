# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 02:08:33 2016

@author: Eugene

Script plots HVC(RA) membrane potential and its inhibitory input
"""

import matplotlib.pyplot as plt
import reading

filename = "/home/eugene/hodgkinData/gabaMaturation300117/membraneTraces/RA/RA268_trial1.bin"	#	datafile name

(t, Vs, _, _, _, Vd, _, _, _, _, Gexc_d, Ginh_d, _, _, _, _, _, _) = reading.read_hh2(filename)


fig1 = plt.figure(1)
ax1 = fig1.add_subplot(411)
ax1.plot(t, Vs, 'r', linewidth = 2.0)
ax1.set_ylabel("Vs (mV)")
plt.title("HVC(RA) neuron")

ax2 = fig1.add_subplot(412)
ax2.plot(t, Vd, 'r', linewidth = 2.0)
ax2.set_ylabel("Vd (mV)")

ax3 = fig1.add_subplot(413)
ax3.plot(t, Gexc_d, 'r', linewidth = 2.0)
ax3.set_ylabel("$G_{exc}^{d}$")
ax3.set_xlabel("t (ms)")

ax4 = fig1.add_subplot(414)
ax4.plot(t, Ginh_d, 'r', linewidth = 2.0)
ax4.set_ylabel("$G_{inh}^{d}$")
ax4.set_xlabel("t (ms)")

plt.grid(True)

plt.show()