# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 17:43:01 2017

@author: jingroup

Script plots membrane traces for same HVC(RA) neuron
"""
import reading
import matplotlib.pyplot as plt

file1 = "/home/eugene/hodgkinData/gabaMaturation010217/membraneTraces/RA207_trial1.bin"
file2 = "/home/eugene/hodgkinData/gabaMaturation010217/membraneTraces/RA207_trial2.bin"
file3 = "/home/eugene/hodgkinData/gabaMaturation010217/membraneTraces/RA207_trial3.bin"


(t_1, Vs_1, _, _, _, Vd_1, _, _, _, _, Gexc_d_1, Ginh_d_1, _, _, _, _, _, _) = reading.read_hh2(file1)
(t_2, Vs_2, _, _, _, Vd_2, _, _, _, _, Gexc_d_2, Ginh_d_2, _, _, _, _, _, _) = reading.read_hh2(file2)
(t_3, Vs_3, _, _, _, Vd_3, _, _, _, _, Gexc_d_3, Ginh_d_3, _, _, _, _, _, _) = reading.read_hh2(file3)


f = plt.figure()
ax1 = f.add_subplot(211)

ax1.plot(t_1, Vs_1, 'r')
ax1.plot(t_2, Vs_2, 'b')
ax1.plot(t_3, Vs_3, 'm')
ax1.set_ylabel("Vs (mV)")

ax2 = f.add_subplot(212)

ax2.plot(t_1, Gexc_d_1, 'r')
ax2.plot(t_2, Gexc_d_2, 'b')
ax2.plot(t_3, Gexc_d_3, 'm')
ax2.set_ylabel("Gexc_d")
ax2.set_xlabel("t (ms)")

plt.show()