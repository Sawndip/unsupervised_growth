# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 17:55:39 2017

@author: jingroup

Script plots repsonse of HVC(RA) neuron to inhibitory kick
"""
import reading
import matplotlib.pyplot as plt

fileRA = "/home/eugene/forFigures/RA_inh_input.bin"
fileI = "/home/eugene/forFigures/I.bin"


(t, Vs, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _) = reading.read_hh2(fileRA)
(t, v, _, _, _, _, _, _, _, _, _) = reading.read_hhi(fileI)

f = plt.figure()

ax = f.add_subplot(111)
ax.plot(t, Vs, 'k')
ax.set_xlim([200, 330])

f = plt.figure()

ax = f.add_subplot(111)
ax.plot(t, v, 'k')
ax.set_xlim([200, 330])


plt.show()
