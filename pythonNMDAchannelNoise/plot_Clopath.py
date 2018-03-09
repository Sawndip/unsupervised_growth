# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 14:53:01 2018

@author: jingroup

Sctipt plots results of clopath STDP model simulation
"""
import reading
import matplotlib.pyplot as plt


filename = "/home/eugene/Output/clopath/preBurstPostSpike/10.bin"
#filename = "/home/eugene/Output/clopath/frequency/preSpikePostSpike/10.bin"

time, vd, u_minus, u_plus, x, w = reading.read_clopath_test(filename)

THETA_MINUS = -75
THETA_PLUS = -30


f = plt.figure()

ax1 = f.add_subplot(311)
ax1.plot(time, vd, label='$v_{d}$', color='k')
ax1.plot(time, u_minus, label='$u_{-}$', color='m')
ax1.plot(time, u_plus, label='$u_{+}$', color='b')
ax1.axhline(y=THETA_MINUS, linestyle='dashed', color='m')
ax1.axhline(y=THETA_PLUS, linestyle='dashed', color='b')

ax1.set_ylabel('V (mV)')
ax1.legend()

ax2 = f.add_subplot(312)
ax2.plot(time, x, color='k')
ax2.set_ylabel('x')

ax3 = f.add_subplot(313)
ax3.plot(time, w, color='k')
ax3.set_ylabel('w')

plt.show()



