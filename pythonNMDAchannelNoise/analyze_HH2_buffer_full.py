#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 13:19:31 2017

@author: jingroup

Script plots full activity of HH2_buffer neuron
"""
import reading
import matplotlib.pyplot as plt
import numpy as np

#filename = "/home/eugene/Output/neuronTest/saturatedInhibitionResponse/Ginh_5.00.bin"
#filename = "/home/eugene/Output/neuronTest/kickResponse/Ginh_10.00_8.bin"
#filename = "/home/eugene/Output/neuronTest/response.bin"
filename = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/noImmatureOut1/RA/RA40.bin"
#filename = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays5/RA/RA5.bin"
#filename = "/home/eugene/Output/neuronTest/inhAndExcInputsResponse/RA.bin"
#filename = "/home/eugene/Output/neuronTest/modelStability/RA29.bin"
#filename = "/mnt/hodgkin_home/eugene/Output/tuneHVCRA/Ad1000/Rc55/RA26.bin"
#filename = "/mnt/hodgkin_home/eugene/Output/tuneHVCRA/Ad1000/Rc200/RA23.bin"

#filename = "/mnt/hodgkin_home/eugene/Output/neuronTest/inhAndExcInputsResponse/RA13.bin"

t, Vs, Vd, Gexc_d, Ginh_d, n, h, r, c, Ca = reading.read_hh2_buffer_full(filename)

#print t
#print Vs

print Vs[-1]
print Vd[-1]
print n[-1]
print h[-1]
print r[-1]
print c[-1]
print Ca[-1]

print "Average Vs = ", np.mean(Vs)
print "Std Vs = ", np.std(Vs)

#print t
#print Vs
tmin = 30
tmax = 200

# membrane potentials
f = plt.figure()

ax1 = f.add_subplot(211)
ax1.plot(t, Vs)
ax1.set_ylabel("Vs (mV)")
ax1.set_xlim([tmin, tmax])

ax2 = f.add_subplot(212)
ax2.plot(t, Vd)
ax2.set_ylabel("Vd (mV)")
ax2.set_xlabel("Time (ms)")
ax2.set_xlim([tmin, tmax])

# conductances
f = plt.figure()

ax1 = f.add_subplot(211)
ax1.plot(t, Gexc_d)
ax1.set_ylabel("Gexc_d (mS/cm^2)")
ax1.set_xlim([tmin, tmax])

ax2 = f.add_subplot(212)
ax2.plot(t, Ginh_d)
ax2.set_ylabel("Ginh_d (mS/cm^2)")
ax2.set_xlabel("Time (ms)")
ax2.set_xlim([tmin, tmax])

# dendritic compartment variables
f = plt.figure()

ax1 = f.add_subplot(311)
ax1.plot(t, r)
ax1.set_ylabel("r")
ax1.set_xlim([tmin, tmax])

ax2 = f.add_subplot(312)
ax2.plot(t, Ca)
ax2.set_ylabel("Ca")
ax2.set_xlim([tmin, tmax])

ax3 = f.add_subplot(313)
ax3.plot(t, c)
ax3.set_ylabel("c")

ax3.set_xlabel("Time (ms)")
ax3.set_xlim([tmin, tmax])

# somatic compartment variables
f = plt.figure()

ax1 = f.add_subplot(211)
ax1.plot(t, n)
ax1.set_ylabel("n")
ax1.set_xlim([tmin, tmax])

ax2 = f.add_subplot(212)
ax2.plot(t, h)
ax2.set_ylabel("h")
ax2.set_xlabel("Time (ms)")
ax2.set_xlim([tmin, tmax])

plt.show()