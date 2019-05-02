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
#filename = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/maturationTransition1/RA/RA570.bin"
#filename = "/home/eugene/Output/neuronTest/noise/noise_immature_4.bin"

#filename = "/mnt/hodgkin/eugene/Output/tuneHVCRA/finalModel/KLTH3.5Rest55mVExamples/age0/RA99.bin"
filename = "/mnt/hodgkin/eugene/Output/tuneHVCRA/finalModel/KLTH3.5Rest55mVExamples/age1.0/RA45.bin"

#filename = "/mnt/hodgkin/eugene/Output/tuneHVCRA/finalModel/KLTH/age100/RA7.bin"


#filename = "/mnt/hodgkin/eugene/Output/tuneHVCRA/finalModel/noKLTHtest/age100/RA40.bin"


#filename = "/mnt/hodgkin/eugene/Output/tuneHVCRA/finalModel/noKLTHCalciumSpike20msLeak0.75/age100/RA10.bin"


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

t = t - 52.7 + 1.4 - 1

tmin = -5
tmax = 30

#tmin = 48
#tmax = tmin + 35

#tmin = 0
#tmax = 500

f = plt.figure()

ax1 = f.add_subplot(111)
ax1.plot(t, Vs)
ax1.set_ylabel("Vs (mV)")
ax1.set_xlim([tmin, tmax])
ax1.set_xlabel('Time (ms)')

#==============================================================================
# # membrane potentials
# f = plt.figure()
# 
# ax1 = f.add_subplot(211)
# ax1.plot(t, Vs)
# ax1.set_ylabel("Vs (mV)")
# ax1.set_xlim([tmin, tmax])
# 
# ax2 = f.add_subplot(212)
# ax2.plot(t, Vd)
# ax2.set_ylabel("Vd (mV)")
# ax2.set_xlabel("Time (ms)")
# ax2.set_xlim([tmin, tmax])
# 
# # conductances
# f = plt.figure()
# 
# ax1 = f.add_subplot(211)
# ax1.plot(t, Gexc_d)
# ax1.set_ylabel("Gexc_d (mS/cm^2)")
# ax1.set_xlim([tmin, tmax])
# 
# ax2 = f.add_subplot(212)
# ax2.plot(t, Ginh_d)
# ax2.set_ylabel("Ginh_d (mS/cm^2)")
# ax2.set_xlabel("Time (ms)")
# ax2.set_xlim([tmin, tmax])
# 
# # dendritic compartment variables
# f = plt.figure()
# 
# ax1 = f.add_subplot(311)
# ax1.plot(t, r)
# ax1.set_ylabel("r")
# ax1.set_xlim([tmin, tmax])
# 
# ax2 = f.add_subplot(312)
# ax2.plot(t, Ca)
# ax2.set_ylabel("Ca")
# ax2.set_xlim([tmin, tmax])
# 
# ax3 = f.add_subplot(313)
# ax3.plot(t, c)
# ax3.set_ylabel("c")
# 
# ax3.set_xlabel("Time (ms)")
# ax3.set_xlim([tmin, tmax])
# 
# # somatic compartment variables
# f = plt.figure()
# 
# ax1 = f.add_subplot(211)
# ax1.plot(t, n)
# ax1.set_ylabel("n")
# ax1.set_xlim([tmin, tmax])
# 
# ax2 = f.add_subplot(212)
# ax2.plot(t, h)
# ax2.set_ylabel("h")
# ax2.set_xlabel("Time (ms)")
# ax2.set_xlim([tmin, tmax])
#==============================================================================

plt.show()