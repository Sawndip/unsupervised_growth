# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 13:19:39 2018

@author: jingroup

Script compares responses of different HVC-RA neurons
"""

import reading
import matplotlib.pyplot as plt
import numpy as np

#filename = "/home/eugene/Output/neuronTest/saturatedInhibitionResponse/Ginh_5.00.bin"
#filename = "/home/eugene/Output/neuronTest/kickResponse/Ginh_10.00_8.bin"
#filename = "/home/eugene/Output/neuronTest/response.bin"
#filename = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/test1/RA/RA231.bin"
#filename = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays5/RA/RA5.bin"
#filename = "/home/eugene/Output/neuronTest/inhAndExcInputsResponse/RA.bin"
#filename = "/home/eugene/Output/neuronTest/modelStability/RA29.bin"
filename1 = "/mnt/hodgkin_home/eugene/Output/tuneHVCRA/Ad1000/Rc200/RA30.bin"
filename2 = "/mnt/hodgkin_home/eugene/Output/tuneHVCRA/Ad1000/Rc200/RA25.bin"
#filename2 = "/mnt/hodgkin_home/eugene/Output/tuneHVCRA/passiveDendrite/NaKchange/Na80K8/RA23.bin"

#filename = "/mnt/hodgkin_home/eugene/Output/neuronTest/inhAndExcInputsResponse/RA13.bin"

t1, Vs1, Vd1, Gexc_d1, Ginh_d1, n1, h1, r1, c1, Ca1 = reading.read_hh2_buffer_full(filename1)
t2, Vs2, Vd2, Gexc_d2, Ginh_d2, n2, h2, r2, c2, Ca2 = reading.read_hh2_buffer_full(filename2)

#print t
#print Vs

label1 = 'Rc200 30 Ge1.05'
label2 = 'Rc200 25 Ge1.5'


tmin = 40
tmax = 100

# membrane potentials
f = plt.figure()

ax1 = f.add_subplot(211)
ax1.plot(t1, Vs1, label=label1)
ax1.plot(t2, Vs2, label=label2)
ax1.legend(loc=1)

ax1.set_ylabel("Vs (mV)")
ax1.set_xlim([tmin, tmax])

ax2 = f.add_subplot(212)
ax2.plot(t1, Vd1, label=label1)
ax2.plot(t2, Vd2, label=label2)
ax2.legend(loc=1)

ax2.set_ylabel("Vd (mV)")
ax2.set_xlabel("Time (ms)")
ax2.set_xlim([tmin, tmax])

# conductances
f = plt.figure()

ax1 = f.add_subplot(211)
ax1.plot(t1, Gexc_d1, label=label1)
ax1.plot(t2, Gexc_d2, label=label2)
ax1.legend()
ax1.set_ylabel("Gexc_d (mS/cm^2)")
ax1.set_xlim([tmin, tmax])

ax2 = f.add_subplot(212)
ax2.plot(t1, Ginh_d1, label=label1)
ax2.plot(t2, Ginh_d2, label=label2)
ax2.legend()
ax2.set_ylabel("Ginh_d (mS/cm^2)")
ax2.set_xlabel("Time (ms)")
ax2.set_xlim([tmin, tmax])

# dendritic compartment variables
f = plt.figure()

ax1 = f.add_subplot(311)
ax1.plot(t1, r1, label=label1)
ax1.plot(t2, r2, label=label2)
ax1.legend()
ax1.set_ylabel("r")
ax1.set_xlim([tmin, tmax])

ax2 = f.add_subplot(312)
ax2.plot(t1, Ca1, label=label1)
ax2.plot(t2, Ca2, label=label2)
ax2.legend()
ax2.set_ylabel("Ca")
ax2.set_xlim([tmin, tmax])

ax3 = f.add_subplot(313)
ax3.plot(t1, c1, label=label1)
ax3.plot(t2, c2, label=label2)
ax3.legend()
ax3.set_ylabel("c")

ax3.set_xlabel("Time (ms)")
ax3.set_xlim([tmin, tmax])

# somatic compartment variables
f = plt.figure()

ax1 = f.add_subplot(211)
ax1.plot(t1, n1, label=label1)
ax1.plot(t2, n2, label=label2)
ax1.legend(loc=1)
ax1.set_ylabel("n")
ax1.set_xlim([tmin, tmax])

ax2 = f.add_subplot(212)
ax2.plot(t1, h1, label=label1)
ax2.plot(t2, h2, label=label2)
ax2.legend(loc=4)
ax2.set_ylabel("h")
ax2.set_xlabel("Time (ms)")
ax2.set_xlim([tmin, tmax])

plt.show()

