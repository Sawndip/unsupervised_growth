# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 02:08:33 2016

@author: Eugene
"""

#	This program reads binary data file that contains membrane potential points and
#	and corresponding time points. First number in the is the number of points
#	in datafile. All numbers are stored as doubles, which means each number needs 8 bytes
#	to be stored.
import os
import struct
import matplotlib.pyplot as plt
#from matplotlib.ticker import FuncFormatter
#from mpldatacursor import datacursor
import numpy as np
import reading

simFileAbs = "/home/eugene/Output/RAneurons/RA22.bin"	#	datafile name
#simFileAbs = "/home/eugene/Output/RA.bin"	#	datafile name

#testFileRel = "1.bin"
SIZE_OF_DOUBLE = 8
SIZE_OF_INT = 4

def print_steady_params():
    print "Steady state parameters:"
    print "Vs = ",Vs[-1]
    print "n = ",n[-1]
    print "h = ",h[-1]
    print "Vd = ",Vd[-1]
    print "r = ",r[-1]
    print "c = ",c[-1]
    print "Ca = ",Ca[-1]

    

#	open file and read all its content to fileContent variable


#script_dir = os.path.dirname(__file__)

#simFileAbs = os.path.join(script_dir, simFileRel)
#testFileAbs = os.path.join(script_dir, testFileRel)
#print repr(fileContent)


#print dataPointsNumber[0]
#print len(fileContent[8:])


#print voltage
#print time

#print len(voltage)
#print len(time)

(t, Vs, Is, n, h, Vd, Id, r, c, Ca, Gexc_d, Ginh_d, Gexc_s, Ginh_s, Ei, flag, Nsoma, Ndend) = reading.read_hh2(simFileAbs)

print Nsoma
#print "Vs", Vs

#print "s_soma", s_s
#print "s_dend", s_d

print "Ei = ", Ei[100]
print_steady_params()


#print(len(Vd))
#print(n[0])
#print(n[11999])
#print(h[0])
#rint(h[11999])
#print(r[0])
#print(r[11999])
#print(c[0])
#print(c[11999])
#rint(Ca[0])
#rint(Ca[11999])
for i in range(Is.size):
    Is[i] /= 1000
    Id[i] /= 1000

#print Is

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(511)
ax1.plot(t, Vs, 'r', linewidth = 2.0)
ax1.set_ylabel("Vs (mV)")
#ax1.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.2f')%(x)))
plt.title("Hodgkin-Huxley 2-compartment Neuron Model")

ax2 = fig1.add_subplot(512)
ax2.plot(t, Is, 'b', linewidth = 2.0)
ax2.set_ylabel("Is (nA)")

ax3 = fig1.add_subplot(513)
ax3.plot(t, Vd, 'r', linewidth = 2.0)
#ax3.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.2f')%(x)))
ax3.set_ylabel("Vd (mV)")

ax4 = fig1.add_subplot(514)
ax4.plot(t, Id, 'b', linewidth = 2.0)
ax4.set_ylabel("Id, ext (nA)")
ax4.set_xlabel("t (ms)")

ax5 = fig1.add_subplot(515)
ax5.plot(t, flag, 'm', linewidth = 2.0)
ax5.set_ylabel("flag")
ax5.set_xlabel("t (ms)")


plt.grid(True)

fig2 = plt.figure(2)
ax1 = fig2.add_subplot(511)
ax1.plot(t, n, 'r', linewidth = 2.0)
ax1.set_ylabel("n")
plt.title("Gating variables"),

ax2 = fig2.add_subplot(512)
ax2.plot(t, h, 'b', linewidth = 2.0)
ax2.set_ylabel("h")

ax3 = fig2.add_subplot(513)
ax3.plot(t, r, 'r', linewidth = 2.0)
ax3.set_ylabel("r")

ax4 = fig2.add_subplot(514)
ax4.plot(t, c, 'b', linewidth = 2.0)
ax4.set_ylabel("c")

ax5 = fig2.add_subplot(515)
ax5.plot(t, Ca, 'r', linewidth = 2.0)
ax5.set_ylabel("[Ca]")

ax5.set_xlabel("t (ms)")

plt.grid(True)


fig3 = plt.figure(3)
ax1 = fig3.add_subplot(411)
ax1.plot(t, Gexc_d, 'r', linewidth = 2.0)
ax1.set_ylabel("$G_{exc}^{d}$")
plt.title('Synaptic conductances. s - refers to soma, d - to dendrite compartments.')


ax3 = fig3.add_subplot(412)
ax3.plot(t, Ginh_d, 'm', linewidth = 2.0)
ax3.set_ylabel("$G_{inh}^{d}$")

ax4 = fig3.add_subplot(413)
ax4.plot(t, Gexc_s, 'y', linewidth = 2.0)
ax4.set_ylabel("$G_{exc}^{s}$")

ax5 = fig3.add_subplot(414)
ax5.plot(t, Ginh_s, 'g', linewidth = 2.0)
ax5.set_ylabel("$G_{inh}^{s}$")

ax4.set_xlabel("t (ms)")

plt.grid(True)

#datacursor(formatter="x: {x:.3f}\ny: {y:.3f}".format)
plt.show()


#print vSimulatedData[60000]
#print nSim[60000]
#print mSim[60000]
#print hSim[60000]
