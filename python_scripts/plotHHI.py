#	This program reads binary data file that contains membrane potential points and
#	and corresponding time points. First number in the is the number of points
#	in datafile. All numbers are stored as doubles, which means each number needs 8 bytes
#	to be stored.
import os
import struct
import matplotlib.pyplot as plt
import numpy as np
import reading

FileRel = "/mnt/hodgkin_home/eugene/Output/networks/networkTest/Ineurons/I1_trial_52_.bin"	#	datafile name
#testFileRel = "1.bin"
SIZE_OF_DOUBLE = 8
SIZE_OF_INT = 4
#	open file and read all its content to fileContent variable

script_dir = os.path.dirname(__file__)

FileAbs = os.path.join(script_dir, FileRel)

(t, v, I, n, m, h, w, Ge, Gi, flag, Nspikes) = reading.read_hhi(FileAbs)

#print I.size

def print_time_state(ind):
    print "Variables at time:", ind
    print "V = ",v[ind]
    print "n = ",n[ind]
    print "h = ",h[ind]
    print "m = ",m[ind]
    print "Ge = ",Ge[ind]
    print "Gi = ",Gi[ind]
    

for i in range(I.size):
   I[i] *= 0.06


print I

print v
#print I
#print I[100]
print Nspikes
print Nspikes / (t[-1] / 1000)

print_time_state(-1)

fig1 = plt.figure(1)

ax1 = fig1.add_subplot(311)
ax1.plot(t, v, 'r', linewidth = 2.0)
ax1.set_ylabel("v (mV)")
plt.title("Hodgkin-Huxley Neuron Model")

ax2 = fig1.add_subplot(312)
ax2.plot(t, I, 'b', linewidth = 2.0)
ax2.set_ylabel("Iext (nA/mm^2)")

ax3 = fig1.add_subplot(313)
ax3.plot(t, flag, 'm', linewidth = 2.0)
ax3.set_ylabel("flag")
ax3.set_xlabel("t (ms)")

plt.grid(True)

fig2 = plt.figure()

ax1 = fig2.add_subplot(411)
ax1.plot(t, n, 'r', linewidth = 2.0)
ax1.set_ylabel('n')

ax2 = fig2.add_subplot(412)
ax2.plot(t, m, 'm', linewidth = 2.0)
ax2.set_ylabel('m')

ax3 = fig2.add_subplot(413)
ax3.plot(t, h, 'b', linewidth = 2.0)
ax3.set_ylabel('h')

ax4 = fig2.add_subplot(414)
ax4.plot(t, w, 'g', linewidth = 2.0)
ax4.set_ylabel('w')
ax4.set_xlabel('t (ms)')

plt.grid(True)

fig3 = plt.figure()

ax1 = fig3.add_subplot(211)
ax1.plot(t, Ge, 'r', linewidth = 2.0)
ax1.set_ylabel('Gexc')

ax2 = fig3.add_subplot(212)
ax2.plot(t, Gi, 'b', linewidth = 2.0)
ax2.set_ylabel('Ginh')
ax2.set_xlabel('t (ms)')

plt.grid(True)

plt.show()
