# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 13:57:03 2018

@author: jingroup

Script plots channel properties of HVC_RA neuron
"""
import matplotlib.pyplot as plt
import numpy as np

v = np.linspace(-80, +80, 1000)

print v
nInf = 1.0 / (1.0 + np.exp(-(v+35.0) / 10.0))
hInf = 1.0 / (1.0 + np.exp((v + 45.0) / 7.0))
mInf = 1.0 / (1.0 + np.exp(-(v + 30.0) / 9.5))
rInf = 1.0 / (1.0 + np.exp(-(v + 5.0) / 10.0))
cInf = 1.0 / (1.0 + np.exp(-(v - 10.0) / 7.0))

tauN = 0.1 + 0.5 / (1.0 + np.exp((v + 27.0) / 15.0))
tauH = 0.1 + 0.75 / (1.0 + np.exp((v + 40.5) / 6.0))

print nInf
print hInf
print mInf

print tauN
print tauH

f = plt.figure()

#plt.title("Na channel")

ax1 = f.add_subplot(311)
ax1.set_title('Activation')
ax1.plot(v, mInf)
ax1.set_ylabel("$m_{\infty}$")

ax2 = f.add_subplot(312)
ax2.set_title('Inactivation')
ax2.plot(v, hInf)
ax2.set_ylabel("$h_{\infty}$")

ax3 = f.add_subplot(313)
ax3.set_title('Inactivation time constant')
ax3.plot(v, tauH)
ax3.set_ylabel(r'$\tau_{h}$ (ms)')
ax3.set_xlabel('V (mV)')


f = plt.figure()

ax1 = f.add_subplot(211)
ax1.plot(v, nInf)
ax1.set_ylabel("$n_{\infty}$")

ax2 = f.add_subplot(212)
ax2.plot(v, tauN)
ax2.set_ylabel(r'$\tau_{n}$ (ms)')
ax2.set_xlabel('V (mV)')


# Calcium variables
f = plt.figure()

ax1 = f.add_subplot(211)
ax1.plot(v, rInf)
ax1.set_ylabel("$r_{\infty}$")

ax2 = f.add_subplot(212)
ax2.plot(v, cInf)
ax2.set_ylabel(r'$c_{\infty}$')
ax2.set_xlabel('V (mV)')

plt.show()

#static double nInf(double v){return 1 / (1 + exp(-(v + 35) / 10));} 
#	static double tauN(double v){return 0.1 + 0.5 / (1 + exp((v + 27) / 15));} 
#	static double hInf(double v){return 1 / (1 + exp((v + 45) / 7));} 
#	static double tauH(double v){return 0.1 + 0.75 / (1 + exp((v + 40.5) / 6));} 
#	static double mInf(double v){return 1 / (1 + exp(-(v + 30) / 9.5));} 

