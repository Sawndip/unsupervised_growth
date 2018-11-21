#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 15:17:10 2017

@author: jingroup

script plots results of f-I tuning curve for HVC-RA neuron
"""
import reading
import matplotlib.pyplot as plt
import os
import numpy as np

dirname = "/mnt/hodgkin/eugene/Output/tuneHVCRA/finalModel/KLTH3.5Rest55mVLargeIrange/age0/"


file_soma = os.path.join(dirname, "fI_soma.bin")
file_dend = os.path.join(dirname, "fI_dend.bin")


ampl_soma, num_spikes_soma, _, = reading.read_fI_HVCRA(file_soma)
ampl_dend, num_spikes_dend, _, = reading.read_fI_HVCRA(file_dend)

plt.figure()
plt.plot(ampl_soma, num_spikes_soma, '-bo')
plt.xlabel('I (nA)')
plt.ylabel('# of spikes')
plt.title('Injection to soma')
plt.legend()

plt.figure()
plt.plot(ampl_dend, num_spikes_dend, '-bo')
plt.xlabel('I (nA)')
plt.ylabel('# of spikes')
plt.title('Injection to dend')
plt.legend()

plt.show()