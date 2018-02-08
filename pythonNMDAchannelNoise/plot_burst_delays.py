# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 13:21:35 2018

@author: jingroup

Script plots burst delay responses to inhibitory kicks
"""
import matplotlib.pyplot as plt

kicks = [0, 0.75, 1.5, 2.25, 3.0, 3.75, 4.5, 5.25, 6.0, 6.75, 7.5, 8.25, 9.0, 9.75, 10.5, 11.25, 12.0, 12.75, 13.5, 14.25]

probability = [0.001, 0.013, 0.097, 0.121, 0.154, 0.143, 0.133, 0.151, 0.134, 0.141, 0.117, 0.126, 0.130,\
                0.125, 0.124, 0.132, 0.129, 0.102, 0.141, 0.109]
                
burst_delays = [39.48, 21.04, 19.61, 19.73, 20.36, 20.91, 21.46, 22.76, 22.91, 23.85, 23.29, 24.38, 24.13,\
                25.95, 25.21, 25.72, 25.73, 26.63, 26.83, 27.17]



plt.figure()

plt.plot(kicks, probability)
plt.xlabel('Gi (mS/cm^2)')
plt.ylabel('probability to burst')

plt.figure()
plt.plot(kicks, burst_delays)
plt.xlabel('Gi (mS/cm^2)')
plt.ylabel('Burst delay (ms)')

plt.show()
