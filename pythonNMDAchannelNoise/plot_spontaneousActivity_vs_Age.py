#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 12:53:52 2018

@author: jingroup

Script plots spontaneous activity versus age
"""
import matplotlib.pyplot as plt

firingRate_Exp = [0.56, 0.32, 0.2, 0.08, 0.01, 0]
age_Exp = [0, 0.05, 0.1, 0.2, 0.3, 0.5]

firingRate_Linear= [0.56, 0.31, 0.19, 0.06, 0.01, 0]
age_Linear = [0, 0.05, 0.1, 0.2, 0.3, 0.5]

plt.figure()

plt.plot(age_Exp, firingRate_Exp, '-bo', label="$\propto$ exp(-age)")
plt.plot(age_Linear, firingRate_Linear, '-ro', label="$\propto$ age")

plt.xlabel('Age')
plt.ylabel('Spontaneous firing rate (Hz)')

plt.legend()

plt.show()

