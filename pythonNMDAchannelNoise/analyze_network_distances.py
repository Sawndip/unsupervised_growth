#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 13:43:36 2018

@author: jingroup

Script plots distribution of shortest distances
"""
import numpy as np
import matplotlib.pyplot as plt

num_pairs_cont = [4166, 11233, 12271, 10807, 9769, 8541, 7534, 6432, 5489, 4360, 3244, 2108, 953, 58]
num_pairs_discont = [1875, 2659, 2679, 2252, 2069, 1824, 1523, 1338, 1026, 694, 397, 81]

num_pairs_cont_subset = [1776, 4321, 4115, 2828, 1814, 761, 61]


plt.figure()

plt.plot(range(1,len(num_pairs_cont_subset)+1), np.array(num_pairs_cont_subset).astype(float) / np.sum(num_pairs_cont_subset), label='continuous')
plt.plot(range(1,len(num_pairs_discont)+1), np.array(num_pairs_discont).astype(float) / np.sum(num_pairs_discont), label='discontinuous')

plt.xlabel('Shortest distance')
plt.ylabel('Fraction of paths')

plt.legend()

plt.show()
