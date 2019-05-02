#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 16 13:58:01 2018

@author: jingroup
"""
import matplotlib.pyplot as plt
import numpy as np

count_1 = [3995, 12903, 15566, 13271, 11022, 8581, 6182, 4199, 2466, 785]
label_1 = "matTrans31"

count_2 = [3640, 10008, 11240, 9732, 8362, 6929, 5438, 4122, 2827, 1699, 694, 11]
label_2 = "matTrans36"


plt.figure()
plt.plot(range(1,len(count_1)+1), np.array(count_1).astype(float) / float(sum(count_1)), label=label_1)
plt.plot(range(1,len(count_2)+1), np.array(count_2).astype(float) / float(sum(count_2)), label=label_2)

plt.legend()
plt.xlabel('Distance')
plt.ylabel('Norm count')

plt.show()