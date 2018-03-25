# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 09:37:24 2018

@author: jingroup

Script plots input weights to certain neurons
"""
import reading
import matplotlib.pyplot as plt
import utils

dirname = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/noImmatureOut2/"
trial_number = 2500
trialStep = 50

source_neurons = [125]
target_neurons = [897, 803, 465, 203, 913, 210, 186, 785, 125, 927]

input_weights = utils.get_source_input_weight_time_sequence(dirname, trial_number, trialStep, source_neurons, target_neurons)


num_targets = len(target_neurons)

time = [i*trialStep for i in range(trial_number / trialStep +1)]

f, axarr = plt.subplots(num_targets, sharex=True)
#axarr[i].set_ylim([-0.1,1.1])

for i in range(num_targets):
    axarr[i].plot(time, input_weights[i], label='target {0}'.format(target_neurons[i]))
    axarr[i].legend(loc=2)
    
axarr[num_targets - 1].set_xlabel('time in # of trials')
    
    
    
plt.show()