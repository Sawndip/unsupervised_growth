# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 20:41:16 2018

@author: jingroup

Script analyzes HVC-RA neuron activity in the network
"""
import matplotlib.pyplot as plt
import reading
import os
import numpy as np
import utils


dirname = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays/"
trial_number = 450
fileActivityHistory = os.path.join(dirname, "activity_history_" + str(trial_number) + ".bin")
fileTraining = os.path.join(dirname, "training_neurons.bin")

(_, _, activity_history) = reading.read_activity_history(fileActivityHistory) 
training_neurons = reading.read_training_neurons(fileTraining)

np.set_printoptions(threshold=np.nan)

#print activity_history[0][0:10]
#print activity_history[0:100,0:100]

#print activity_history

# check that activity is reasonable
print "Negative numbers in activity: ",np.any(activity_history < 0)
print "Large positive numbers in activity: ",np.any(activity_history > 10)


num_neurons = activity_history.shape[0]
history_size = activity_history.shape[1]
window_size = 20
FIRING_RATE_THRESHOLD = 0.1
SPIKE_PROBABILITY_THRESHOLD = 0.1

fire_indicators = np.ones_like(activity_history, np.float32)
fire_indicators[activity_history == 0] = 0

#print fire_indicators[156]
#==============================================================================
# if trial_number >= history_size:    
#     num_points = trial_number / window_size
# else:
#     num_points = history_size / window_size
#==============================================================================
    
firing_rate = np.empty((num_neurons, history_size / window_size), np.float32)
probability_to_spike = np.empty((num_neurons, history_size / window_size), np.float32)

time = np.array([trial_number - window_size/2 - i*window_size for i in range(history_size / window_size)])

for i in range(num_neurons):
    for j in range(history_size / window_size):
        firing_rate[i][j] = float(np.sum(activity_history[i][j*window_size:(j+1)*window_size])) / float(window_size)
        probability_to_spike[i][j] = float(np.sum(fire_indicators[i][j*window_size:(j+1)*window_size])) / float(window_size)
        
average_firing_rate_in_previous_trials = np.empty(num_neurons, np.float32)
average_probability_to_spike_in_previous_trials = np.empty(num_neurons, np.float32)

if trial_number < history_size:
    for i in range(num_neurons):
        average_firing_rate_in_previous_trials[i] = np.mean(firing_rate[i][0:trial_number/window_size])
        average_probability_to_spike_in_previous_trials[i] = np.mean(probability_to_spike[i][0:trial_number/window_size])
else:
    for i in range(num_neurons):
        average_firing_rate_in_previous_trials[i] = np.mean(firing_rate[i][:])
        average_probability_to_spike_in_previous_trials[i] = np.mean(probability_to_spike[i][:])

pool_neurons = [i for i in range(num_neurons) if i not in training_neurons]

population_firing_rate_of_pool_neurons = firing_rate[pool_neurons,:].sum(axis=0)

plt.figure()

plt.plot(time, population_firing_rate_of_pool_neurons)
plt.xlabel('Time (# of trials)')
plt.ylabel('Population firing rate (# of spikes per trial)')
plt.xlim([0, trial_number + 10])
#==============================================================================
# if trial_number < history_size:
#     for i in range(num_neurons):
#         average_firing_rate_in_previous_trials[i] = np.mean(activity_history.astype(float)[i][0:trial_number])
#         average_probability_to_spike_in_previous_trials[i] = np.mean(fire_indicators[i][0:trial_number])
# else:
#     for i in range(num_neurons):
#         average_firing_rate_in_previous_trials[i] = np.mean(activity_history.astype(float)[i][:])
#         average_probability_to_spike_in_previous_trials[i] = np.mean(fire_indicators[i][:])
#==============================================================================

#print average_firing_rate_in_previous_trials[156]
#print average_firing_rate_in_previous_trials[0]

#print average_firing_rate_in_previous_trials

active_neurons = np.where(average_probability_to_spike_in_previous_trials > SPIKE_PROBABILITY_THRESHOLD)[0]
active_neurons_not_training = [i for i in active_neurons if i not in training_neurons]



#print training_neurons
print active_neurons_not_training

#print activity_history[:][2496]
#print activity_history[4][:]

max_num_neurons_to_plot = 10

if ( len(active_neurons_not_training) > 0 ):
    if ( len(active_neurons_not_training) >= max_num_neurons_to_plot ):
        num_neurons_to_plot = max_num_neurons_to_plot
    else:
        num_neurons_to_plot = len(active_neurons_not_training)
    
    print "Neurons plotted: ",active_neurons_not_training[:num_neurons_to_plot]

    f, axarr = plt.subplots(num_neurons_to_plot, sharex=True)
    for i in range(num_neurons_to_plot):
        axarr[i].plot(time, firing_rate[active_neurons_not_training[i]], label='neuron {0}'.format(active_neurons_not_training[i]))
        
    axarr[0].set_xlim([0, trial_number + 10])
    
    # add a big axes, hide frame
    f.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    plt.xlabel("Time in trial numbers")
    plt.ylabel("Firing rate (# of spikes per trial)")
    
    f, axarr = plt.subplots(num_neurons_to_plot, sharex=True)
    for i in range(num_neurons_to_plot):
        axarr[i].plot(time, probability_to_spike[active_neurons_not_training[i]], label='neuron {0}'.format(active_neurons_not_training[i]))
        #axarr[i].plot(time, fire_indicators[active_neurons_not_training[i]], label='neuron {0}'.format(active_neurons_not_training[i]))
        axarr[i].legend(loc=2)
        axarr[i].set_ylim([-0.1,1.1])
    axarr[0].set_xlim([0, trial_number + 10])
    
    # add a big axes, hide frame
    f.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    plt.xlabel("Time in trial numbers")
    plt.ylabel("Probability to spike")
    
    trialStep = 50
    
    input_weights = utils.get_input_weight_time_sequence(dirname, trial_number, trialStep, training_neurons, active_neurons_not_training[:num_neurons_to_plot])
    
    time_for_weights = [i*trialStep for i in range(input_weights.shape[1])]
    
    f, axarr = plt.subplots(num_neurons_to_plot, sharex=True)
    for i in range(num_neurons_to_plot):
        axarr[i].plot(time_for_weights, input_weights[i], label='neuron {0}'.format(active_neurons_not_training[i]))
        #axarr[i].plot(time, fire_indicators[active_neurons_not_training[i]], label='neuron {0}'.format(active_neurons_not_training[i]))
        axarr[i].legend(loc=2)
        #axarr[i].set_ylim([-0.1,1.1])
    axarr[0].set_xlim([0, trial_number + 10])
    
    # add a big axes, hide frame
    f.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    plt.xlabel("Time in trial numbers")
    plt.ylabel("Input weights from training neurons")
    
    plt.show()