# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 16:18:21 2017

@author: jingroup

Script analyzes spike patterns of mature chain 
"""
import matplotlib.pyplot as plt
import reading
import numpy as np
from matplotlib import cm
import matplotlib as mpl
from matplotlib import ticker



file_chain_test = "/home/eugene/results/noDelays/replacement/clustered/matureTest/190617_lionx_3/mature_chain_test.bin"
file_RA_super = "/home/eugene/results/noDelays/replacement/clustered/190617_lionx_3/RA_RA_super_connections_189000_.bin"

CURRENT_INJECTION_TIME = 100
ROBUST_FIRING_THRESHOLD = 0.25

N_RA, num_trials, firing_robustness, average_num_dend_spikes_in_trial, \
average_num_soma_spikes_in_trial, mean_burst_time, std_burst_time = reading.read_chain_test(file_chain_test)

print N_RA
print num_trials
#print firing_robustness
#print average_num_dend_spikes_in_trial
#print average_num_soma_spikes_in_trial
#print mean_burst_time
#print std_burst_time


(N_RA, RA_super_targets, RA_super_targets_G) = reading.read_connections(file_RA_super)

# select neurons connected by supersynapses

superTargets = set([element for sublist in RA_super_targets for element in sublist])
superSources = set([i for i in xrange(len(RA_super_targets)) if len(RA_super_targets[i]) > 0])
stronglyConnected = superTargets | superSources

num_soma_spikes_stronglyConnected = []
num_dend_spikes_stronglyConnected = []
mean_burst_time_stronglyConnected = []
std_burst_time_stronglyConnected = []
firing_robustness_stronglyConnected = []
id_stronglyConnected = []
id_stronglyConnectedNotFired = []

for i in xrange(N_RA):
    if i in stronglyConnected:
        if mean_burst_time[i] != -1:
            id_stronglyConnected.append(i)
            firing_robustness_stronglyConnected.append(firing_robustness[i]) 
            mean_burst_time_stronglyConnected.append(mean_burst_time[i] - CURRENT_INJECTION_TIME) 
            std_burst_time_stronglyConnected.append(std_burst_time[i])       
            num_soma_spikes_stronglyConnected.append(average_num_soma_spikes_in_trial[i])
            num_dend_spikes_stronglyConnected.append(average_num_dend_spikes_in_trial[i])
            
        else:
            id_stronglyConnectedNotFired.append(i)


mean_burst_time_stronglyConnected, id_stronglyConnected, firing_robustness_stronglyConnected, \
std_burst_time_stronglyConnected, num_soma_spikes_stronglyConnected, \
num_dend_spikes_stronglyConnected = \
    zip(*sorted(zip(mean_burst_time_stronglyConnected, id_stronglyConnected, \
    firing_robustness_stronglyConnected, std_burst_time_stronglyConnected, num_soma_spikes_stronglyConnected, \
    num_dend_spikes_stronglyConnected)))

# plot spike times for strongly connected neurons
f1 = plt.figure()
ax1 = f1.add_subplot(111)

# normalize colors so that min firing rate is at the bottom and 1.0 is on top
norm = mpl.colors.Normalize(vmin=ROBUST_FIRING_THRESHOLD, vmax=1.0)
cmap = cm.coolwarm

m = cm.ScalarMappable(norm=norm, cmap=cmap)

colors = m.to_rgba(np.array(firing_robustness_stronglyConnected))
# plot fake scatterplot to conveniently attach colorbar
plot = ax1.scatter(firing_robustness_stronglyConnected, firing_robustness_stronglyConnected, c=firing_robustness_stronglyConnected, cmap = 'coolwarm')
cb = f1.colorbar(plot)
plot.set_clim(vmin=ROBUST_FIRING_THRESHOLD, vmax=1.0)
ax1.clear()



# (generate plot here)
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()

for i in range(len(mean_burst_time_stronglyConnected)):
    sigma = std_burst_time_stronglyConnected[i]
    mean = mean_burst_time_stronglyConnected[i]
    
    ax1.vlines(mean, i-0.5, i+0.5, lw=2, color=colors[i], cmap = 'coolwarm')
    ax1.hlines(i, mean - sigma, mean + sigma)
    ax1.vlines(mean-sigma, i-0.25, i+0.25)
    ax1.vlines(mean+sigma, i-0.25, i+0.25)

#plot.set_clim(vmin=0.9, vmax=1)        
plt.yticks(xrange(len(mean_burst_time_stronglyConnected)), id_stronglyConnected)
ax1.set_xlim([-5, mean_burst_time_stronglyConnected[-1] + max(std_burst_time_stronglyConnected) + 5])
ax1.set_ylabel("neuron")
ax1.set_xlabel("spike time(ms)")
ax1.set_title("Statistics of spike times for strongly connected neurons")
#ax1.grid(True)

# plot jitter for strongly connected neurons
f2 = plt.figure()
ax2 = f2.add_subplot(111)

ax2.plot(range(0, len(mean_burst_time_stronglyConnected)), std_burst_time_stronglyConnected)
ax2.set_title("Jitter for strongly connected neurons")
ax2.set_ylabel("standard deviation (ms)")
ax2.set_xlabel("neuron ID along the chain")

# plot average number of spikes produced
f3 = plt.figure()
plt.suptitle("Average number of spikes in trial for strongly connected neurons")
ax3 = f3.add_subplot(121)

ax3.plot(range(0, len(num_soma_spikes_stronglyConnected)), num_soma_spikes_stronglyConnected)
ax3.set_title("Somatic spikes")
ax3.set_ylabel("# spikes")
ax3.set_xlabel("neuron ID along the chain")
ax3.set_ylim([min(num_soma_spikes_stronglyConnected) - 1, max(num_soma_spikes_stronglyConnected) + 1])

ax4 = f3.add_subplot(122)

ax4.plot(range(0, len(num_dend_spikes_stronglyConnected)), num_dend_spikes_stronglyConnected)
ax4.set_title("Dendritic spikes")
ax4.set_ylabel("# spikes")
ax4.set_xlabel("neuron ID along the chain")
ax4.set_ylim([min(num_dend_spikes_stronglyConnected) - 1, max(num_dend_spikes_stronglyConnected) + 1])


# select neurons with high firing robustness
num_soma_spikes_robust = []
num_dend_spikes_robust = []
mean_burst_time_robust = []
std_burst_time_robust = []
firing_robustness_robust = []
id_robust = []

for i in xrange(N_RA):
    if firing_robustness[i] >= ROBUST_FIRING_THRESHOLD:
        id_robust.append(i)
        firing_robustness_robust.append(firing_robustness[i]) 
        mean_burst_time_robust.append(mean_burst_time[i] - CURRENT_INJECTION_TIME) 
        std_burst_time_robust.append(std_burst_time[i])  
        num_soma_spikes_robust.append(average_num_soma_spikes_in_trial[i])
        num_dend_spikes_robust.append(average_num_dend_spikes_in_trial[i])

id_robust_not_stronglyConnected = [] # neurons that fire robustly but are not strongly connected
      
for i in id_robust:
    if i not in id_stronglyConnected:
        id_robust_not_stronglyConnected.append(i)
        
# if there are neurons that fire robustly but are not strongly connected, plot 
# spikes times of all robust neurons
if len(id_robust_not_stronglyConnected) > 0:
    print "Neurons that fire robustly, but are not strongly connected: ",id_robust_not_stronglyConnected
                
    mean_burst_time_robust, id_robust, firing_robustness_robust, \
    std_burst_time_robust, num_soma_spikes_robust, \
    num_dend_spikes_robust = \
    zip(*sorted(zip(mean_burst_time_robust, id_robust, \
    firing_robustness_robust, std_burst_time_robust, num_soma_spikes_robust, \
    num_dend_spikes_robust)))
    
    f4 = plt.figure()
    ax5 = f4.add_subplot(111)

    # normalize colors so that min firing rate is at the bottom and 1.0 is on top
    norm = mpl.colors.Normalize(vmin=min(firing_robustness_robust), vmax=1.0)
    cmap = cm.coolwarm
    
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    
    colors = m.to_rgba(np.array(firing_robustness_robust))
    
    plot = ax5.scatter(firing_robustness_robust, firing_robustness_robust, c = firing_robustness_robust, cmap = 'coolwarm')
    ax5.clear()
    f4.colorbar(plot)
    
    for i in range(len(mean_burst_time_robust)):
        if (mean_burst_time_robust[i] > 0 ):
            sigma = std_burst_time_robust[i]
            mean = mean_burst_time_robust[i]
            
            ax3.vlines(mean, i-0.5, i+0.5, lw=2, color=colors[i])
            ax3.hlines(i, mean - sigma, mean + sigma)
            ax3.vlines(mean-sigma, i-0.25, i+0.25)
            ax3.vlines(mean+sigma, i-0.25, i+0.25)
    
    plt.yticks(xrange(len(mean_burst_time_robust)), id_robust)
    ax5.set_xlim([-5, mean_burst_time_robust[-1] + max(std_burst_time_robust) + 5])
    ax5.set_ylabel("neuron")
    ax5.set_xlabel("spike time(ms)")
    ax5.set_title("Statistics of spike times for robustly firing neurons")
    ax5.grid(True)

    # plot jitter for robustly firing neurons
    f5 = plt.figure()
    ax6 = f5.add_subplot(111)
    
    ax6.plot(range(0, len(mean_burst_time_robust)), std_burst_time_robust)
    ax6.set_title("Jitter for robustly connected neurons")
    ax6.set_ylabel("standard deviation (ms)")
    ax6.set_xlabel("neuron ID along the chain")
    
    # plot average number of spikes produced
    f6 = plt.figure()
    plt.suptitle("Average number of spikes in trial for robustly firing neurons")
    ax7 = f6.add_subplot(121)
    
    ax7.plot(range(0, len(num_soma_spikes_robust)), num_soma_spikes_robust)
    ax7.set_title("Somatic spikes")
    ax7.set_ylabel("# spikes")
    ax7.set_xlabel("neuron ID along the chain")
    ax7.set_ylim([min(num_soma_spikes_robust) - 1, max(num_soma_spikes_robust) + 1])
    
    ax8 = f6.add_subplot(122)
    
    ax8.plot(range(0, len(num_dend_spikes_robust)), num_dend_spikes_robust)
    ax8.set_title("Dendritic spikes")
    ax8.set_ylabel("# spikes")
    ax8.set_xlabel("neuron ID along the chain")
    ax8.set_ylim([min(num_dend_spikes_robust) - 1, max(num_dend_spikes_robust) + 1])
    
    
plt.show()
