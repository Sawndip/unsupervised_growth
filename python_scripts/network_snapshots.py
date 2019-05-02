# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 14:04:44 2017

@author: jingroup

Script generates pajek net files with neurons fired during different frames
"""
import reading
import mature_chain_analysis as analysis
import os

SIDE = 100.0

def write_snapshot(chain_neurons, group, (xx,yy), filename):
    """
    Generate pajek net file with network snapshot
    
    Input: 
    chain_neurons: neurons in the chain;
    group: neurons that fired in the same frame
    (xx, yy): coordinates of HVC(RA) neurons
    filename: filename to which write the data
    """
    with open(filename, 'w') as f:
        f.write("*Vertices %s \n" % len(chain_neurons))
    
        for i, n in enumerate(chain_neurons):
            if n in group:   
                f.write('{0} "{1}" {2} {3} ic Black\n'.format(i+1, n, xx[n] / SIDE, yy[n] / SIDE))
            else:
                f.write('{0} "{1}" {2} {3} ic White\n'.format(i+1, n, xx[n] / SIDE, yy[n] / SIDE))
            
        
        f.write("*Arcs\n")
        	
    

dirname = "/home/eugene/results/noDelays/replacement/dispersed/"

simname = "190617_lionx_2"
trial_number = 181000

outdir = os.path.join(dirname, "matureTest/" + simname + "/snapshots")

file_chain_test = os.path.join(dirname, "matureTest/" + simname + "/mature_chain_test.bin")
file_RA_super = os.path.join(dirname, simname + "/RA_RA_super_connections_" + str(trial_number) + "_.bin")
file_RA_xy = os.path.join(dirname, simname + "/RA_xy_" + str(trial_number) + "_.bin")

_RA, num_trials, firing_robustness, average_num_dend_spikes_in_trial, \
average_num_soma_spikes_in_trial, mean_burst_time, std_burst_time = reading.read_chain_test(file_chain_test)

mean_burst_time_stronglyConnected, id_stronglyConnected, firing_robustness_stronglyConnected, \
        std_burst_time_stronglyConnected, num_soma_spikes_stronglyConnected, \
        num_dend_spikes_stronglyConnected = analysis.get_bursts_of_chain_neurons(mean_burst_time, std_burst_time,
                                firing_robustness, average_num_soma_spikes_in_trial, average_num_dend_spikes_in_trial, file_RA_super)



groups, burst_times_in_groups = analysis.split_in_frames(mean_burst_time_stronglyConnected, id_stronglyConnected)




print groups
print burst_times_in_groups

print sum([len(group) for group in groups])
print len(id_stronglyConnected)

coordinates = reading.read_coordinates(file_RA_xy)

for i, group in enumerate(groups):
    filename = os.path.join(outdir, "snapshot" + str(i+1) + ".net")
    write_snapshot(id_stronglyConnected, group, coordinates, filename)
