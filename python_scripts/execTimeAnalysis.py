# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 17:01:29 2017

@author: jingroup

Script analyzes execution time results
"""
import reading
import os
import matplotlib.pyplot as plt

dirNoSave = '/home/eugene/Output/timing/npScaling/noSave/'
dirSave = '/home/eugene/Output/timing/npScaling/save/np/'


def get_execution_time_data(dirname):
    """
    Reads all files for dirname containing info on trial execution time
    """
    names = os.listdir(dirname)
    
    num_processes = []
    execution_time = []
    
    for n in names:
        filename = os.path.join(dirname, n)
        
        (N_RA, N_I, np, timestep, network_update_frequency, average_execution_time) = reading.read_execution_time(filename)
        
        num_processes.append(np)
        execution_time.append(average_execution_time)
        

    num_processes, execution_time = zip(*sorted(zip(num_processes, execution_time)))
    
    return (N_RA, N_I, num_processes, timestep, network_update_frequency, execution_time)
 
(N_RA, N_I, num_processes, timestep, network_update_frequency, execution_time_no_save) = get_execution_time_data(dirNoSave)
(N_RA, N_I, num_processes, timestep, network_update_frequency, execution_time_save) = get_execution_time_data(dirSave)

print network_update_frequency

f = plt.figure()
ax = f.add_subplot(111)
ax.set_title("Average execution time of one trial")
ax.plot(num_processes, execution_time_no_save, 'r', label="no writing")
ax.plot(num_processes, execution_time_save, 'b', label="with writing")

ax.set_xlabel("# processes")
ax.set_ylabel("execution time (s)")
#ax.set_yscale("log")
#ax.set_xscale("log")

plt.legend()


# analyze network update data
dirname = "/home/eugene/Output/timing/netUpdateScaling_2/"

names = os.listdir(dirname)
    
net_update_frequency = []

execution_time = []
    
for n in names:
    filename = os.path.join(dirname, n)
        
    (N_RA, N_I, np, timestep, network_update_frequency, average_execution_time) = reading.read_execution_time(filename)
    
    #print network_update_frequency 
    #print timestep
    #print N_RA
    
    net_update_frequency.append(network_update_frequency)
    execution_time.append(average_execution_time)
        

net_update_frequency, execution_time = zip(*sorted(zip(net_update_frequency, execution_time)))

print net_update_frequency

f = plt.figure()
ax = f.add_subplot(111)
ax.set_title("Average execution time of one trial")
ax.plot(net_update_frequency, execution_time, 'r')


ax.set_xlabel("network update frequency (ms)")
ax.set_ylabel("execution time (s)")


plt.show()