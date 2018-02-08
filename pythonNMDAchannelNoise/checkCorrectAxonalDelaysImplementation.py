# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 16:45:16 2018

@author: jingroup

Script compares results for trial with no axonal delays implemented and
trial with axonal delays set to zero
"""
import reading
import numpy as np

fileDend_noDelays = "/home/eugene/Output/networks/chainGrowth/testGrowthNoDelays/spike_times_dend.bin"
fileSoma_noDelays = "/home/eugene/Output/networks/chainGrowth/testGrowthNoDelays/spike_times_soma.bin"
fileInterneuron_noDelays = "/home/eugene/Output/networks/chainGrowth/testGrowthNoDelays/spike_times_interneuron.bin"
fileWeights_noDelays = "/home/eugene/Output/networks/chainGrowth/testGrowthNoDelays/weights_10.bin"
fileActiveSynapses_noDelays = "/home/eugene/Output/networks/chainGrowth/testGrowthNoDelays/RA_RA_active_connections.bin"

fileDend_zeroDelays = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays/spike_times_dend.bin"
fileSoma_zeroDelays = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays/spike_times_soma.bin"
fileInterneuron_zeroDelays = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays/spike_times_interneuron.bin"
fileWeights_zeroDelays = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays/weights_10.bin"
fileActiveSynapses_zeroDelays = "/home/eugene/Output/networks/chainGrowth/testGrowthDelays/RA_RA_active_connections.bin"

(_, _, spike_times_dend_noDelays, neuron_fired_dend_noDelays) = reading.read_time_info(fileDend_noDelays)
(_, _, spike_times_soma_noDelays, neuron_fired_soma_noDelays) = reading.read_time_info(fileSoma_noDelays)
(_, _, spike_times_interneuron_noDelays, interneuron_fired_noDelays) = reading.read_time_info(fileInterneuron_noDelays)

(_, _, spike_times_dend_zeroDelays, neuron_fired_dend_zeroDelays) = reading.read_time_info(fileDend_zeroDelays)
(_, _, spike_times_soma_zeroDelays, neuron_fired_soma_zeroDelays) = reading.read_time_info(fileSoma_zeroDelays)
(_, _, spike_times_interneuron_zeroDelays, interneuron_fired_zeroDelays) = reading.read_time_info(fileInterneuron_zeroDelays)

(_, _, weights_noDelays) = reading.read_weights(fileWeights_noDelays) 
(_, _, weights_zeroDelays) = reading.read_weights(fileWeights_zeroDelays)

(_, _, active_synapses_noDelays) = reading.read_synapses(fileActiveSynapses_noDelays)
(_, _, active_synapses_zeroDelays) = reading.read_synapses(fileActiveSynapses_zeroDelays)

print np.array_equal(np.array(spike_times_dend_noDelays), np.array(spike_times_dend_zeroDelays))
print np.array_equal(np.array(neuron_fired_dend_noDelays), np.array(neuron_fired_dend_zeroDelays))

print np.array_equal(np.array(spike_times_soma_noDelays), np.array(spike_times_soma_zeroDelays))
print np.array_equal(np.array(neuron_fired_soma_noDelays), np.array(neuron_fired_soma_zeroDelays))

print np.array_equal(np.array(spike_times_interneuron_noDelays), np.array(spike_times_interneuron_zeroDelays))
print np.array_equal(np.array(interneuron_fired_noDelays), np.array(interneuron_fired_zeroDelays))

print np.array_equal(weights_noDelays, weights_zeroDelays)

print np.array_equal(np.array(active_synapses_noDelays), np.array(active_synapses_zeroDelays))

#print active_synapses_noDelays
#print weights_noDelays[860][850:1000]

#print weights_zeroDelays[860][850:1000]


#print spike_times_dend_noDelays
#print spike_times_dend_zeroDelays
