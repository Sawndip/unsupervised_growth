# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 19:30:06 2018

@author: jingroup

Script checks that all neuron replacement was done correctly
"""
import matplotlib.pyplot as plt
import reading
import os
import numpy as np

starting_trial = 5

dirname = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/maturationTransition3/"

fileCoordBefore = os.path.join(dirname, "RA_xy_" + str(starting_trial) + "beforeReplacement.bin")
fileActiveBefore = os.path.join(dirname, "RA_RA_active_connections_" + str(starting_trial) + "beforeReplacement.bin")
fileSuperBefore = os.path.join(dirname, "RA_RA_super_connections_" + str(starting_trial) + "beforeReplacement.bin")
fileRA2IBefore = os.path.join(dirname, "RA_I_connections_" + str(starting_trial) + "beforeReplacement.bin")
fileI2RABefore = os.path.join(dirname, "I_RA_connections_" + str(starting_trial) + "beforeReplacement.bin")
fileWeightsBefore = os.path.join(dirname, "weights_" + str(starting_trial) + "beforeReplacement.bin")
fileAxonalDelaysRA2IBefore = os.path.join(dirname, "axonal_delays_RA2I_" + str(starting_trial) + "beforeReplacement.bin")
fileAxonalDelaysRA2RABefore = os.path.join(dirname, "axonal_delays_RA2RA_" + str(starting_trial) + "beforeReplacement.bin")
fileAxonalDelaysI2RABefore = os.path.join(dirname, "axonal_delays_I2RA_" + str(starting_trial) + "beforeReplacement.bin")
fileActivityHistoryBefore = os.path.join(dirname, "activity_history_" + str(starting_trial) + "beforeReplacement.bin")
fileReplacementHistoryBefore = os.path.join(dirname, "replacement_history_" + str(starting_trial) + "beforeReplacement.bin")
fileRemodeledIndicatorsBefore = os.path.join(dirname, "remodeled_indicators_" + str(starting_trial) + "beforeReplacement.bin")
fileMatureIndicatorsBefore = os.path.join(dirname, "mature_" + str(starting_trial) + "beforeReplacement.bin")
fileMaturationPropertiesBefore = os.path.join(dirname, "maturation_properties_" + str(starting_trial) + "beforeReplacement.bin")


fileCoordAfter = os.path.join(dirname, "RA_xy_" + str(starting_trial) + "afterReplacement.bin")
fileActiveAfter = os.path.join(dirname, "RA_RA_active_connections_" + str(starting_trial) + "afterReplacement.bin")
fileSuperAfter = os.path.join(dirname, "RA_RA_super_connections_" + str(starting_trial) + "afterReplacement.bin")
fileRA2IAfter = os.path.join(dirname, "RA_I_connections_" + str(starting_trial) + "afterReplacement.bin")
fileI2RAAfter = os.path.join(dirname, "I_RA_connections_" + str(starting_trial) + "afterReplacement.bin")
fileWeightsAfter = os.path.join(dirname, "weights_" + str(starting_trial) + "afterReplacement.bin")
fileAxonalDelaysRA2IAfter = os.path.join(dirname, "axonal_delays_RA2I_" + str(starting_trial) + "afterReplacement.bin")
fileAxonalDelaysRA2RAAfter = os.path.join(dirname, "axonal_delays_RA2RA_" + str(starting_trial) + "afterReplacement.bin")
fileAxonalDelaysI2RAAfter = os.path.join(dirname, "axonal_delays_I2RA_" + str(starting_trial) + "afterReplacement.bin")
fileActivityHistoryAfter = os.path.join(dirname, "activity_history_" + str(starting_trial) + "afterReplacement.bin")
fileReplacementHistoryAfter = os.path.join(dirname, "replacement_history_" + str(starting_trial) + "afterReplacement.bin")
fileRemodeledIndicatorsAfter = os.path.join(dirname, "remodeled_indicators_" + str(starting_trial) + "afterReplacement.bin")
fileMatureIndicatorsAfter = os.path.join(dirname, "mature_" + str(starting_trial) + "afterReplacement.bin")
fileMaturationPropertiesAfter = os.path.join(dirname, "maturation_properties_" + str(starting_trial) + "afterReplacement.bin")

fileReplaced = os.path.join(dirname, "replaced_neurons.bin")

###########################
#### Read data from files
###########################
N = 1000
N_I = 200

replaced_neurons = reading.read_replaced_neurons(fileReplaced)

replaced_neurons_trial = set(replaced_neurons[starting_trial])
not_replaced_trial = [i for i in range(N) if i not in replaced_neurons_trial]

#print replaced_neurons_trial

######### coordinates ##############
coordBefore = reading.read_coordinates(fileCoordBefore)
coordAfter = reading.read_coordinates(fileCoordAfter)

######### active synapses ##############
(_, _, active_synapses_before) = reading.read_synapses(fileActiveBefore)
(_, _, active_synapses_after) = reading.read_synapses(fileActiveAfter)

######### super synapses ###############
(_, _, super_synapses_before) = reading.read_synapses(fileSuperBefore)
(_, _, super_synapses_after) = reading.read_synapses(fileSuperAfter)

############# HVC-RA -> HVC-I connections #############
(_, targets_id_RA2I_before, weights_RA2I_before, \
    syn_lengths_RA2I_before, axonal_delays_RA2I_before) = reading.read_connections(fileRA2IBefore)

(_, targets_id_RA2I_after, weights_RA2I_after, \
    syn_lengths_RA2I_after, axonal_delays_RA2I_after) = reading.read_connections(fileRA2IAfter)

############# HVC-I -> HVC-RA connections #############
(_, targets_id_I2RA_before, weights_I2RA_before, \
    syn_lengths_I2RA_before, axonal_delays_I2RA_before) = reading.read_connections(fileI2RABefore)

(_, targets_id_I2RA_after, weights_I2RA_after, \
    syn_lengths_I2RA_after, axonal_delays_I2RA_after) = reading.read_connections(fileI2RAAfter)

############ HVC-RA -> HVC-RA connections #############
(_, _, weights_before) = reading.read_weights(fileWeightsBefore) 
(_, _, weights_after) = reading.read_weights(fileWeightsAfter)

########## Axonal delays ######################

(_, _, axonal_delays_I2RA_from_delay_file_before) = reading.read_axonal_delays(fileAxonalDelaysI2RABefore) 
(_, _, axonal_delays_I2RA_from_delay_file_after) = reading.read_axonal_delays(fileAxonalDelaysI2RAAfter) 

(_, _, axonal_delays_RA2RA_from_delay_file_before) = reading.read_axonal_delays(fileAxonalDelaysRA2RABefore) 
(_, _, axonal_delays_RA2RA_from_delay_file_after) = reading.read_axonal_delays(fileAxonalDelaysRA2RAAfter) 

(_, _, axonal_delays_RA2I_from_delay_file_before) = reading.read_axonal_delays(fileAxonalDelaysRA2IBefore) 
(_, _, axonal_delays_RA2I_from_delay_file_after) = reading.read_axonal_delays(fileAxonalDelaysRA2IAfter) 


################## Activity History ###################
(_, _, activity_history_before) = reading.read_activity_history(fileActivityHistoryBefore) 
(_, _, activity_history_after) = reading.read_activity_history(fileActivityHistoryAfter) 

################# Replacement History ##################
(_, _, replacement_history_before) = reading.read_replacement_history(fileReplacementHistoryBefore) 
(_, _, replacement_history_after) = reading.read_replacement_history(fileReplacementHistoryAfter) 

################# Axon-remodeling indicators ###########
(_, _, remodeled_indicators_before) = reading.read_remodeled_indicators(fileRemodeledIndicatorsBefore) 
(_, _, remodeled_indicators_after) = reading.read_remodeled_indicators(fileRemodeledIndicatorsAfter) 

################# Maturation indicators ###########
#(_, _, mature_indicators_before) = reading.read_mature_indicators(fileMatureIndicatorsBefore)
#(_, _, mature_indicators_after) = reading.read_mature_indicators(fileMatureIndicatorsAfter)

################# Maturation properties ###########
(_, _, mature_indicators_before, maturation_rate_before, Erest_before, GCa_before) = reading.read_maturation_properties(fileMaturationPropertiesBefore)
(_, _, mature_indicators_after, maturation_rate_after, Erest_after, GCa_after) = reading.read_maturation_properties(fileMaturationPropertiesAfter)


#print active_synapses

#print "Active before reading: ",active_synapses_before
#print "Active after reading: ",active_synapses_after

#######################################################
#### Check that data before and after reading matches
#######################################################

print weights_before[24]
print weights_after[24]

print activity_history_before[24]

print activity_history_after[24]

#print replacement_history_before
#print replacement_history_after

#### check that neurons that changed their coordinates are those and only those that were replaced ###
set_diffCoord = set(list(np.where(np.all(coordBefore != coordAfter, axis=1))[0]))

print set_diffCoord == replaced_neurons_trial

#### check that replaced neurons have zero active and super outputs and not replaced have no connections
#### to replaced

for i in replaced_neurons_trial:
    if len(active_synapses_after[i]) > 0:
        print "Neuron {0} has active outputs".format(i)
    
    if len(super_synapses_after[i]) > 0:
        print "Neuron {0} has super outputs".format(i)
        
print np.where(weights_after[list(replaced_neurons_trial)] > 0)[0]
print np.where(weights_after[not_replaced_trial][:,list(replaced_neurons_trial)] > 0)[0]
print np.array_equal(weights_before[not_replaced_trial][:,not_replaced_trial], weights_after[not_replaced_trial][:,not_replaced_trial])
 

for i in not_replaced_trial:
    for active in active_synapses_after[i]:
        if active in replaced_neurons_trial:
            print "Neuron {0} has active output to replaced neuron".format(active)

    for super in super_synapses_after[i]:
        if super in replaced_neurons_trial:
            print "Neuron {0} has active output to replaced neuron".format(super)

    # check that difference in active and super connections is only in replaced neurons

    difference = set(active_synapses_before[i]).symmetric_difference(set(active_synapses_after[i]))
    
    for j in difference:
        if j not in replaced_neurons_trial:
            print "Active synapses for neuron {0} differ more than in replaced neurons".format(i)
            
    difference = set(super_synapses_before[i]).symmetric_difference(set(super_synapses_after[i]))
    
    for j in difference:
        if j not in replaced_neurons_trial:
            print "Super synapses for neuron {0} differ more than in replaced neurons".format(i)
            
#### check that not replaced neurons have same connections to interneurons
print np.array_equal(np.array(targets_id_RA2I_before)[not_replaced_trial], np.array(targets_id_RA2I_after)[not_replaced_trial])
print np.array_equal(np.array(weights_RA2I_before)[not_replaced_trial], np.array(weights_RA2I_after)[not_replaced_trial])
print np.array_equal(np.array(syn_lengths_RA2I_before)[not_replaced_trial], np.array(syn_lengths_RA2I_after)[not_replaced_trial])
print np.array_equal(np.array(axonal_delays_RA2I_before)[not_replaced_trial], np.array(axonal_delays_RA2I_after)[not_replaced_trial])

#### and replaced neurons have different connections to interneurons
print np.array_equal(np.array(targets_id_RA2I_before)[list(replaced_neurons_trial)], np.array(targets_id_RA2I_after)[list(replaced_neurons_trial)])
print np.array_equal(np.array(weights_RA2I_before)[list(replaced_neurons_trial)], np.array(weights_RA2I_after)[list(replaced_neurons_trial)])
print np.array_equal(np.array(syn_lengths_RA2I_before)[list(replaced_neurons_trial)], np.array(syn_lengths_RA2I_after)[list(replaced_neurons_trial)])
print np.array_equal(np.array(axonal_delays_RA2I_before)[list(replaced_neurons_trial)], np.array(axonal_delays_RA2I_after)[list(replaced_neurons_trial)])

#### check that connections from interneurons to HVC-RA differ only in replaced neurons 
for i in range(N_I):
    difference = set(targets_id_I2RA_before[i]).symmetric_difference(set(targets_id_I2RA_after[i]))
    
    for j in difference:
        if j not in replaced_neurons_trial:
            print "Connections from interneuron {0} differ more than in replaced neurons".format(i)


### check that histories did not change for not replaced neurons and some histories
### changes for replaced neurons

print np.array_equal(activity_history_before[not_replaced_trial], activity_history_after[not_replaced_trial])
print np.logical_and.reduce(np.logical_and.reduce(activity_history_after[list(replaced_neurons_trial)] == 0))

print np.array_equal(replacement_history_before[not_replaced_trial], replacement_history_after[not_replaced_trial])
print np.logical_and.reduce(replacement_history_after[list(replaced_neurons_trial)] == 0)

print np.array_equal(remodeled_indicators_before[not_replaced_trial], remodeled_indicators_after[not_replaced_trial])
print np.logical_and.reduce(remodeled_indicators_after[list(replaced_neurons_trial)] == 0)

print np.array_equal(mature_indicators_before[not_replaced_trial], mature_indicators_after[not_replaced_trial])
print np.logical_and.reduce(mature_indicators_after[list(replaced_neurons_trial)] == 0)

MATURATION_RATE_SPONTANEOUS = 100000
E_REST_IMMATURE = -55.0
GCA_IMMATURE = 0.0

print np.array_equal(maturation_rate_before[not_replaced_trial], maturation_rate_after[not_replaced_trial])
print np.logical_and.reduce(maturation_rate_after[list(replaced_neurons_trial)] == MATURATION_RATE_SPONTANEOUS)

print np.array_equal(Erest_before[not_replaced_trial], Erest_after[not_replaced_trial])
print np.logical_and.reduce(Erest_after[list(replaced_neurons_trial)] == E_REST_IMMATURE)

print np.array_equal(GCa_before[not_replaced_trial], GCa_after[not_replaced_trial])
print np.logical_and.reduce(GCa_after[list(replaced_neurons_trial)] == GCA_IMMATURE)
