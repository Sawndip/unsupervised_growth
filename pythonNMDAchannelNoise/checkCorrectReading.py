# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 19:30:06 2018

@author: jingroup

Script checks that all data files were correctly read by network:
Files after reading match the ones before.
"""
import matplotlib.pyplot as plt
import reading
import os
import numpy as np

starting_trial = 200

dirname_before = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/maturationTransition3/"
dirname_after = "/home/eugene/Output/networks/chainGrowth/passiveDendrite/maturationTransition3/"


fileCoordBefore = os.path.join(dirname_before, "RA_xy_" + str(starting_trial) + ".bin")
fileActiveBefore = os.path.join(dirname_before, "RA_RA_active_connections_" + str(starting_trial) + ".bin")
fileSuperBefore = os.path.join(dirname_before, "RA_RA_super_connections_" + str(starting_trial) + ".bin")
fileRA2IBefore = os.path.join(dirname_before, "RA_I_connections_" + str(starting_trial) + ".bin")
fileI2RABefore = os.path.join(dirname_before, "I_RA_connections_" + str(starting_trial) + ".bin")
fileWeightsBefore = os.path.join(dirname_before, "weights_" + str(starting_trial) + ".bin")
fileAxonalDelaysRA2IBefore = os.path.join(dirname_before, "axonal_delays_RA2I_" + str(starting_trial) + ".bin")
fileAxonalDelaysRA2RABefore = os.path.join(dirname_before, "axonal_delays_RA2RA_" + str(starting_trial) + ".bin")
fileAxonalDelaysI2RABefore = os.path.join(dirname_before, "axonal_delays_I2RA_" + str(starting_trial) + ".bin")
fileActivityHistoryBefore = os.path.join(dirname_before, "activity_history_" + str(starting_trial) + ".bin")
fileReplacementHistoryBefore = os.path.join(dirname_before, "replacement_history_" + str(starting_trial) + ".bin")
fileRemodeledIndicatorsBefore = os.path.join(dirname_before, "remodeled_indicators_" + str(starting_trial) + ".bin")
fileMatureIndicatorsBefore = os.path.join(dirname_before, "mature_" + str(starting_trial) + ".bin")
fileMaturationPropertiesBefore = os.path.join(dirname_before, "maturation_properties_" + str(starting_trial) + ".bin")

fileCoordAfter = os.path.join(dirname_after, "RA_xy_" + str(starting_trial) + "afterReading.bin")
fileActiveAfter = os.path.join(dirname_after, "RA_RA_active_connections_" + str(starting_trial) + "afterReading.bin")
fileSuperAfter = os.path.join(dirname_after, "RA_RA_super_connections_" + str(starting_trial) + "afterReading.bin")
fileRA2IAfter = os.path.join(dirname_after, "RA_I_connections_" + str(starting_trial) + "afterReading.bin")
fileI2RAAfter = os.path.join(dirname_after, "I_RA_connections_" + str(starting_trial) + "afterReading.bin")
fileWeightsAfter = os.path.join(dirname_after, "weights_" + str(starting_trial) + "afterReading.bin")
fileAxonalDelaysRA2IAfter = os.path.join(dirname_after, "axonal_delays_RA2I_" + str(starting_trial) + "afterReading.bin")
fileAxonalDelaysRA2RAAfter = os.path.join(dirname_after, "axonal_delays_RA2RA_" + str(starting_trial) + "afterReading.bin")
fileAxonalDelaysI2RAAfter = os.path.join(dirname_after, "axonal_delays_I2RA_" + str(starting_trial) + "afterReading.bin")
fileActivityHistoryAfter = os.path.join(dirname_after, "activity_history_" + str(starting_trial) + "afterReading.bin")
fileReplacementHistoryAfter = os.path.join(dirname_after, "replacement_history_" + str(starting_trial) + "afterReading.bin")
fileRemodeledIndicatorsAfter = os.path.join(dirname_after, "remodeled_indicators_" + str(starting_trial) + "afterReading.bin")
fileMatureIndicatorsAfter = os.path.join(dirname_after, "mature_" + str(starting_trial) + "afterReading.bin")
fileMaturationPropertiesAfter = os.path.join(dirname_before, "maturation_properties_" + str(starting_trial) + "afterReading.bin")

###########################
#### Read data from files
###########################

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

print np.array_equal(coordBefore, coordAfter)

print np.array_equal(np.array(active_synapses_before), np.array(active_synapses_after))
print np.array_equal(np.array(super_synapses_before), np.array(super_synapses_after))

print np.array_equal(np.array(targets_id_RA2I_before), np.array(targets_id_RA2I_after))
print np.array_equal(np.array(weights_RA2I_before), np.array(weights_RA2I_after))
print np.array_equal(np.array(syn_lengths_RA2I_before), np.array(syn_lengths_RA2I_after))
print np.array_equal(np.array(axonal_delays_RA2I_before), np.array(axonal_delays_RA2I_after))

print np.array_equal(np.array(targets_id_I2RA_before), np.array(targets_id_I2RA_after))
print np.array_equal(np.array(weights_I2RA_before), np.array(weights_I2RA_after))
print np.array_equal(np.array(syn_lengths_I2RA_before), np.array(syn_lengths_I2RA_after))
print np.array_equal(np.array(axonal_delays_I2RA_before), np.array(axonal_delays_I2RA_after))

print np.array_equal(weights_before, weights_after)
print np.array_equal(np.array(axonal_delays_I2RA_from_delay_file_before), np.array(axonal_delays_I2RA_from_delay_file_after))
print np.array_equal(np.array(axonal_delays_RA2RA_from_delay_file_before), np.array(axonal_delays_RA2RA_from_delay_file_after))
print np.array_equal(np.array(axonal_delays_RA2I_from_delay_file_before), np.array(axonal_delays_RA2I_from_delay_file_after))

print np.array_equal(activity_history_before, activity_history_after)
print np.array_equal(replacement_history_before, replacement_history_after)
print np.array_equal(remodeled_indicators_before, remodeled_indicators_after)

print np.array_equal(mature_indicators_before, mature_indicators_after)
print np.array_equal(maturation_rate_before, maturation_rate_after)
print np.array_equal(Erest_before, Erest_after)
print np.array_equal(GCa_before, GCa_after)
