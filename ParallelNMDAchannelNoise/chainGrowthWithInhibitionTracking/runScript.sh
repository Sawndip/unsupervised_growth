#!/bin/sh
# path to network directory
# path to training file
# path to growth configuration file
# path to output directory
# number of simulation trials
# time resolution of inhibitory conductance
mpirun -np 8 ./main /mnt/hodgkin/eugene/Output/networks/chainGrowth/network2000RA550I_v2/ /mnt/hodgkin/eugene/Output/networks/chainGrowth/network2000RA550I_v2/training_neurons_10.bin /mnt/hodgkin/eugene/Output/networks/chainGrowth/matTrans84/growth_parameters.cfg /mnt/hodgkin/eugene/Output/networks/chainGrowth/matTrans84/ 10 0.25
