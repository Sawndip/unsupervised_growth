#!/bin/sh
# path to network directory
# path to training file
# path to growth configuration file
# path to output directory
mpirun -np 23 ./main /home/eugene/Output/networks/chainGrowth/network2000RA550I/ /home/eugene/Output/networks/chainGrowth/network2000RA550I/training_neurons_20.bin /home/eugene/Output/networks/chainGrowth/passiveDendrite/matTrans4_network2000RA550I/growth_parameters.cfg /home/eugene/Output/networks/chainGrowth/passiveDendrite/matTrans4_network2000RA550I/
