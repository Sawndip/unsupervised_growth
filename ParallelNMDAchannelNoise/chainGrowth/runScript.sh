#!/bin/sh
# path to network directory
# path to training file
# path to growth configuration file
# path to output directory
mpirun -np 23 ./main /home/eugene/Output/networks/chainGrowth/network2000/ /home/eugene/Output/networks/chainGrowth/network2000/training_neurons_20.bin /home/eugene/Output/networks/chainGrowth/passiveDendrite/matTrans3_network2000/growth_parameters.cfg /home/eugene/Output/networks/chainGrowth/passiveDendrite/matTrans3_network2000/
