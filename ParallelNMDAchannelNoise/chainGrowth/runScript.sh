#!/bin/sh
# path to network directory
# path to training file
# path to growth configuration file
# path to output directory
mpirun -np 23 ./main /home/eugene/Output/networks/chainGrowth/test/ /home/eugene/Output/networks/chainGrowth/test/training_neurons_random.bin /home/eugene/Output/networks/chainGrowth/testGrowthDelays/growth_parameters.cfg /home/eugene/Output/networks/chainGrowth/testGrowthDelays/
