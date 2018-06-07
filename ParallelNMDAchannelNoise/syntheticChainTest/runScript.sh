#!/bin/sh
# number of layers in synfire chain
# number of neurons in each layer
# probability to connect neuron in the next group
# strength of connection between neurons in the chain
# path to network directory
# path to growth configuration file
# number of testing trials
# path to output directory
mpirun -np 23 ./main 50 10 0.5 800.0 /home/eugene/Output/networks/chainGrowth/network2000RA550I/ /home/eugene/Output/networks/chainGrowth/passiveDendrite/chain_0.5_delay_3.0/growth_parameters.cfg 20 /home/eugene/Output/networks/chainGrowth/passiveDendrite/chain_0.5_delay_3.0/
