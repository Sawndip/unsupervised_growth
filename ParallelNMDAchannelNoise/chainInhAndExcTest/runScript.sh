#!/bin/sh
# path to network directory
# file with training neurons
# path to output directory
mpirun -np 23 ./main /home/eugene/Output/networks/chainGrowth/test/ /home/eugene/Output/networks/chainGrowth/test/training_neurons_random.bin /home/eugene/Output/networks/chainGrowth/passiveDendrite/inhAndExcTest/
