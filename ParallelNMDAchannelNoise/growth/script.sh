#!/bin/sh
# path to network directory
# path to training file
# path to growth configuration file
# path to output directory
mpirun -np 23 ./growth /home/eugene/Output/networks/sphere/ /home/eugene/Output/networks/sphere/training/clustered_training.bin /home/eugene/Output/networks/sphereTest/ /home/eugene/Output/networks/sphereTest/
