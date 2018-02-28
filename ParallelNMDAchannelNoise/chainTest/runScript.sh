#!/bin/sh
# path to network directory
# starting trial
# number of testing trials
# path to output directory
mpirun -np 23 ./main /home/eugene/Output/networks/chainGrowth/passiveDendrite/events1/ 4300 50 /home/eugene/Output/networks/chainGrowth/passiveDendrite/test1/
