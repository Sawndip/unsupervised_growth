#!/bin/sh
# path to network directory
# starting trial
# number of testing trials
# path to output directory
mpirun -np 23 ./main /home/eugene/Output/networks/chainGrowth/passiveDendrite/noImmatureOut4/ 23500 10 /home/eugene/Output/networks/chainGrowth/passiveDendrite/test/noImmatureOut4/
