#!/bin/sh
# path to network directory
# starting trial
# number of testing trials
# path to output directory
#mpirun -np 8 ./main /mnt/hodgkin/eugene/results/immature/clusters/matTrans78/ 86400 10 /mnt/hodgkin/eugene/results/immature/clusters/test/matTrans78/trial86400/
mpirun -np 8 ./main /mnt/hodgkin/eugene/Output/networks/chainGrowth/matTrans85/ 2000 10 /mnt/hodgkin/eugene/Output/networks/chainGrowth/matTrans85/test/
