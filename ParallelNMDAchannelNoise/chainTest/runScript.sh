#!/bin/sh
# path to network directory
# starting trial
# number of testing trials
# path to output directory
mpirun -np 8 ./main /mnt/hodgkin/eugene/results/immature/clusters/matTrans63/ 25200 10 /mnt/hodgkin/eugene/results/immature/clusters/test/matTrans63/trial25200/
