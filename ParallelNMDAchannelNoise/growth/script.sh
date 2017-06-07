#!/bin/sh
# "path to configuration file"
# path to network directory
mpirun -np 8 ./growth /home/eugene/Output/networks/generateTest/ /home/eugene/Output/networks/networkTest/ /home/eugene/Output/networks/networkTest/
