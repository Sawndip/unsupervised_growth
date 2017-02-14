#!/bin/sh
# run program with following command line arguments:
# N_RA, N_I
# network_update_frequency
# outputFile, outputDirectory
mpirun -np 24 ./out 300 100 1.0 exec_time_1.bin /home/eugene/Output/
