#!/bin/sh
# "path to configuration file"
# path to data directory; path to output directory; starting trial; fraction of chain neurons to kill
mpirun -np 23 ./recovery /mnt/hodgkin_home/eugene/Output/networks/networkTest/ /mnt/hodgkin_home/eugene/Output/networks/networkTest/ 500 0.1
