#!/bin/sh
# run program with following command line arguments: 
# mu_soma, sigma_soma, mu_dend, sigma_dend,
# N_RA, num_inh_clusters_in_row, num_inh_in_cluster, N_ss, N_TR, num_trials, outputFolder, filenumber
mpirun -np 23 ./chainTest 0 0 0 150 300 10 1 4 4 100 /home/eugene/Output/ 163
