#!/bin/sh
# run program with following command line arguments: Gie_mean, BETA, BETA_SUPERSYNAPSE, A_P, A_D, ACTIVATION, 
#                                                   SUPERSYNAPSE_THRESHOLD, G_MAX, N_RA, num_inh_clusters_in_row, num_inh_in_cluster, N_ss, N_TR, outputFolder
mpirun -np 5 ./out 0.5 0.995 0.999 0.0015 0.00010 0.003 0.01 0.070 60 4 2 2 2 Output/
