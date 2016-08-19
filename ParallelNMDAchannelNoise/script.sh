#!/bin/sh
# run program with following command line arguments: Gie_mean, BETA, BETA_SUPERSYNAPSE, T_P, T_D, TAU_P, TAU_D,  A_P, A_D, A_P_SUPER, A_D_SUPER, ACTIVATION, 
#                                                   SUPERSYNAPSE_THRESHOLD, G_MAX, N_RA, num_inh_clusters_in_row, num_inh_in_cluster, N_ss, N_TR, outputFolder,
#													reading, filenumber, testing
mpirun -np 5 ./out 0.5 0.995 0.992 10 10 20 20 0.0015 0.00010 0.0040 0.00025 0.003 0.01 0.080 300 10 1 2 2 /home/eugene/Output/ 0 0 0
