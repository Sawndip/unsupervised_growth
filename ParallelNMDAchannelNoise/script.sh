#!/bin/sh
# run program with following command line arguments: sigma_soma, sigma_dend, Gie_mean, BETA, BETA_SUPERSYNAPSE, T_P, T_D, TAU_P, TAU_D,  
#													A_P, A_D, A_P_SUPER, A_D_SUPER, ACTIVATION, 
#                                                   SUPERSYNAPSE_THRESHOLD, G_MAX, N_RA, num_inh_clusters_in_row, num_inh_in_cluster, N_ss, N_TR, outputFolder,
#													reading, filenumber, testingm network_update_frequency
mpirun -np 6 ./out 100 150 0.5 0.9985 0.992 10 10 30 30 0.0015 0.00010 0.0040 0.00025 0.003 0.01 0.080 100 6 1 2 2 /home/eugene/Output/ 0 1 0 1.0
