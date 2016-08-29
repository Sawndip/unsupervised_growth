#!/bin/sh
# run program with following command line arguments: 
# sigma_soma, sigma_dend,
# Gie_mean, BETA, BETA_SUPERSYNAPSE, T_P, T_D, TAU_P, TAU_D,  A_P, A_D, A_P_SUPER, A_D_SUPER, ACTIVATION, 
# SUPERSYNAPSE_THRESHOLD, G_MAX, N_RA, num_inh_clusters_in_row, num_inh_in_cluster, N_ss, N_TR, outputFolder,
# reading, filenumber, testing, network_update_frequency, A_RA2I, B_RA2I, LAMBDA_RA2I, MEAN_RA2I, SIGMA_RA2I, C_I2RA, LAMBDA_I2RA 
mpirun -np 23 ./out 100 150 0.5 0.9980 0.996 10 10 20 20 0.0015 0.00010 0.0015 0.00010 0.003 0.01 0.060 300 10 1 4 4 /home/eugene/Output/ 0 1 0 1.0 4.0 0.025 2.0 30.0 20.0 4.0 2.0
