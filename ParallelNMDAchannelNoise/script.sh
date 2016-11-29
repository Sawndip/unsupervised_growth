#!/bin/sh
# run program with following command line arguments: 
# mu_soma, sigma_soma, mu_dend, sigma_dend,
# Gie_mean, Ei, BETA, BETA_SUPERSYNAPSE, T_P, T_D, TAU_P, TAU_D,  A_P, A_D, A_P_SUPER, A_D_SUPER, f0, ACTIVATION, 
# MATURATION_THRESHOLD, SUPERSYNAPSE_THRESHOLD, G_MAX, N_RA, num_inh_clusters_in_row, num_inh_in_cluster, N_ss, N_TR, outputFolder,
# reading, filenumber, testing, training,  network_update_frequency, A_RA2I, SIGMA_RA2I, B_I2RA, SIGMA_I2RA
mpirun -np 23 ./out 0 0 0 170 0.35 -55 0.997 0.975 10 10 30 30 0.00005 0.000002 0.0045 0.00030 1.2 0.0001 0.0015 0.160 0.045 300 10 1 4 4 /home/eugene/Output/ 0 1 0 1 1.0 5.0 8.0 5.0 6.0
