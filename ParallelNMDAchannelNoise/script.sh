#!/bin/sh
# run program with following command line arguments: 
# mu_soma, sigma_soma, mu_dend, sigma_dend,
# Gie_mean, Ei, BETA, BETA_SUPERSYNAPSE, T_P, T_D, TAU_P, TAU_D,  A_P, A_D, A_P_SUPER, A_D_SUPER, f0, ACTIVATION, 
# SUPERSYNAPSE_THRESHOLD, G_MAX, GABA_DOWN, N_RA, N_I, N_ss, N_TR, outputFolder,
# reading, filenumber, training,  network_update_frequency, A_RA2I, SIGMA_RA2I, B_I2RA, SIGMA_I2RA
mpirun -np 23 ./out 0 0 0 170 0.35 -55 0.995 0.9995 10 10 30 30 0.00010 0.000010 0.00010 0.000010 1.2 0.0001 0.0015 0.030 0.025 300 100 4 4 /home/eugene/Output/ 0 1 1 1.0 5.0 8.0 5.0 6.0
