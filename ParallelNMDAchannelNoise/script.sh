#!/bin/sh
# run program with following command line arguments:
# N_TR, N_RA, N_I, N_ss,
# Gie_mean, A_RA2I, SIGMA_RA2I, B_I2RA, SIGMA_I2RA,
# mu_soma, sigma_soma, mu_dend, sigma_dend,
# A_P, A_D, f0, BETA, BETA_SUPERSYNAPSE, ACTIVATION, 
# SUPERSYNAPSE_THRESHOLD, G_MAX,
# GABA_RATE_THRESHOLD, MATURATION_RATE_THRESHOLD, DEATH_RATE_THRESHOLD,
# RATE_WINDOW_SHORT, RATE_WINDOW_LONG, GABA_DOWN
# outputDirectory
mpirun -np 20 ./out 4 300 100 4 0.20 5.0 8.0 5.0 6.0 0 100 0 200 0.000025 0.000025 1.2 0.9999 0.9999 0.0005 0.015 0.050 0.05 0.9 0.01 50 500 0.065 /home/eugene/Output/
