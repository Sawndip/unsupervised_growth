#ifndef STRUCTURES_H
#define STRUCTURES_H

struct SynapticParameters
{
    double R;

    double A_P;
    double A_D;
    double F_0;
    double T_P;
    double T_D;
    double TAU_P;
    double TAU_D;
    
    double BETA;
    double BETA_SUPERSYNAPSE;

    double ACTIVATION;
    double SUPERSYNAPSE_THRESHOLD;
    double G_MAX;

    double STDP_WINDOW;
};

struct SpatialConnectivityParameters
{
    double SIDE;
    double MIN_INTERNEURON_DISTANCE;

    double A_RA2I;
    double SIGMA_RA2I;
    double B_I2RA;
    double SIGMA_I2RA;

    double Gei_mean;
    double Gei_std;
    
    double Gie_mean;
    double Gie_std;
};

struct GabaParameters
{
    double E_GABA_MATURE;
    double E_GABA_IMMATURE;

    double GABA_RATE_THRESHOLD;
    double MATURATION_RATE_THRESHOLD;
    double DEATH_RATE_THRESHOLD;

    int RATE_WINDOW_SHORT;
    int RATE_WINDOW_LONG;

    double GABA_DOWN;
};

struct TimeParameters
{
    double timestep;
    double network_update_frequency;
    double trial_duration;
    double WAITING_TIME;
};

struct NoiseParameters
{
    double mean_soma;
    double std_soma;
    
    double mean_dend;
    double std_dend;
};

#endif
