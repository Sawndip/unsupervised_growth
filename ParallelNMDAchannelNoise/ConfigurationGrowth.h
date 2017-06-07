#ifndef CONFIGURATION_GROWTH_H
#define CONFIGURATION_GROWTH_H

#include <libconfig.h++>
#include <string>
#include "Configuration.h"

using namespace libconfig;

// structure with synaptic parameters
struct SynapticParameters
{
    int Nss; // maximum number of allowed outgoing supersynapses per neuron

    double R; // learning rate (default value is 1)

    double A_P; // constant for LTP
    double A_D; // constant for LTD
    double T_0; // shift of STDP rule to the right
    double T_P; // best time for LTP response
    double T_D; // best time for LTD response
    double TAU_P; // time scale of LTP efficiency
    double TAU_D; // time scale of LTD efficiency
    
    double BETA; // potentiation decay constant for active and silent synapses
    double BETA_SUPERSYNAPSE; // potentiation decay constant for super synapses

    double ACTIVATION_THRESHOLD; // threshold for synapse activation
    double SUPERSYNAPSE_THRESHOLD; // threshold for supersynapse activation
    double WEIGHT_MAX; // maximum synaptic weight

    double STDP_WINDOW; // time window size in which neuronal spikes still matter for STDP rules
    
};

// structure with maturation parameters
struct GabaParameters
{
    double E_GABA_MATURE; // GABA reverse potential for mature neurons
    double E_GABA_IMMATURE; // GABA reverse potential for newborn neurons

    double GABA_RATE_THRESHOLD; // burst rate threshold for maturation. if exceeded, GABA reverse potential decreases by GABA_DOWN
    double DEATH_RATE_THRESHOLD; // burst rate threshold for neuronal replacement

    int RATE_WINDOW_SHORT; // time window in which "instantaneous" firing rate for a short period is measured
    int RATE_WINDOW_LONG; // time window in which "average" firing rate for a long period is measured

    double GABA_DOWN; // decrease in GABA potential
};

// structure with paramters of neuronal dynamics and simulation trials
struct TimeParameters
{
    double timestep; // timestep for solving neuronal dynamics in ms
    double network_update_frequency; // frequency of communication between processes in ms
    double trial_duration; // trial duration in s
    double WAITING_TIME; // time in the beginning of trial before which training current cannot be injected
};

// structure with noise parameters
struct NoiseParameters
{
    double white_noise_mean_soma; // mean of white noise current injected to soma
    double white_noise_std_soma; // standard deviation of white noise current injected to soma
    
    double white_noise_mean_dend; // mean of white noise current injected to dendrite
    double white_noise_std_dend; // standard deviation of white noise cuurent injected to dendrite
};

// structure with inhibitory strength parameters
struct InhibitoryParameters
{
	double Gie_mean; // mean HVC(I) -> HVC(RA) connection strength
    double Gie_std; // standard deviation of HVC(I) -> HVC(RA) connection strength
};


class ConfigurationGrowth : public Configuration
{
public:
    void read_configuration(const char* filename); // read configuration from the file 
    void write_configuration(const char* filename) const; // write configuration to the file
    void print_configuration() const; // print configuration on the screeen
    
    struct SynapticParameters get_synaptic_parameters() const {return synaptic_params;};
    struct GabaParameters get_gaba_parameters() const {return gaba_params;};
    struct TimeParameters get_time_parameters() const {return time_params;};
    struct NoiseParameters get_noise_parameters() const {return noise_params;};
    struct InhibitoryParameters get_inhibitory_parameters() const {return inhibitory_params;};
    
private:
    struct SynapticParameters synaptic_params;
    struct GabaParameters gaba_params;
    struct TimeParameters time_params;
    struct NoiseParameters noise_params;
    struct InhibitoryParameters inhibitory_params;
    
};

#endif
