#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <libconfig.h++>
#include <string>

using namespace libconfig;

// structure with network parameters
struct NetworkParameters
{
    int N_RA; // number of HVC(RA) neurons
    int N_TR; // number of training HVC(RA) neurons
    int N_I; // number of HVC(I) neurons
};

// structure with synaptic parameters
struct SynapticParameters
{
    int Nss; // maximum number of allowed outgoing supersynapses per neuron

    double R; // learning rate (default value is 1)

    double A_P; // constant for LTP
    double A_D; // constant for LTD
    double F_0; // punishment for synchronous neuronal firing
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

// structure with spatial parameters
struct SpatialParameters
{
    double SIDE; // side of the square modeling HVC
    double MIN_INTERNEURON_DISTANCE; // minimum allowed distance between two neurons

    // spatial connectivity distributions (modeled by Gaussians)
    double A_RA2I; // constant for HVC(RA) -> HVC(I) connections
    double SIGMA_RA2I; // standard deviation for HVC(RA) -> HVC(I) connections
    double B_I2RA; // constant for HVC(I) -> HVC(RA) connections
    double SIGMA_I2RA; // standard deviation for HVC(I) -> HVC(RA) connections

    // distributions for connection strength
    double Gei_mean; // mean HVC(RA) -> HVC(I) connection strength
    double Gei_std; // standard deviation of HVC(RA) -> HVC(I) connection strength
    
    double Gie_mean; // mean HVC(I) -> HVC(RA) connection strength
    double Gie_std; // standard deviation of HVC(I) -> HVC(RA) connection strength
};

// structure with maturation parameters
struct GabaParameters
{
    double E_GABA_MATURE; // GABA reverse potential for mature neurons
    double E_GABA_IMMATURE; // GABA reverse potential for newborn neurons

    double GABA_RATE_THRESHOLD; // burst rate threshold for maturation. if exceeded, GABA reverse potential decreases by GABA_DOWN
    double MATURATION_RATE_THRESHOLD; // burst rate threshold for neuronal maturation
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


class Configuration
{
public:
    void read_configuration(const char* filename); // read configuration from the file 
    void write_configuration(const char* filename) const; // write configuration to the file
    void print_configuration(); // print configuration on the screeen
    

    struct NetworkParameters get_network_parameters() const {return network_params;};
    struct SynapticParameters get_synaptic_parameters() const {return synaptic_params;};
    struct SpatialParameters get_spatial_parameters() const {return spatial_params;};
    struct GabaParameters get_gaba_parameters() const {return gaba_params;};
    struct TimeParameters get_time_parameters() const {return time_params;};
    struct NoiseParameters get_noise_parameters() const {return noise_params;};

    std::string get_output_directory() const {return outputDirectory;};
private:
    struct NetworkParameters network_params;
    struct SynapticParameters synaptic_params;
    struct SpatialParameters spatial_params;
    struct GabaParameters gaba_params;
    struct TimeParameters time_params;
    struct NoiseParameters noise_params;

    std::string outputDirectory; // directory for data output
};

#endif
