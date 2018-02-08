#ifndef CONFIGURATION_NETWORK_GROWTH_H
#define CONFIGURATION_NETWORK_GROWTH_H

#include <libconfig.h++>
#include <string>
#include "Configuration.h"

using namespace libconfig;

// structure with synaptic parameters
struct SynapticParameters
{
    int Nss; // maximum number of allowed outgoing supersynapses per neuron

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
    
    SynapticParameters() : Nss(0), A_P(0.0), A_D(0.0), T_0(0.0), T_P(0.0),
						   T_D(0.0), TAU_P(0.0), TAU_D(0.0), BETA(0.0), BETA_SUPERSYNAPSE(0.0),
						   ACTIVATION_THRESHOLD(0.0), SUPERSYNAPSE_THRESHOLD(0.0), WEIGHT_MAX(0.0)
						   {};
};

// structure with maturation parameters
struct MaturationParameters
{
	double E_REST_MATURE; // resting potential for mature neurons
    double E_REST_IMMATURE; // resting potential for newborn neurons

    double E_GABA_MATURE; // GABA reverse potential for mature neurons
    double E_GABA_IMMATURE; // GABA reverse potential for newborn neurons

	double AD_MATURE; // dendritic area for mature neurons
    double AD_IMMATURE; // dendritic area for newborn neurons

	double GK_MATURE; // potassium somatic conductance for mature neurons
    double GK_IMMATURE; // potassium somatic conductance for newborn neurons
    
    double GNA_MATURE; // sodium somatic conductance for mature neurons
    double GNA_IMMATURE; // sodium somatic conductance for newborn neurons

    double RC_MATURE; // coupling resistance for mature neurons
    double RC_IMMATURE; // coupling resistance for newborn neurons

    double DEATH_RATE_THRESHOLD; // burst rate threshold for neuronal replacement

    int RATE_WINDOW_LONG; // time window in which "average" firing rate for a long period is measured
    
    MaturationParameters() : E_REST_MATURE(-80.0), E_REST_IMMATURE(-65.0), E_GABA_MATURE(-80.0), E_GABA_IMMATURE(-55.0),
							 AD_MATURE(10000.0), AD_IMMATURE(1000.0), GK_MATURE(8.0), GK_IMMATURE(16.0), 
							 GNA_MATURE(60.0), GNA_IMMATURE(40.0), RC_MATURE(55.0), RC_IMMATURE(1.0), 
							 DEATH_RATE_THRESHOLD(0.01), RATE_WINDOW_LONG(2000)
							 {};
};

// structure with noise parameters
struct NoiseParameters
{
    double white_noise_mean_soma; // mean of white noise current injected to soma
    double white_noise_std_soma; // standard deviation of white noise current injected to soma
    
    double white_noise_mean_dend; // mean of white noise current injected to dendrite
    double white_noise_std_dend; // standard deviation of white noise cuurent injected to dendrite
    
    NoiseParameters() : white_noise_mean_soma(0.0), white_noise_std_soma(0.01),
						white_noise_mean_dend(0.0), white_noise_std_dend(0.025)
						{};
};

// structure with network parameters
struct ConnectionParameters
{
    double Gei_max; // max HVC(RA) -> HVC(I)  connection strength
	double Gie_max; // max HVC(I)  -> HVC(RA) connection strength
    
	double delay_constant; // constant which converts distances between neurons into axonal delays

	ConnectionParameters() : Gei_max(0.25), Gie_max(0.0), delay_constant(0.0) {};
};


class ConfigurationNetworkGrowth : public Configuration
{
public:
    void read_configuration(const char* filename); // read configuration from the file 
    void write_configuration(const char* filename) const; // write configuration to the file
    void print_configuration() const; // print configuration on the screeen
    
    struct SynapticParameters get_synaptic_parameters() const {return synaptic_params;};
    struct MaturationParameters get_maturation_parameters() const {return maturation_params;};
    struct NoiseParameters get_noise_parameters() const {return noise_params;};
    struct ConnectionParameters get_connection_parameters() const {return connection_params;};
    
private:
    struct SynapticParameters synaptic_params;
    struct MaturationParameters maturation_params;
    struct NoiseParameters noise_params;
    struct ConnectionParameters connection_params;
    
};

#endif
