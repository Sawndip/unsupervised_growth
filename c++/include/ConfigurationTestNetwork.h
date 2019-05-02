#ifndef CONFIGURATION_TEST_NETWORK_H
#define CONFIGURATION_TEST_NETWORK_H

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

class ConfigurationTestNetwork : public Configuration
{
public:
    void read_configuration(const char* filename); // read configuration from the file 
    void write_configuration(const char* filename) const; // write configuration to the file
    void print_configuration() const; // print configuration on the screeen
    
    struct SynapticParameters get_synaptic_parameters() const {return synaptic_params;};
    struct GabaParameters get_gaba_parameters() const {return gaba_params;};
    
private:
    struct SynapticParameters synaptic_params;
    struct GabaParameters gaba_params;
};

#endif
