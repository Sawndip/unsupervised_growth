#ifndef CONFIGURATION_NETWORK_TOPOLOGY_H
#define CONFIGURATION_NETWORK_TOPOLOGY_H

#include <libconfig.h++>
#include <string>
#include "Configuration.h"

using namespace libconfig;

// structure with spatial parameters
struct TopologyParameters
{
	std::string arrangement; // arrangement of neurons; can be square; sphere; cube

    // spatial connectivity distributions (modeled by Gaussians)
    double A_RA2I; // constant for HVC(RA) -> HVC(I) connections
    double SIGMA_RA2I; // standard deviation for HVC(RA) -> HVC(I) connections
    double B_I2RA; // constant for HVC(I) -> HVC(RA) connections
    double SIGMA_I2RA; // standard deviation for HVC(I) -> HVC(RA) connections

    TopologyParameters() : A_RA2I(0.0), SIGMA_RA2I(0.0), B_I2RA(0.0), SIGMA_I2RA(0.0) {};
};

class ConfigurationNetworkTopology : public Configuration
{
public:
    virtual void read_configuration(const char* filename); // read configuration from the file 
    virtual void write_configuration(const char* filename) const; // write configuration to the file
    virtual void print_configuration() const; // print configuration on the screeen
    
    struct TopologyParameters get_topology_parameters() const {return topology_params;};
	void set_topology_parameters(const struct TopologyParameters &params) {topology_params = params;};
private:
    struct TopologyParameters topology_params;
};

#endif
