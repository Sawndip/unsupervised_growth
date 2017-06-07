#ifndef CONFIGURATION_NETWORK_GENERATOR_H
#define CONFIGURATION_NETWORK_GENERATOR_H

#include <libconfig.h++>
#include <string>
#include "Configuration.h"

using namespace libconfig;

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
};

class ConfigurationNetworkGenerator : public Configuration
{
public:
    virtual void read_configuration(const char* filename); // read configuration from the file 
    virtual void write_configuration(const char* filename) const; // write configuration to the file
    virtual void print_configuration() const; // print configuration on the screeen
    
    struct SpatialParameters get_spatial_parameters() const {return spatial_params;};

private:
    struct SpatialParameters spatial_params;
};

#endif
