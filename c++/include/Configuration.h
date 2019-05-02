#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <libconfig.h++>

using namespace libconfig;

class Configuration
{
public:
    virtual void read_configuration(const char* filename) = 0; // read configuration from the file 
    virtual void write_configuration(const char* filename) const = 0; // write configuration to the file
    virtual void print_configuration() const = 0; // print configuration on the screeen
    
};

#endif
