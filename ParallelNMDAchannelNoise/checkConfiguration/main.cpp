#include "ConfigurationNetworkGrowth.h"
#include "ConfigurationNetworkTopology.h"
#include <string>
#include <iostream>

int main(int argc, char **argv)
{
    std::string configGrowthFilename;
    std::string configTopologyFilename;

    if (argc > 2)
    {
        configGrowthFilename = argv[1];
        configTopologyFilename = argv[2];
    }
    else 
    {
        std::cerr << "Usage: <configuration network growth file name> <configuration network topology file name>" << std::endl;
        return -1;
    }

    ConfigurationNetworkGrowth configGrowth;
    ConfigurationNetworkTopology configTopology;

    configGrowth.read_configuration(configGrowthFilename.c_str());
    configGrowth.print_configuration();

    configTopology.read_configuration(configTopologyFilename.c_str());
    configTopology.print_configuration();
    
    return 0;
}
