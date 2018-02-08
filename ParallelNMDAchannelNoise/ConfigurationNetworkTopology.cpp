#include "ConfigurationNetworkTopology.h"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

void ConfigurationNetworkTopology::read_configuration(const char* filename)
{
    
    Config cfg;

    // Read the file. If there is an error, report it and exit.
    try
    {
        cfg.readFile(filename);
    }
    catch(const FileIOException &fioex)
    {
        std::cerr << "I/O error while reading file." << std::endl;
    }
    catch(const ParseException &pex)
    {
        std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                  << " - " << pex.getError() << std::endl;
    }
    

    const Setting& root = cfg.getRoot();
    
    // Get the spatial parameters
    try
    {
        const Setting& topology_params_setting = root["spatialParameters"];

        // Only output the record if all of the expected fields are present.
        topology_params_setting.lookupValue("arrangement", topology_params.arrangement);
          
        topology_params_setting.lookupValue("A_RA2I", topology_params.A_RA2I);
        topology_params_setting.lookupValue("SIGMA_RA2I", topology_params.SIGMA_RA2I);
        
        topology_params_setting.lookupValue("B_I2RA", topology_params.B_I2RA);
        topology_params_setting.lookupValue("SIGMA_I2RA", topology_params.SIGMA_I2RA);
    }
    catch(const SettingNotFoundException &nfex)
    {
        std::cerr << "Topology parameter setting is not found in configuration file!" << std::endl;
    }

}

void ConfigurationNetworkTopology::print_configuration() const
{
    // display all spatial parameters read from file
    std::cout << std::endl << "Topology parameters read from configuration file: " << std::endl << std::endl;
    
    std::cout << "arrangement = " << topology_params.arrangement << std::endl;
 
    std::cout << "A_RA2I = " << topology_params.A_RA2I << std::endl;
    std::cout << "SIGMA_RA2I = " << topology_params.SIGMA_RA2I << std::endl;
    std::cout << "B_I2RA = " << topology_params.B_I2RA << std::endl;
    std::cout << "SIGMA_I2RA = " << topology_params.SIGMA_I2RA << std::endl << std::endl;
}

void ConfigurationNetworkTopology::write_configuration(const char* filename) const
{
    std::ofstream out;

    out.open(filename, std::ios::out);

    // write spatial parameters
    out << "spatialParameters = " << std::endl;
    out << "{" << std::endl;
    out << "\tarrangement = \"" << topology_params.arrangement << "\";" << std::endl;
   
    out << "\tA_RA2I = " << std::fixed << std::setprecision(6) << topology_params.A_RA2I << ";" << std::endl;
    out << "\tSIGMA_RA2I = " << topology_params.SIGMA_RA2I << ";" << std::endl;
    out << "\tB_I2RA = " << topology_params.B_I2RA << ";" << std::endl;
    out << "\tSIGMA_I2RA = " << topology_params.SIGMA_I2RA << ";" << std::endl << std::endl;
    out << "};" << std::endl << std::endl;

    out.close();
}
