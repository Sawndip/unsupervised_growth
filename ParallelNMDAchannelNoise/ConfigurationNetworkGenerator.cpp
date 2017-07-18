#include "ConfigurationNetworkGenerator.h"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

void ConfigurationNetworkGenerator::read_configuration(const char* filename)
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
        const Setting& spatial_params_setting = root["spatialParameters"];

        // Only output the record if all of the expected fields are present.
        spatial_params_setting.lookupValue("dimensionality", spatial_params.dimensionality);
          
        spatial_params_setting.lookupValue("MIN_INTERNEURON_DISTANCE", spatial_params.MIN_INTERNEURON_DISTANCE);
        
        
        spatial_params_setting.lookupValue("A_RA2I", spatial_params.A_RA2I);
        spatial_params_setting.lookupValue("SIGMA_RA2I", spatial_params.SIGMA_RA2I);
        
        spatial_params_setting.lookupValue("B_I2RA", spatial_params.B_I2RA);
        spatial_params_setting.lookupValue("SIGMA_I2RA", spatial_params.SIGMA_I2RA);
        
        spatial_params_setting.lookupValue("Gei_mean", spatial_params.Gei_mean);
        spatial_params_setting.lookupValue("Gei_std", spatial_params.Gei_std);
        
    }
    catch(const SettingNotFoundException &nfex)
    {
        std::cerr << "Spatial parameter setting is not found in configuration file!" << std::endl;
    }

}

void ConfigurationNetworkGenerator::print_configuration() const
{
    // display all spatial parameters read from file
    std::cout << std::endl << "Spatial parameters read from configuration file: " << std::endl << std::endl;
    
    std::cout << "dimensionality = " << spatial_params.dimensionality << std::endl;
    std::cout << "MIN_INTERNEURON_DISTANCE = " << spatial_params.MIN_INTERNEURON_DISTANCE << std::endl << std::endl;
    
    std::cout << "A_RA2I = " << spatial_params.A_RA2I << std::endl;
    std::cout << "SIGMA_RA2I = " << spatial_params.SIGMA_RA2I << std::endl;
    std::cout << "B_I2RA = " << spatial_params.B_I2RA << std::endl;
    std::cout << "SIGMA_I2RA = " << spatial_params.SIGMA_I2RA << std::endl << std::endl;

    std::cout << "Gei_mean = " << spatial_params.Gei_mean << std::endl;
    std::cout << "Gei_std = " << spatial_params.Gei_std << std::endl;
}

void ConfigurationNetworkGenerator::write_configuration(const char* filename) const
{
    std::ofstream out;

    out.open(filename, std::ios::out);

    // write spatial parameters
    out << "spatialParameters = " << std::endl;
    out << "{" << std::endl;
    out << "\tdimensionality = " << spatial_params.dimensionality << ";" << std::endl;
    out << "\tMIN_INTERNEURON_DISTANCE = " << std::fixed << std::setprecision(6) << spatial_params.MIN_INTERNEURON_DISTANCE << ";" << std::endl << std::endl;
    
    out << "\tA_RA2I = " << spatial_params.A_RA2I << ";" << std::endl;
    out << "\tSIGMA_RA2I = " << spatial_params.SIGMA_RA2I << ";" << std::endl;
    out << "\tB_I2RA = " << spatial_params.B_I2RA << ";" << std::endl;
    out << "\tSIGMA_I2RA = " << spatial_params.SIGMA_I2RA << ";" << std::endl << std::endl;
    
    out << "\tGei_mean = " << spatial_params.Gei_mean << ";" << std::endl;
    out << "\tGei_std = " << spatial_params.Gei_std << ";" << std::endl;
    out << "};" << std::endl << std::endl;

    out.close();
}
