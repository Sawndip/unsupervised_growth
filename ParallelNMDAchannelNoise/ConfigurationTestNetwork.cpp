#include "ConfigurationTestNetwork.h"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

void ConfigurationTestNetwork::read_configuration(const char* filename)
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
    
    // Get the synaptic parameters
    try
    {
        const Setting& synaptic_params_setting = root["synapticParameters"];

        synaptic_params_setting.lookupValue("Nss", synaptic_params.Nss);
          
        synaptic_params_setting.lookupValue("A_P", synaptic_params.A_P);
        synaptic_params_setting.lookupValue("A_D", synaptic_params.A_D);
        synaptic_params_setting.lookupValue("T_0", synaptic_params.T_0);
        synaptic_params_setting.lookupValue("T_P", synaptic_params.T_P);
        synaptic_params_setting.lookupValue("T_D", synaptic_params.T_D);
        synaptic_params_setting.lookupValue("TAU_P", synaptic_params.TAU_P);
        synaptic_params_setting.lookupValue("TAU_D", synaptic_params.TAU_D);
          
        synaptic_params_setting.lookupValue("BETA", synaptic_params.BETA);
        synaptic_params_setting.lookupValue("BETA_SUPERSYNAPSE", synaptic_params.BETA_SUPERSYNAPSE);
          
        synaptic_params_setting.lookupValue("ACTIVATION_THRESHOLD", synaptic_params.ACTIVATION_THRESHOLD);
        synaptic_params_setting.lookupValue("SUPERSYNAPSE_THRESHOLD", synaptic_params.SUPERSYNAPSE_THRESHOLD);
        synaptic_params_setting.lookupValue("WEIGHT_MAX", synaptic_params.WEIGHT_MAX);
    }
    catch(const SettingNotFoundException &nfex)
    {
        std::cerr << "Synaptic parameter setting is not found in configuration file!" << std::endl;
    }

    // Get the maturation parameters
    try
    {
        const Setting& gaba_params_setting = root["gabaParameters"];

        gaba_params_setting.lookupValue("E_GABA_MATURE", gaba_params.E_GABA_MATURE);
        gaba_params_setting.lookupValue("E_GABA_IMMATURE", gaba_params.E_GABA_IMMATURE);
        
        gaba_params_setting.lookupValue("GABA_RATE_THRESHOLD", gaba_params.GABA_RATE_THRESHOLD);
        gaba_params_setting.lookupValue("DEATH_RATE_THRESHOLD", gaba_params.DEATH_RATE_THRESHOLD);
        
        gaba_params_setting.lookupValue("RATE_WINDOW_SHORT", gaba_params.RATE_WINDOW_SHORT);
        gaba_params_setting.lookupValue("RATE_WINDOW_LONG", gaba_params.RATE_WINDOW_LONG);
        
        gaba_params_setting.lookupValue("GABA_DOWN", gaba_params.GABA_DOWN);
        
    }
    catch(const SettingNotFoundException &nfex)
    {
        std::cerr << "Maturation parameter setting is not found in configuration file!" << std::endl;
    }
}

void ConfigurationTestNetwork::print_configuration() const
{
    // display all synaptic parameters read from file
    std::cout << std::endl << "Synaptic parameters read from configuration file: " << std::endl << std::endl;
 
    std::cout << "Nss = " << synaptic_params.Nss << std::endl << std::endl;
    
    std::cout << "A_P = " << synaptic_params.A_P << std::endl;
    std::cout << "A_D = " << synaptic_params.A_D << std::endl;
    std::cout << "T_0 = " << synaptic_params.T_0 << std::endl;
    std::cout << "T_P = " << synaptic_params.T_P << std::endl;
    std::cout << "T_D = " << synaptic_params.T_D << std::endl;
    std::cout << "TAU_P = " << synaptic_params.TAU_P << std::endl;
    std::cout << "TAU_D = " << synaptic_params.TAU_D << std::endl << std::endl;
    
    std::cout << "BETA = " << synaptic_params.BETA << std::endl;
    std::cout << "BETA_SUPERSYNAPSE = " << synaptic_params.BETA_SUPERSYNAPSE << std::endl << std::endl;
    
    std::cout << "ACTIVATION_THRESHOLD = " << synaptic_params.ACTIVATION_THRESHOLD << std::endl;
    std::cout << "SUPERSYNAPSE_THRESHOLD = " << synaptic_params.SUPERSYNAPSE_THRESHOLD << std::endl;
    std::cout << "WEIGHT_MAX = " << synaptic_params.WEIGHT_MAX << std::endl << std::endl;

    // display all maturation parameters read from file
    std::cout << std::endl << "Maturation parameters read from configuration file: " << std::endl << std::endl;
    
    std::cout << "E_GABA_MATURE = " << gaba_params.E_GABA_MATURE << std::endl;
    std::cout << "E_GABA_IMMATURE = " << gaba_params.E_GABA_IMMATURE << std::endl << std::endl;
    
    std::cout << "GABA_RATE_THRESHOLD = " << gaba_params.GABA_RATE_THRESHOLD << std::endl;
    std::cout << "DEATH_RATE_THRESHOLD = " << gaba_params.DEATH_RATE_THRESHOLD << std::endl << std::endl;

    std::cout << "RATE_WINDOW_SHORT = " << gaba_params.RATE_WINDOW_SHORT << std::endl;
    std::cout << "RATE_WINDOW_LONG = " << gaba_params.RATE_WINDOW_LONG << std::endl << std::endl;

    std::cout << "GABA_DOWN = " << gaba_params.GABA_DOWN << std::endl; 
}

void ConfigurationTestNetwork::write_configuration(const char* filename) const
{
    std::ofstream out;

    out.open(filename, std::ios::out);
    
    // write synaptic parameters
    out << "synapticParameters = " << std::endl;
    out << "{" << std::endl;
    out << "\tNss = " << synaptic_params.Nss << ";" << std::endl << std::endl;
    
    out << "\tA_P = " << synaptic_params.A_P << ";" << std::endl;
    out << "\tA_D = " << synaptic_params.A_D << ";" << std::endl;
    out << "\tT_0 = " << synaptic_params.T_0 << ";" << std::endl;
    out << "\tT_P = " << synaptic_params.T_P << ";" << std::endl;
    out << "\tT_D = " << synaptic_params.T_D << ";" << std::endl;
    out << "\tTAU_P = " << synaptic_params.TAU_P << ";" << std::endl;
    out << "\tTAU_D = " << synaptic_params.TAU_D << ";" << std::endl << std::endl;
    
    out << "\tBETA = " << synaptic_params.BETA << ";" << std::endl;
    out << "\tBETA_SUPERSYNAPSE = " << synaptic_params.BETA_SUPERSYNAPSE << ";" << std::endl << std::endl;

    out << "\tACTIVATION_THRESHOLD = " << synaptic_params.ACTIVATION_THRESHOLD << ";" << std::endl;
    out << "\tSUPERSYNAPSE_THRESHOLD = " << synaptic_params.SUPERSYNAPSE_THRESHOLD << ";" << std::endl;
    out << "\tWEIGHT_MAX = " << synaptic_params.WEIGHT_MAX << ";" << std::endl << std::endl;

    out << "};" << std::endl << std::endl;
    
    // write maturation parameters
    out << "gabaParameters = " << std::endl;
    out << "{" << std::endl;
    out << "\tE_GABA_MATURE = " << gaba_params.E_GABA_MATURE << ";" << std::endl;
    out << "\tE_GABA_IMMATURE = " << gaba_params.E_GABA_IMMATURE << ";" << std::endl << std::endl;
    
    out << "\tGABA_RATE_THRESHOLD = " << gaba_params.GABA_RATE_THRESHOLD << ";" << std::endl;
    out << "\tDEATH_RATE_THRESHOLD = " << gaba_params.DEATH_RATE_THRESHOLD << ";" << std::endl << std::endl;
    
    out << "\tRATE_WINDOW_SHORT = " << gaba_params.RATE_WINDOW_SHORT << ";" << std::endl;
    out << "\tRATE_WINDOW_LONG = " << gaba_params.RATE_WINDOW_LONG << ";" << std::endl << std::endl;
    
    out << "\tGABA_DOWN = " << gaba_params.GABA_DOWN << ";" << std::endl;
    out << "};" << std::endl << std::endl;
    
    out.close();
}
