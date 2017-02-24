#include "Configuration.h"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

void Configuration::read_configuration(const char* filename)
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
    
    // Get directory for data output
    if (!root.lookupValue("outputDirectory", outputDirectory))
        std::cerr << "No 'outputDirectory' setting in configuration file!" << std::endl;
    
    // Get the network parameters
    try
    {
        const Setting& network_params_setting = root["networkParameters"];

        network_params_setting.lookupValue("N_RA", network_params.N_RA);
        network_params_setting.lookupValue("N_TR", network_params.N_TR);
        network_params_setting.lookupValue("N_I", network_params.N_I);
    }
    catch(const SettingNotFoundException &nfex)
    {
        std::cerr << "Network parameter setting is not found in configuration file!" << std::endl;
    }
    // Get the synaptic parameters
    try
    {
        const Setting& synaptic_params_setting = root["synapticParameters"];

        synaptic_params_setting.lookupValue("Nss", synaptic_params.Nss);
        
        synaptic_params_setting.lookupValue("R", synaptic_params.R);
          
        synaptic_params_setting.lookupValue("A_P", synaptic_params.A_P);
        synaptic_params_setting.lookupValue("A_D", synaptic_params.A_D);
        synaptic_params_setting.lookupValue("F_0", synaptic_params.F_0);
        synaptic_params_setting.lookupValue("T_P", synaptic_params.T_P);
        synaptic_params_setting.lookupValue("T_D", synaptic_params.T_D);
        synaptic_params_setting.lookupValue("TAU_P", synaptic_params.TAU_P);
        synaptic_params_setting.lookupValue("TAU_D", synaptic_params.TAU_D);
          
        synaptic_params_setting.lookupValue("BETA", synaptic_params.BETA);
        synaptic_params_setting.lookupValue("BETA_SUPERSYNAPSE", synaptic_params.BETA_SUPERSYNAPSE);
          
        synaptic_params_setting.lookupValue("ACTIVATION_THRESHOLD", synaptic_params.ACTIVATION_THRESHOLD);
        synaptic_params_setting.lookupValue("SUPERSYNAPSE_THRESHOLD", synaptic_params.SUPERSYNAPSE_THRESHOLD);
        synaptic_params_setting.lookupValue("WEIGHT_MAX", synaptic_params.WEIGHT_MAX);
          
        synaptic_params_setting.lookupValue("STDP_WINDOW", synaptic_params.STDP_WINDOW);
    }
    catch(const SettingNotFoundException &nfex)
    {
        std::cerr << "Synaptic parameter setting is not found in configuration file!" << std::endl;
    }


    // Get the spatial parameters
    try
    {
        const Setting& spatial_params_setting = root["spatialParameters"];

        // Only output the record if all of the expected fields are present.
        spatial_params_setting.lookupValue("SIDE", spatial_params.SIDE);
          
        spatial_params_setting.lookupValue("MIN_INTERNEURON_DISTANCE", spatial_params.MIN_INTERNEURON_DISTANCE);
        
        
        spatial_params_setting.lookupValue("A_RA2I", spatial_params.A_RA2I);
        spatial_params_setting.lookupValue("SIGMA_RA2I", spatial_params.SIGMA_RA2I);
        
        spatial_params_setting.lookupValue("B_I2RA", spatial_params.B_I2RA);
        spatial_params_setting.lookupValue("SIGMA_I2RA", spatial_params.SIGMA_I2RA);
        
        spatial_params_setting.lookupValue("Gei_mean", spatial_params.Gei_mean);
        spatial_params_setting.lookupValue("Gei_std", spatial_params.Gei_std);
        spatial_params_setting.lookupValue("Gie_mean", spatial_params.Gie_mean);
        spatial_params_setting.lookupValue("Gie_std", spatial_params.Gie_std);
        
    }
    catch(const SettingNotFoundException &nfex)
    {
        std::cerr << "Spatial parameter setting is not found in configuration file!" << std::endl;
    }

    // Get the maturation parameters
    try
    {
        const Setting& gaba_params_setting = root["gabaParameters"];

        gaba_params_setting.lookupValue("E_GABA_MATURE", gaba_params.E_GABA_MATURE);
        gaba_params_setting.lookupValue("E_GABA_IMMATURE", gaba_params.E_GABA_IMMATURE);
        
        gaba_params_setting.lookupValue("GABA_RATE_THRESHOLD", gaba_params.GABA_RATE_THRESHOLD);
        gaba_params_setting.lookupValue("MATURATION_RATE_THRESHOLD", gaba_params.MATURATION_RATE_THRESHOLD);
        gaba_params_setting.lookupValue("DEATH_RATE_THRESHOLD", gaba_params.DEATH_RATE_THRESHOLD);
        
        gaba_params_setting.lookupValue("RATE_WINDOW_SHORT", gaba_params.RATE_WINDOW_SHORT);
        gaba_params_setting.lookupValue("RATE_WINDOW_LONG", gaba_params.RATE_WINDOW_LONG);
        
        gaba_params_setting.lookupValue("GABA_DOWN", gaba_params.GABA_DOWN);
        
    }
    catch(const SettingNotFoundException &nfex)
    {
        std::cerr << "Maturation parameter setting is not found in configuration file!" << std::endl;
    }
    
    // Get the time parameters
    try
    {
        const Setting& time_params_setting = root["timeParameters"];

        time_params_setting.lookupValue("timestep", time_params.timestep);
        time_params_setting.lookupValue("network_update_frequency", time_params.network_update_frequency);
        time_params_setting.lookupValue("trial_duration", time_params.trial_duration);
        time_params_setting.lookupValue("WAITING_TIME", time_params.WAITING_TIME);
        
    }
    catch(const SettingNotFoundException &nfex)
    {
        std::cerr << "Time parameter setting is not found in configuration file!" << std::endl;
    }
    
    // Get the noise parameters
    try
    {
        const Setting& noise_params_setting = root["noiseParameters"];

        noise_params_setting.lookupValue("white_noise_mean_soma", noise_params.white_noise_mean_soma);
        noise_params_setting.lookupValue("white_noise_std_soma", noise_params.white_noise_std_soma);
        noise_params_setting.lookupValue("white_noise_mean_dend", noise_params.white_noise_mean_dend);
        noise_params_setting.lookupValue("white_noise_std_dend", noise_params.white_noise_std_dend);
        
    }
    catch(const SettingNotFoundException &nfex)
    {
        std::cerr << "Noise parameter setting is not found in configuration file!" << std::endl;
    }
}

void Configuration::print_configuration()
{
    std::cout << std::endl << "Directory for data output = " << outputDirectory << std::endl; 

    // display all network parameters read from file
    std::cout << std::endl << "Network parameters read from configuration file: " << std::endl << std::endl;
    
    std::cout << "N_RA = " << network_params.N_RA << std::endl;
    std::cout << "N_TR = " << network_params.N_TR << std::endl;
    std::cout << "N_I = " << network_params.N_I << std::endl;
    
    // display all synaptic parameters read from file
    std::cout << std::endl << "Synaptic parameters read from configuration file: " << std::endl << std::endl;
 
    std::cout << "Nss = " << synaptic_params.Nss << std::endl << std::endl;

    std::cout << "R = " << synaptic_params.R << std::endl << std::endl;
    
    std::cout << "A_P = " << synaptic_params.A_P << std::endl;
    std::cout << "A_D = " << synaptic_params.A_D << std::endl;
    std::cout << "F_0 = " << synaptic_params.F_0 << std::endl;
    std::cout << "T_P = " << synaptic_params.T_P << std::endl;
    std::cout << "T_D = " << synaptic_params.T_D << std::endl;
    std::cout << "TAU_P = " << synaptic_params.TAU_P << std::endl;
    std::cout << "TAU_D = " << synaptic_params.TAU_D << std::endl << std::endl;
    
    std::cout << "BETA = " << synaptic_params.BETA << std::endl;
    std::cout << "BETA_SUPERSYNAPSE = " << synaptic_params.BETA_SUPERSYNAPSE << std::endl << std::endl;
    
    std::cout << "ACTIVATION_THRESHOLD = " << synaptic_params.ACTIVATION_THRESHOLD << std::endl;
    std::cout << "SUPERSYNAPSE_THRESHOLD = " << synaptic_params.SUPERSYNAPSE_THRESHOLD << std::endl;
    std::cout << "WEIGHT_MAX = " << synaptic_params.WEIGHT_MAX << std::endl << std::endl;

    std::cout << "STDP_WINDOW = " << synaptic_params.STDP_WINDOW << std::endl;

    // display all spatial parameters read from file
    std::cout << std::endl << "Spatial parameters read from configuration file: " << std::endl << std::endl;
    
    std::cout << "SIDE = " << spatial_params.SIDE << std::endl;
    std::cout << "MIN_INTERNEURON_DISTANCE = " << spatial_params.MIN_INTERNEURON_DISTANCE << std::endl << std::endl;
    
    std::cout << "A_RA2I = " << spatial_params.A_RA2I << std::endl;
    std::cout << "SIGMA_RA2I = " << spatial_params.SIGMA_RA2I << std::endl;
    std::cout << "B_I2RA = " << spatial_params.B_I2RA << std::endl;
    std::cout << "SIGMA_I2RA = " << spatial_params.SIGMA_I2RA << std::endl << std::endl;

    std::cout << "Gei_mean = " << spatial_params.Gei_mean << std::endl;
    std::cout << "Gei_std = " << spatial_params.Gei_std << std::endl;
    std::cout << "Gie_mean = " << spatial_params.Gie_mean << std::endl;
    std::cout << "Gie_std = " << spatial_params.Gie_std << std::endl;
  
    // display all maturation parameters read from file
    std::cout << std::endl << "Maturation parameters read from configuration file: " << std::endl << std::endl;
    
    std::cout << "E_GABA_MATURE = " << gaba_params.E_GABA_MATURE << std::endl;
    std::cout << "E_GABA_IMMATURE = " << gaba_params.E_GABA_IMMATURE << std::endl << std::endl;
    
    std::cout << "GABA_RATE_THRESHOLD = " << gaba_params.GABA_RATE_THRESHOLD << std::endl;
    std::cout << "MATURATION_RATE_THRESHOLD = " << gaba_params.MATURATION_RATE_THRESHOLD << std::endl;
    std::cout << "DEATH_RATE_THRESHOLD = " << gaba_params.DEATH_RATE_THRESHOLD << std::endl << std::endl;

    std::cout << "RATE_WINDOW_SHORT = " << gaba_params.RATE_WINDOW_SHORT << std::endl;
    std::cout << "RATE_WINDOW_LONG = " << gaba_params.RATE_WINDOW_LONG << std::endl << std::endl;

    std::cout << "GABA_DOWN = " << gaba_params.GABA_DOWN << std::endl;
  
    // display all time parameters read from file
    std::cout << std::endl << "Time parameters read from configuration file: " << std::endl << std::endl;
    
    std::cout << "timestep = " << time_params.timestep << std::endl;
    std::cout << "network_update_frequency = " << time_params.network_update_frequency << std::endl;
    std::cout << "trial_duration = " << time_params.trial_duration << std::endl;
    std::cout << "WAITING_TIME = " << time_params.WAITING_TIME << std::endl;
  
    // display all noise parameters read from file
    std::cout << std::endl << "Noise parameters read from configuration file: " << std::endl << std::endl;
    
    std::cout << "white_noise_mean_soma = " << noise_params.white_noise_mean_soma << std::endl;
    std::cout << "white_noise_std_soma = " << noise_params.white_noise_std_soma << std::endl;
    std::cout << "white_noise_mean_dend = " << noise_params.white_noise_mean_dend << std::endl;
    std::cout << "white_noise_std_dend = " << noise_params.white_noise_std_dend << std::endl;
  
}

void Configuration::write_configuration(const char* filename) const
{
    std::ofstream out;

    out.open(filename, std::ios::out);
    
    // write output directory
    out << "outputDirectory = \"" << outputDirectory << "\";" << std::endl << std::endl;

    // write network parameters
    out << "networkParameters = " << std::endl;
    out << "{" << std::endl;
    out << "\tN_RA = " << network_params.N_RA << ";" << std::endl;
    out << "\tN_TR = " << network_params.N_TR << ";" << std::endl;
    out << "\tN_I = " << network_params.N_I << ";" << std::endl;
    out << "};" << std::endl << std::endl;
    
    // write synaptic parameters
    out << "synapticParameters = " << std::endl;
    out << "{" << std::endl;
    out << "\tNss = " << synaptic_params.Nss << ";" << std::endl << std::endl;
    
    out << "\tR = " << std::fixed << std::setprecision(6) << synaptic_params.R << ";" << std::endl << std::endl;

    out << "\tA_P = " << synaptic_params.A_P << ";" << std::endl;
    out << "\tA_D = " << synaptic_params.A_D << ";" << std::endl;
    out << "\tF_0 = " << synaptic_params.F_0 << ";" << std::endl;
    out << "\tT_P = " << synaptic_params.T_P << ";" << std::endl;
    out << "\tT_D = " << synaptic_params.T_D << ";" << std::endl;
    out << "\tTAU_P = " << synaptic_params.TAU_P << ";" << std::endl;
    out << "\tTAU_D = " << synaptic_params.TAU_D << ";" << std::endl << std::endl;
    
    out << "\tBETA = " << synaptic_params.BETA << ";" << std::endl;
    out << "\tBETA_SUPERSYNAPSE = " << synaptic_params.BETA_SUPERSYNAPSE << ";" << std::endl << std::endl;

    out << "\tACTIVATION_THRESHOLD = " << synaptic_params.ACTIVATION_THRESHOLD << ";" << std::endl;
    out << "\tSUPERSYNAPSE_THRESHOLD = " << synaptic_params.SUPERSYNAPSE_THRESHOLD << ";" << std::endl;
    out << "\tWEIGHT_MAX = " << synaptic_params.WEIGHT_MAX << ";" << std::endl << std::endl;

    out << "\tSTDP_WINDOW = " << synaptic_params.STDP_WINDOW << ";" << std::endl;
    out << "};" << std::endl << std::endl;

    // write spatial parameters
    out << "spatialParameters = " << std::endl;
    out << "{" << std::endl;
    out << "\tSIDE = " << spatial_params.SIDE << ";" << std::endl;
    out << "\tMIN_INTERNEURON_DISTANCE = " << spatial_params.MIN_INTERNEURON_DISTANCE << ";" << std::endl << std::endl;
    
    out << "\tA_RA2I = " << spatial_params.A_RA2I << ";" << std::endl;
    out << "\tSIGMA_RA2I = " << spatial_params.SIGMA_RA2I << ";" << std::endl;
    out << "\tB_I2RA = " << spatial_params.B_I2RA << ";" << std::endl;
    out << "\tSIGMA_I2RA = " << spatial_params.SIGMA_I2RA << ";" << std::endl << std::endl;
    
    out << "\tGei_mean = " << spatial_params.Gei_mean << ";" << std::endl;
    out << "\tGei_std = " << spatial_params.Gei_std << ";" << std::endl;
    out << "\tGie_mean = " << spatial_params.Gie_mean << ";" << std::endl;
    out << "\tGie_std = " << spatial_params.Gie_std << ";" << std::endl << std::endl;
    out << "};" << std::endl << std::endl;
    
    // write maturation parameters
    out << "gabaParameters = " << std::endl;
    out << "{" << std::endl;
    out << "\tE_GABA_MATURE = " << gaba_params.E_GABA_MATURE << ";" << std::endl;
    out << "\tE_GABA_IMMATURE = " << gaba_params.E_GABA_IMMATURE << ";" << std::endl << std::endl;
    
    out << "\tGABA_RATE_THRESHOLD = " << gaba_params.GABA_RATE_THRESHOLD << ";" << std::endl;
    out << "\tMATURATION_RATE_THRESHOLD = " << gaba_params.MATURATION_RATE_THRESHOLD << ";" << std::endl;
    out << "\tDEATH_RATE_THRESHOLD = " << gaba_params.DEATH_RATE_THRESHOLD << ";" << std::endl << std::endl;
    
    out << "\tRATE_WINDOW_SHORT = " << gaba_params.RATE_WINDOW_SHORT << ";" << std::endl;
    out << "\tRATE_WINDOW_LONG = " << gaba_params.RATE_WINDOW_LONG << ";" << std::endl << std::endl;
    
    out << "\tGABA_DOWN = " << gaba_params.GABA_DOWN << ";" << std::endl;
    out << "};" << std::endl << std::endl;
    
    // write time parameters
    out << "timeParameters = " << std::endl;
    out << "{" << std::endl;
    out << "\ttimestep = " << time_params.timestep << ";" << std::endl;
    out << "\tnetwork_update_frequency = " << time_params.network_update_frequency << ";" << std::endl;
    out << "\ttrial_duration = " << time_params.trial_duration << ";" << std::endl;
    out << "\tWAITING_TIME = " << time_params.WAITING_TIME << ";" << std::endl;
    out << "};" << std::endl << std::endl;

    // write noise parameters
    out << "noiseParameters = " << std::endl;
    out << "{" << std::endl;
    out << "\twhite_noise_mean_soma = " << noise_params.white_noise_mean_soma << ";" << std::endl;
    out << "\twhite_noise_std_soma = " << noise_params.white_noise_std_soma << ";" << std::endl;
    out << "\twhite_noise_mean_dend = " << noise_params.white_noise_mean_dend << ";" << std::endl;
    out << "\twhite_noise_std_dend = " << noise_params.white_noise_std_dend << ";" << std::endl;
    out << "};" << std::endl << std::endl;

    out.close();
}