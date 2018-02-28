#include "ConfigurationNetworkGrowth.h"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

void ConfigurationNetworkGrowth::read_configuration(const char* filename)
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
        const Setting& maturation_params_setting = root["maturationParameters"];

        maturation_params_setting.lookupValue("E_GABA_MATURE", maturation_params.E_GABA_MATURE);
        maturation_params_setting.lookupValue("E_GABA_IMMATURE", maturation_params.E_GABA_IMMATURE);
        
        maturation_params_setting.lookupValue("E_REST_MATURE", maturation_params.E_REST_MATURE);
        maturation_params_setting.lookupValue("E_REST_IMMATURE", maturation_params.E_REST_IMMATURE);
        
        maturation_params_setting.lookupValue("AD_MATURE", maturation_params.AD_MATURE);
        maturation_params_setting.lookupValue("AD_IMMATURE", maturation_params.AD_IMMATURE);
        
        maturation_params_setting.lookupValue("GK_MATURE", maturation_params.GK_MATURE);
        maturation_params_setting.lookupValue("GK_IMMATURE", maturation_params.GK_IMMATURE);
        
        maturation_params_setting.lookupValue("GNA_MATURE", maturation_params.GNA_MATURE);
        maturation_params_setting.lookupValue("GNA_IMMATURE", maturation_params.GNA_IMMATURE);
        
        maturation_params_setting.lookupValue("GCA_MATURE", maturation_params.GCA_MATURE);
        maturation_params_setting.lookupValue("GCA_IMMATURE", maturation_params.GCA_IMMATURE);
        
        maturation_params_setting.lookupValue("GCAK_MATURE", maturation_params.GCAK_MATURE);
        maturation_params_setting.lookupValue("GCAK_IMMATURE", maturation_params.GCAK_IMMATURE);
        
        maturation_params_setting.lookupValue("GSL_MATURE", maturation_params.GSL_MATURE);
        maturation_params_setting.lookupValue("GSL_IMMATURE", maturation_params.GSL_IMMATURE);
        
        maturation_params_setting.lookupValue("GDL_MATURE", maturation_params.GDL_MATURE);
        maturation_params_setting.lookupValue("GDL_IMMATURE", maturation_params.GDL_IMMATURE);
        
        maturation_params_setting.lookupValue("RC_MATURE", maturation_params.RC_MATURE);
        maturation_params_setting.lookupValue("RC_IMMATURE", maturation_params.RC_IMMATURE);
        
        
        maturation_params_setting.lookupValue("DEATH_RATE_THRESHOLD", maturation_params.DEATH_RATE_THRESHOLD);
        
        maturation_params_setting.lookupValue("RATE_WINDOW_LONG", maturation_params.RATE_WINDOW_LONG);
    }
    catch(const SettingNotFoundException &nfex)
    {
        std::cerr << "Maturation parameter setting is not found in configuration file!" << std::endl;
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
    
    // Get the connection parameters
    try
    {
        const Setting& connection_params_setting = root["connectionParameters"];

        connection_params_setting.lookupValue("Gie_max", connection_params.Gie_max);
        connection_params_setting.lookupValue("Gei_max", connection_params.Gei_max);
        connection_params_setting.lookupValue("delay_constant", connection_params.delay_constant);
    }
    catch(const SettingNotFoundException &nfex)
    {
        std::cerr << "Connection parameter setting is not found in configuration file!" << std::endl;
    }
}

void ConfigurationNetworkGrowth::print_configuration() const
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
    
    std::cout << "E_GABA_MATURE = " << maturation_params.E_GABA_MATURE << std::endl;
    std::cout << "E_GABA_IMMATURE = " << maturation_params.E_GABA_IMMATURE << std::endl << std::endl;
    
    std::cout << "E_REST_MATURE = " << maturation_params.E_REST_MATURE << std::endl;
    std::cout << "E_REST_IMMATURE = " << maturation_params.E_REST_IMMATURE << std::endl << std::endl;
    
    std::cout << "AD_MATURE = " << maturation_params.AD_MATURE << std::endl;
    std::cout << "AD_IMMATURE = " << maturation_params.AD_IMMATURE << std::endl << std::endl;
    
    std::cout << "GK_MATURE = " << maturation_params.GK_MATURE << std::endl;
    std::cout << "GK_IMMATURE = " << maturation_params.GK_IMMATURE << std::endl << std::endl;
    
    std::cout << "GNA_MATURE = " << maturation_params.GNA_MATURE << std::endl;
    std::cout << "GNA_IMMATURE = " << maturation_params.GNA_IMMATURE << std::endl << std::endl;
    
    std::cout << "GCA_MATURE = " << maturation_params.GCA_MATURE << std::endl;
    std::cout << "GCA_IMMATURE = " << maturation_params.GCA_IMMATURE << std::endl << std::endl;
    
    std::cout << "GCAK_MATURE = " << maturation_params.GCAK_MATURE << std::endl;
    std::cout << "GCAL_IMMATURE = " << maturation_params.GCAK_IMMATURE << std::endl << std::endl;
    
    std::cout << "GSL_MATURE = " << maturation_params.GSL_MATURE << std::endl;
    std::cout << "GSL_IMMATURE = " << maturation_params.GSL_IMMATURE << std::endl << std::endl;
    
    std::cout << "GDL_MATURE = " << maturation_params.GDL_MATURE << std::endl;
    std::cout << "GDL_IMMATURE = " << maturation_params.GDL_IMMATURE << std::endl << std::endl;
    
    std::cout << "RC_MATURE = " << maturation_params.RC_MATURE << std::endl;
    std::cout << "RC_IMMATURE = " << maturation_params.RC_IMMATURE << std::endl << std::endl;
    
    std::cout << "DEATH_RATE_THRESHOLD = " << maturation_params.DEATH_RATE_THRESHOLD << std::endl << std::endl;

    std::cout << "RATE_WINDOW_LONG = " << maturation_params.RATE_WINDOW_LONG << std::endl;
    
    
    // display all noise parameters read from file
    std::cout << std::endl << "Noise parameters read from configuration file: " << std::endl << std::endl;
    
    std::cout << "white_noise_mean_soma = " << noise_params.white_noise_mean_soma << std::endl;
    std::cout << "white_noise_std_soma = " << noise_params.white_noise_std_soma << std::endl;
    std::cout << "white_noise_mean_dend = " << noise_params.white_noise_mean_dend << std::endl;
    std::cout << "white_noise_std_dend = " << noise_params.white_noise_std_dend << std::endl;
    
    // display connection parameters read from file
    std::cout << std::endl << "Connection parameters read from configuration file: " << std::endl << std::endl;
    
    std::cout << "Gei_max = " << connection_params.Gei_max << std::endl;
    std::cout << "Gie_max = " << connection_params.Gie_max << std::endl;
    std::cout << "delay_constant = " << connection_params.delay_constant << std::endl;
   
}

void ConfigurationNetworkGrowth::write_configuration(const char* filename) const
{
    std::ofstream out;

    out.open(filename, std::ios::out);
    
    // write synaptic parameters
    out << "synapticParameters = " << std::endl;
    out << "{" << std::endl;
    out << "\tNss = " << synaptic_params.Nss << ";" << std::endl << std::endl;
    
    out << "\tA_P = " << std::fixed << std::setprecision(6) << synaptic_params.A_P << ";" << std::endl;
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
    out << "maturationParameters = " << std::endl;
    out << "{" << std::endl;
    out << "\tE_GABA_MATURE = " << maturation_params.E_GABA_MATURE << ";" << std::endl;
    out << "\tE_GABA_IMMATURE = " << maturation_params.E_GABA_IMMATURE << ";" << std::endl << std::endl;
    
    out << "\tE_REST_MATURE = " << maturation_params.E_REST_MATURE << ";" << std::endl;
    out << "\tE_REST_IMMATURE = " << maturation_params.E_REST_IMMATURE << ";" << std::endl << std::endl;
    
    out << "\tAD_MATURE = " << maturation_params.AD_MATURE << ";" << std::endl;
    out << "\tAD_IMMATURE = " << maturation_params.AD_IMMATURE << ";" << std::endl << std::endl;
    
    out << "\tGK_MATURE = " << maturation_params.GK_MATURE << ";" << std::endl;
    out << "\tGK_IMMATURE = " << maturation_params.GK_IMMATURE << ";" << std::endl << std::endl;
    
    out << "\tGNA_MATURE = " << maturation_params.GNA_MATURE << ";" << std::endl;
    out << "\tGNA_IMMATURE = " << maturation_params.GNA_IMMATURE << ";" << std::endl << std::endl;
    
    out << "\tGCA_MATURE = " << maturation_params.GCA_MATURE << ";" << std::endl;
    out << "\tGCA_IMMATURE = " << maturation_params.GCA_IMMATURE << ";" << std::endl << std::endl;
    
    out << "\tGCAK_MATURE = " << maturation_params.GCAK_MATURE << ";" << std::endl;
    out << "\tGCAK_IMMATURE = " << maturation_params.GCAK_IMMATURE << ";" << std::endl << std::endl;
    
    out << "\tGSL_MATURE = " << maturation_params.GSL_MATURE << ";" << std::endl;
    out << "\tGSL_IMMATURE = " << maturation_params.GSL_IMMATURE << ";" << std::endl << std::endl;
    
    out << "\tGDL_MATURE = " << maturation_params.GDL_MATURE << ";" << std::endl;
    out << "\tGDL_IMMATURE = " << maturation_params.GDL_IMMATURE << ";" << std::endl << std::endl;
    
    out << "\tRC_MATURE = " << maturation_params.RC_MATURE << ";" << std::endl;
    out << "\tRC_IMMATURE = " << maturation_params.RC_IMMATURE << ";" << std::endl << std::endl;
    
    out << "\tDEATH_RATE_THRESHOLD = " << maturation_params.DEATH_RATE_THRESHOLD << ";" << std::endl << std::endl;
    
    out << "\tRATE_WINDOW_LONG = " << maturation_params.RATE_WINDOW_LONG << ";" << std::endl << std::endl;
    
    out << "};" << std::endl << std::endl;
    
    // write noise parameters
    out << "noiseParameters = " << std::endl;
    out << "{" << std::endl;
    out << "\twhite_noise_mean_soma = " << noise_params.white_noise_mean_soma << ";" << std::endl;
    out << "\twhite_noise_std_soma = " << noise_params.white_noise_std_soma << ";" << std::endl;
    out << "\twhite_noise_mean_dend = " << noise_params.white_noise_mean_dend << ";" << std::endl;
    out << "\twhite_noise_std_dend = " << noise_params.white_noise_std_dend << ";" << std::endl;
    out << "};" << std::endl << std::endl;
    
    // write connection parameters
    out << "connectionParameters = " << std::endl;
    out << "{" << std::endl;
    out << "\tGei_max = " << connection_params.Gei_max << ";" << std::endl;
    out << "\tGie_max = " << connection_params.Gie_max << ";" << std::endl;
    out << "\tdelay_constant = " << connection_params.delay_constant << ";" << std::endl;
    
    out << "};" << std::endl << std::endl;

    out.close();
}
