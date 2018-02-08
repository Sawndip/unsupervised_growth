#include "TestNetwork.h"
#include "ConfigurationTestNetwork.h"
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
	std::string configFilename; // configuration file name
	std::string outputDir; // name of the output directory
	
	if ( argc > 2 )
	{
		configFilename = argv[1];
		outputDir = argv[2];
	}
	else
	{
		std::cerr << "Usage: <configFilename> <outputDir>" << std::endl;
		return -1;
	}
	
	// read STDP and Maturation parameters from configuration file 
	ConfigurationTestNetwork cfg;
	cfg.read_configuration(configFilename.c_str());
	
	struct SynapticParameters synaptic_params = cfg.get_synaptic_parameters();
    struct GabaParameters gaba_params = cfg.get_gaba_parameters();
    
    // noise parameters
    double white_noise_mean_soma = 0.0;
    double white_noise_std_soma = 0.15;
    double white_noise_mean_dend = 0.0;
    double white_noise_std_dend = 0.25;
    
    // simulation trial parameters
    double trial_duration = 1000;
    double timestep = 0.02;
	
	TestNetwork testNetwork; // network

	// set all parameters
	testNetwork.set_synaptic_parameters(synaptic_params);
	testNetwork.set_gaba_parameters(gaba_params);
	testNetwork.set_noise_parameters(white_noise_mean_soma, white_noise_std_soma,
									 white_noise_mean_dend, white_noise_std_dend);
	
	int N_TR = 4;
    int num_groups = 4;
    int num_neurons_in_group = 4;
    
	double Gei = 0.25;
	double dGie = 0.75;
	
	int save_freq = 25;

	testNetwork.create_second_layer_growth_sim(N_TR, num_groups, num_neurons_in_group, Gei, dGie);
	testNetwork.simulate(trial_duration, timestep, save_freq, outputDir);
	
	return 0;
}
