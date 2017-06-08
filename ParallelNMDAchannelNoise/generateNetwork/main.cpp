#include "../ConfigurationNetworkGenerator.h"
#include "../NetworkGenerator.h"
#include <string>
#include <iostream>

using namespace std;

int main(int argc, char** argv)
{
	int N_RA = 75; // number of HVC(RA) neurons
	int N_I = 25; // number of HVC(I) neurons
	int N_TR = 4; // number of training HVC(RA) neurons

    std::string configurationFile; // configuration file
    std::string outputDirectory; // directory to which write network information
  
    if (argc > 2)
    {
        configurationFile = argv[1]; // configuration file
		outputDirectory = argv[2];
        
        std::cout << "Path to configuration file: " << configurationFile << std::endl;
		std::cout << "Path to output directory: " << outputDirectory << std::endl;
       
    }
    else
    {
        std::cerr << "Not enough command line arguments were provided!" << std::endl;
		return -1;
	}
    
    int rank; // MPI process rank
    bool training = true; // indicator if training neurons are innervated

    ConfigurationNetworkGenerator cfg;

    cfg.read_configuration(configurationFile.c_str());
    
    cfg.print_configuration();

    Poisson_noise noise_generator;
    
    unsigned seed = 1991;
    
    noise_generator.set_seed(seed);

	NetworkGenerator network_gen(cfg, &noise_generator);

    network_gen.generate_network_with_clustered_training(N_RA, N_TR, N_I); // generate network
    //network_gen.generate_network_with_dispersed_training(); // generate network and write to to directory networkDir
    network_gen.write_invariant_network_to_directory(outputDirectory);
    network_gen.write_alterable_network_to_directory("_initial", outputDirectory); // write initial network to directory
    
	

	return 0;

}
