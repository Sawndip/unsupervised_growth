#include "../PoolParallel.h"
#include <mpi.h>
#include "../ConfigurationNetworkGenerator.h"
#include "../ConfigurationGrowth.h"
#include <string>
#include <iostream>

using namespace std;

int main(int argc, char** argv)
{
	int rank; // MPI process rank
	
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    std::string networkConfigurationDir; // path to directory with network configuration file
    std::string trainingFile; // path to file with training neurons
    
    std::string growthConfigurationDir; // path to directory with growth configuration file
    std::string outputDir; // path to output directory
    
    if (argc > 1)
    {
        networkConfigurationDir = argv[1];
        trainingFile = argv[2];
        growthConfigurationDir = argv[3];
        outputDir = argv[4];
  
		if (rank == 0)
		{
			std::cout << "Path to directory with network configuration file: " << networkConfigurationDir << std::endl;
			std::cout << "Path to file with training neurons: " << trainingFile << std::endl;
			std::cout << "Path to directory with growth configuration file: " << growthConfigurationDir << std::endl;
			std::cout << "Path to output directory: " << outputDir << std::endl;
		}
    }
    else
		if (rank == 0)
			std::cerr << "Not enough command line arguments were not provided!" << std::endl;
    
    
    bool training = true; // indicator if training neurons are innervated
	int save_freq_short = 100; // saving frequency for the graph
	int save_freq_long = 500; // saving frequency for the network state backup
	
	
    ConfigurationGrowth growth_cfg;
    ConfigurationNetworkGenerator network_cfg;
    
    network_cfg.read_configuration((networkConfigurationDir + "network_parameters.cfg").c_str());
	growth_cfg.read_configuration((growthConfigurationDir + "growth_parameters.cfg").c_str());

	if (rank == 0)
	{
		network_cfg.print_configuration();
		growth_cfg.print_configuration();
	}

	PoolParallel pool(network_cfg, growth_cfg, outputDir);

	
    pool.new_chain_growth(networkConfigurationDir, trainingFile, training, save_freq_short, save_freq_long);
	
	MPI_Finalize();


	return 0;

}
