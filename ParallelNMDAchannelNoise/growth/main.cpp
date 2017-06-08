#include "../PoolParallel.h"
#include <mpi.h>
#include "../ConfigurationNetworkGenerator.h"
#include "../ConfigurationGrowth.h"
#include <string>
#include <iostream>

using namespace std;

int main(int argc, char** argv)
{

    std::string growthConfigurationDir; // path to directory with growth configuration file
    std::string networkConfigurationDir; // path to directory with network configuration file
    std::string outputDir; // path to output directory
    
    if (argc > 1)
    {
       
        networkConfigurationDir = argv[1];
        growthConfigurationDir = argv[2];
        outputDir = argv[3];
  
		std::cout << "Path to directory with network configuration file: " << networkConfigurationDir << std::endl;
        std::cout << "Path to directory with growth configuration file: " << growthConfigurationDir << std::endl;
        std::cout << "Path to output directory: " << outputDir << std::endl;
    }
    else
        std::cerr << "Not enough command line arguments were not provided!" << std::endl;
    
    int rank; // MPI process rank
    bool training = true; // indicator if training neurons are innervated
	int save_freq_short = 20; // saving frequency for the graph
	int save_freq_long = 50; // saving frequency for the network state backup
	
	
    ConfigurationGrowth growth_cfg;
    ConfigurationNetworkGenerator network_cfg;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    network_cfg.read_configuration((networkConfigurationDir + "network_parameters.cfg").c_str());
	growth_cfg.read_configuration((growthConfigurationDir + "growth_parameters.cfg").c_str());

	//if (rank == 0)
	//{
	//	network_cfg.print_configuration();
	//	growth_cfg.print_configuration();
	//}

	PoolParallel pool(network_cfg, growth_cfg, outputDir);

	
    pool.new_chain_growth(networkConfigurationDir, training, save_freq_short, save_freq_long);
	
	MPI_Finalize();


	return 0;

}
