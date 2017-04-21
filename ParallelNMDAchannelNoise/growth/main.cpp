#include "../PoolParallel.h"
#include <mpi.h>
#include "../Configuration.h"
#include <string>
#include <iostream>

using namespace std;

int main(int argc, char** argv)
{

    std::string configurationFile; // configuration file
    std::string networkDir; // directory which contains network with fixed connections

    if (argc > 2)
    {
        configurationFile = argv[1]; // configuration file
        networkDir = argv[2]; 
        
        std::cout << "Path to configuration file: " << configurationFile << std::endl;
        std::cout << "Directory with fixed network: " << networkDir << std::endl;
        
    }
    else
        std::cerr << "Not enough command line arguments were not provided!" << std::endl;
    
    int rank; // MPI process rank
    bool training = true; // indicator if training neurons are innervated

    Configuration cfg;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    cfg.read_configuration(configurationFile.c_str());

	PoolParallel pool(cfg);

	pool.print_simulation_parameters();
    
    pool.read_fixed_network(networkDir); // read network from files
    
   	int save_freq_short = 40;
	int save_freq_long = 100;

    pool.chain_growth(training, save_freq_short, save_freq_long);
	
	MPI_Finalize();


	return 0;

}
