#include "../PoolParallel.h"
#include <mpi.h>
#include "../ConfigurationNetworkGenerator.h"
#include "../ConfigurationGrowth.h"
#include <string>

using namespace std;

int main(int argc, char** argv)
{
	std::string outputDir; // output directory
    std::string dataDir; // directory with data
    int starting_trial; // trial id from which to start growth

    if (argc > 2)
    {
        dataDir = argv[1]; // 
        outputDir = argv[2];
        starting_trial = atoi(argv[3]); // 
        
        std::cout << "Data directory: " << dataDir << std::endl;
        std::cout << "Output directory: " << outputDir << std::endl;
        std::cout << "Starting trial: " << starting_trial << std::endl;
    }
    else
        std::cerr << "Not enough command line arguments were not provided!" << std::endl;
    
    
    std::string configurationFile = dataDir + "parameters.cfg"; // configuration file
    
    int rank; // MPI process rank
    
    ConfigurationGrowth growth_cfg;
    ConfigurationNetworkGenerator network_cfg;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    network_cfg.read_configuration((dataDir + "network_parameters.cfg").c_str());
	growth_cfg.read_configuration((dataDir + "growth_parameters.cfg").c_str());


    //if (rank == 0)
    //    cfg.print_configuration();

	PoolParallel pool(network_cfg, growth_cfg, outputDir);


	int save_freq_short = 40; // save frequency for graph update
	int save_freq_long = 100; // save frequency for network state update
	bool training = true;
	
	pool.print_simulation_parameters(); 
		
    pool.continue_chain_growth(dataDir, starting_trial, training, save_freq_short, save_freq_long);
	
	MPI_Finalize();

	return 0;

}
