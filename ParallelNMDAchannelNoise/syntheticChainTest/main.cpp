#include <mpi.h>
#include "NetworkGrowthSimulator.h"
#include "ConfigurationNetworkGrowth.h"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	
	int rank;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	std::string networkDir; // directory with network topology
	std::string outputDir; // directory where to save output data
	std::string configFile; // file name with network growth parameters

	
	int num_layers; // number of layers in synfire chain
	int num_group; // number of neurons in each layer
	double probability; // probability of connecting neuron in the next layer
	double Gee; // strength of connection between neurons in the chain
	
	int num_test_trials; // number of testing trials
	
	if ( argc > 7 )
	{
		num_layers = atoi(argv[1]);
		num_group = atoi(argv[2]);
		probability = atof(argv[3]);
		Gee = atof(argv[4]);
		
		networkDir = argv[5];
		configFile = argv[6];
		
		num_test_trials = atoi(argv[7]);
		
		outputDir = argv[8];
		
	}
	else
	{
		std::cerr << "Usage: <num_layers>, <num_group>, <probability>, <Gee>," << "\n"  
				  <<  "<directory with network topology>, <path to configuration file, " << "\n"
				  << "<num_test_trials>, <directory where to save output data>" << std::endl; 
				  
		return -1;
	}
	
	NetworkGrowthSimulator networkSimulator;
	ConfigurationNetworkGrowth cfg;
	
	
	// read topology parameters on master process
	if (rank == 0)
		cfg.read_configuration(configFile.c_str());
	
	networkSimulator.test_synthetic_chain(cfg, num_layers, num_group, probability, Gee, networkDir, num_test_trials, outputDir);
	
	MPI_Finalize();
	
    return 0;
}
