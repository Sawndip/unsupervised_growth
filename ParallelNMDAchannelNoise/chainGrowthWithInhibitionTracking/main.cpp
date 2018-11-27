#include <mpi.h>
#include "NetworkGrowthSimulator.h"
#include "ConfigurationNetworkGrowth.h"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	
	int rank;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	std::string networkDir; // directory with network topology
	std::string fileTraining; // file name with HVC-RA training neurons
	std::string configFile; // file name with network growth parameters
	std::string outputDir; // directory where to save output data
	int num_trials; // number of trials for simulation

	double time_resolution_conductance; // time resolution of inhibitory conductance
	
	if ( argc > 5 )
	{
		networkDir = argv[1];
		fileTraining = argv[2];
		configFile = argv[3];
		outputDir = argv[4];
		num_trials = atoi(argv[5]);
		time_resolution_conductance = atof(argv[6]);
	}
	else
	{
		std::cerr << "Usage: <directory with network topology>, <file with training neurons>" << "\n"  
				  << "<file with parameters for network growth>, <directory where to save output data>" << "\n"
				  << "<number of trials>, <time resolution for inhibitory conductance>" << std::endl; 
				  
		return -1;
	}
	
	ConfigurationNetworkGrowth cfg;
	
	bool training = true;
	int save_freq_long = 5;
	
	// read topology parameters on master process
	if (rank == 0)
		cfg.read_configuration(configFile.c_str());
	
	NetworkGrowthSimulator networkSimulator;
	
	networkSimulator.new_chain_growth_with_inhibition_tracking(cfg, networkDir, fileTraining, outputDir, training, 
											 save_freq_long, num_trials, time_resolution_conductance);
	
	
	MPI_Finalize();
	
    return 0;
}
