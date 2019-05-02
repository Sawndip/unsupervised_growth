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
	
	if ( argc > 4 )
	{
		networkDir = argv[1];
		fileTraining = argv[2];
		configFile = argv[3];
		outputDir = argv[4];
	}
	else
	{
		std::cerr << "Usage: <directory with network topology>, <file with training neurons>" << "\n"  
				  << "<file with parameters for network growth>, <directory where to save output data>" << std::endl; 
				  
		return -1;
	}
	
	ConfigurationNetworkGrowth cfg;
	
	bool training = true;
	int save_freq_short = 75;
	int save_freq_long = 150;
	
	// read topology parameters on master process
	if (rank == 0)
		cfg.read_configuration(configFile.c_str());
	
	NetworkGrowthSimulator networkSimulator;
	
	networkSimulator.new_chain_growth(cfg, networkDir, fileTraining, outputDir, training, 
											save_freq_short, save_freq_long);
	
	MPI_Finalize();
	
    return 0;
}
