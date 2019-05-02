#include <mpi.h>
#include "NetworkGrowthSimulator.h"
#include "ConfigurationNetworkGrowth.h"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	
	int rank;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	std::string dataDir; // directory with data
	std::string outDir; // output directory
	int starting_trial; // trial from each growth continues
	
	if ( argc > 2 )
	{
		dataDir = argv[1];
		outDir = argv[2];
		starting_trial = atoi(argv[3]);
	}
	else
	{
		std::cerr << "Usage: <directory with data>, <output directory>, <starting trial>" << std::endl; 
				  
		return -1;
	}
	
	bool training = true;
	int save_freq_short = 50;
	int save_freq_long = 100;
	
	NetworkGrowthSimulator networkSimulator;
	
	networkSimulator.continue_chain_growth(dataDir, outDir, starting_trial, training, save_freq_short, save_freq_long);
	
	MPI_Finalize();
	
    return 0;
}
