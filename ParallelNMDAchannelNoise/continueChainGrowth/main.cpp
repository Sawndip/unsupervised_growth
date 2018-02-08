#include <mpi.h>
#include "NetworkGrowthSimulator.h"
#include "ConfigurationNetworkGrowth.h"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	
	int rank;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	std::string dataDir; // directory with data
	int starting_trial; // trial from each growth continues
	
	if ( argc > 2 )
	{
		dataDir = argv[1];
		starting_trial = atoi(argv[2]);
	}
	else
	{
		std::cerr << "Usage: <directory with data>, <starting trial>" << std::endl; 
				  
		return -1;
	}
	
	bool training = true;
	int save_freq_short = 1;
	int save_freq_long = 10;
	
	NetworkGrowthSimulator networkSimulator;
	
	networkSimulator.continue_chain_growth(dataDir, starting_trial, training, save_freq_short, save_freq_long);
	
	MPI_Finalize();
	
    return 0;
}
