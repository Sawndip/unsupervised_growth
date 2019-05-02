#include <mpi.h>
#include "NetworkGrowthSimulator.h"
#include "ConfigurationNetworkGrowth.h"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	
	int rank;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	std::string networkDir; // directory with network topology
	int trial_number; // trial number with data
	int num_test_trials; // number of testing trials
	std::string outputDir; // directory where to save output data
	
	
	if ( argc > 4 )
	{
		networkDir = argv[1];
		trial_number = atoi(argv[2]);
		num_test_trials = atoi(argv[3]);
		outputDir = argv[4];
	}
	else
	{
		std::cerr << "Usage: <directory with network>, <trial number>" << "\n"  
				  <<  "<number of testing trials>,  <directory where to save output data>" << std::endl; 
				  
		return -1;
	}
	
	NetworkGrowthSimulator networkSimulator;
	
	networkSimulator.test_chain(networkDir, trial_number, num_test_trials, outputDir);
	
	MPI_Finalize();
	
    return 0;
}
