#include <mpi.h>
#include "NetworkGrowthSimulator.h"
#include "ConfigurationNetworkTopology.h"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	
	int rank;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if (rank == 0)
	{
		std::string configFile; // file name with network topology parameters
		std::string networkDir; // directory where to save generated network
		
		if ( argc > 2 )
		{
			configFile = argv[1];
			networkDir = argv[2];
		}
		else
		{
			std::cerr << "Usage: <file name with network topology parameters>, <directory where to save generated network>" 
					  << std::endl;
					  
			return -1;
		}
		
		// read topology parameters
		ConfigurationNetworkTopology cfg;
		
		cfg.read_configuration(configFile.c_str());
		
		struct TopologyParameters top_parameters = cfg.get_topology_parameters();
		
		NetworkGrowthSimulator networkSimulator;
		
		int N_RA = 1000;
		int N_I = 200;
		int N_TR = 10;
	
		
		// generate network topology
		networkSimulator.generate_network_topology(N_RA, N_I, N_TR, top_parameters, networkDir);
	}
	
	MPI_Finalize();
	
    return 0;
}
