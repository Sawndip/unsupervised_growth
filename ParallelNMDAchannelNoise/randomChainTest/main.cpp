#include "../PoolParallel.h"
#include <mpi.h>
#include "../Configuration.h"
#include <string>

using namespace std;

int main(int argc, char** argv)
{
    
    std::string configurationFile = "/home/eugene/Output/networks/RandomChainTest220217/parameters.cfg"; // configuration file
    
    int rank; // MPI process rank
    Configuration cfg;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    cfg.read_configuration(configurationFile.c_str());

    //if (rank == 0)
    //    cfg.print_configuration();

	PoolParallel pool(cfg);

	pool.print_simulation_parameters();
    
    double start_time = MPI_Wtime();
   
    int num_trials = 10;
    int num_layers = 6; // number of synfire chain groups

    pool.initialize_coordinates();
    pool.initialize_random_chain_connections(num_layers);
	
    double end_time = MPI_Wtime();

    if (rank == 0)
        printf("Execution time = %f\n", end_time - start_time);

	MPI_Finalize();


	return 0;

}
