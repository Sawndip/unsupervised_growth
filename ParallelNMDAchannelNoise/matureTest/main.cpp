#include "../PoolParallel.h"
#include <mpi.h>
#include "../Configuration.h"
#include <string>

using namespace std;

int main(int argc, char** argv)
{
    
    
    std::string dataDir = "/home/eugene/Output/networks/gabaMaturation130417/"; // directory with data 
    std::string outputDir = "/home/eugene/Output/matureTest/gabaMaturation130417/"; // directory with data 
    
    std::string configurationFile = dataDir + "parameters.cfg"; // configuration file
    
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
    
    int num_trials = 20; // number to trials to perform

    pool.test_grown_chain(num_trials, dataDir, outputDir);
	
    double end_time = MPI_Wtime();

    if (rank == 0)
        printf("Execution time = %f\n", end_time - start_time);

	MPI_Finalize();


	return 0;

}
