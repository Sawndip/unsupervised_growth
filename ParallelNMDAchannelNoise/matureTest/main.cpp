#include "../PoolParallel.h"
#include <mpi.h>
#include "../ConfigurationNetworkGenerator.h"
#include "../ConfigurationGrowth.h"
#include <string>

using namespace std;

int main(int argc, char** argv)
{
    
    
    std::string dataDir = "/mnt/hodgkin_home/eugene/Output/networks/networkTest/"; // directory with data 
    std::string outputDir = "/mnt/hodgkin_home/eugene/Output/matureTest/test/"; // directory with data 
    
    int starting_trial = 500; // trial number defining network state
    
   
    int rank; // MPI process rank
    

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    ConfigurationGrowth growth_cfg;
    ConfigurationNetworkGenerator network_cfg;
    
    
    network_cfg.read_configuration((dataDir + "network_parameters.cfg").c_str());
	growth_cfg.read_configuration((dataDir + "growth_parameters.cfg").c_str());


    //if (rank == 0)
    //    cfg.print_configuration();

	PoolParallel pool(network_cfg, growth_cfg, outputDir);

	pool.print_simulation_parameters();
    
    double start_time = MPI_Wtime();
    
    int num_trials = 20; // number to trials to perform

    pool.test_grown_chain(num_trials, dataDir, starting_trial, outputDir);
	
    double end_time = MPI_Wtime();

    if (rank == 0)
        printf("Execution time = %f\n", end_time - start_time);

	MPI_Finalize();


	return 0;

}
