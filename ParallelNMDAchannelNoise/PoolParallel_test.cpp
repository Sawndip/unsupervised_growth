#include "PoolParallel.h"
#include <mpi.h>
#include "Configuration.h"
#include <string>
#include <iostream>

using namespace std;

int main(int argc, char** argv)
{

    std::string configurationFile; // configuration file

    if (argc > 1)
    {
        configurationFile = argv[1]; // configuration file
        std::cout << "Path to configuration file: " << configurationFile << std::endl;
    }
    else
        std::cerr << "Path to file with configuration was not provided!" << std::endl;
    
    int rank; // MPI process rank
    bool training = true; // indicator if training neurons are innervated

    Configuration cfg;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    cfg.read_configuration(configurationFile.c_str());

    //if (rank == 0)
    //    cfg.print_configuration();

	PoolParallel pool(cfg);

    //pool.initialize_coordinates();
	//pool.initialize_connections();

    //pool.print_invariable_connections();

	pool.print_simulation_parameters();
    
   	int save_freq_short = 40;
	int save_freq_long = 100;

    double start_time = MPI_Wtime();

    //pool.chain_growth_default(save_freq_short, save_freq_long);
    pool.chain_growth_with_clustered_training(training, save_freq_short, save_freq_long);
    //pool.chain_growth_with_no_RA2I_connections(save_freq_short, save_freq_long);
	
    double end_time = MPI_Wtime();

    if (rank == 0)
        printf("Execution time = %f\n", end_time - start_time);

	MPI_Finalize();


	return 0;

}
