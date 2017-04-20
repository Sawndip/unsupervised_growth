#include "../PoolParallel.h"
#include <mpi.h>
#include "../Configuration.h"
#include <string>

using namespace std;

int main(int argc, char** argv)
{
    std::string dataDir; // directory with data

    if (argc > 1)
    {
        dataDir = argv[1]; // 
        std::cout << "Data directory: " << dataDir << std::endl;
    }
    else
        std::cerr << "Data directory was not provided!" << std::endl;
    
    
    std::string configurationFile = dataDir + "parameters.cfg"; // configuration file
    
    int rank; // MPI process rank
    Configuration cfg;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    cfg.read_configuration(configurationFile.c_str());

    //if (rank == 0)
    //    cfg.print_configuration();

	PoolParallel pool(cfg);

	int save_freq_short = 40; // save frequency for graph update
	int save_freq_long = 100; // save frequency for network state update

	pool.print_simulation_parameters(); 
    pool.continue_growth(dataDir, save_freq_short, save_freq_long);
	
	MPI_Finalize();

	return 0;

}
