#include <mpi.h>
#include "NetworkGrowthSimulator.h"
#include "ConfigurationNetworkGrowth.h"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	
	int rank;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	std::string networkDir; // directory with network topology
	std::string fileTraining; // file with training neurons
	std::string outputDir; // directory where to save output data
	
	
	if ( argc > 3 )
	{
		networkDir = argv[1];
		fileTraining = argv[2];
		outputDir = argv[3];
	}
	else
	{
		std::cerr << "Usage: <directory with network>, <trial number>" << "\n"  
				  <<  "<number of testing trials>,  <directory where to save output data>" << std::endl; 
				  
		return -1;
	}
	
	NetworkGrowthSimulator networkSimulator;

	double white_noise_mean_soma = 0.00;
	double white_noise_std_soma = 0.10;
	double white_noise_mean_dend = 0.00;
	double white_noise_std_dend = 0.05;

	double Gie_max = 3.0;
	double Gei_max = 0.35;
	double delay_constant = 1.0;

	struct ConnectionParameters connection_parameters;
	struct NoiseParameters noise_parameters;

	noise_parameters.white_noise_mean_soma = white_noise_mean_soma;
	noise_parameters.white_noise_std_soma = white_noise_std_soma;
	noise_parameters.white_noise_mean_dend = white_noise_mean_dend;
	noise_parameters.white_noise_std_dend = white_noise_std_dend;

	connection_parameters.Gie_max = Gie_max;
	connection_parameters.Gei_max = Gei_max;
	connection_parameters.delay_constant = delay_constant;

	networkSimulator.test_inhAndExc_response(connection_parameters, noise_parameters,
													networkDir, fileTraining, outputDir);
	
	MPI_Finalize();
	
    return 0;
}
