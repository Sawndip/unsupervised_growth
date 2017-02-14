#include "../PoolParallel.h"
#include <mpi.h>
#include "../structures.h"
#include <fstream>

using namespace std;

void write_execution_time(int N_RA, int N_I, int MPI_size, double timestep, double network_update_frequency, double exec_time, const char* filename)
{
    std::ofstream out;
      
    out.open(filename, std::ios::out | std::ios::binary );
        
    out.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));
    out.write(reinterpret_cast<char *>(&N_I), sizeof(N_I));
    out.write(reinterpret_cast<char *>(&MPI_size), sizeof(MPI_size));
    out.write(reinterpret_cast<char *>(&timestep), sizeof(timestep));
    out.write(reinterpret_cast<char *>(&network_update_frequency), sizeof(network_update_frequency));

    out.write(reinterpret_cast<char *>(&exec_time), sizeof(exec_time));

    out.close();
}

int main(int argc, char** argv)
{

    double Gei_mean = 0.25;
	double Gei_std = 0.0;
	double Gie_mean = 0.2;
    double Gie_std = 0.0;

    double trial_duration = 1000;
    double timestep = 0.02; // timestep for solving ODE for neurons

    double mean_soma = 0;
    double std_soma = 100;
    double mean_dend = 0;
    double std_dend = 200;
    
    double R = 1; // learning rate
    double SIDE = 100; // length of the square side of HVC

    double E_GABA_MATURE = -80; // GABA reverse potential of mature HVC(RA) neuron
    double E_GABA_IMMATURE = -55; // GABA reverse potential of newborn immature HVC(RA) neuron

    double A_P = 0.000025;
    double A_D = 0.000025;
    double F_0 = 1.2;

    double T_P = 10;
    double T_D = 10;
    double TAU_P = 30;
    double TAU_D = 30;

    double BETA = 0.9999;
    double BETA_SUPERSYNAPSE = 0.9999;
    
    double ACTIVATION = 0.0005;
    double SUPERSYNAPSE_THRESHOLD = 0.015;
    double G_MAX = 0.050;

    int RATE_WINDOW_SHORT = 50;
    int RATE_WINDOW_LONG = 500;

    double GABA_DOWN = 0.065;
    
    double GABA_RATE_THRESHOLD = 0.05;
    double MATURATION_RATE_THRESHOLD = 0.9;
    double DEATH_RATE_THRESHOLD = 0.01;

    double A_RA2I = 5.0;
    double SIGMA_RA2I = 8.0;
    double B_I2RA = 5.0;
    double SIGMA_I2RA = 6.0;

    double STDP_WINDOW = 100; // window size in which spikes contribute to LTP and LTD
    double MIN_INTERNEURON_DISTANCE = 0.01; // minimum allowed distance between neurons

    double WAITING_TIME = 100; // time in the beginning of the trial before current is injected to training neurons

    int N_RA, N_I;
	int N_ss = 4;
    int N_TR = 4;
    
    string outputDirectory;
    string outputFile;	
    int rank; // MPI process rank
    int np; // number of MPI processes running

    struct SynapticParameters syn_params; // structure with synaptic parameters
    struct SpatialConnectivityParameters spatial_params; // structure with spatial connectivity parameters
    struct GabaParameters gaba_params; // structure with gaba maturation paramters
	struct TimeParameters time_params; // structure with time parameters for simulation
    struct NoiseParameters noise_params; // structure with noise parameters for neurons

    // fill structure fields
    // synapse parameters structure
    syn_params.R = R;
	
    syn_params.T_P = T_P;
	syn_params.T_D = T_D;
	syn_params.TAU_P = TAU_P;
	syn_params.TAU_D = TAU_D;

    syn_params.A_P = A_P;
    syn_params.A_D = A_D;
	syn_params.F_0 = F_0;
    syn_params.BETA = BETA;
    syn_params.BETA_SUPERSYNAPSE = BETA_SUPERSYNAPSE;
        
    syn_params.ACTIVATION = ACTIVATION;
    syn_params.SUPERSYNAPSE_THRESHOLD = SUPERSYNAPSE_THRESHOLD;
    syn_params.G_MAX = G_MAX;
    
    syn_params.STDP_WINDOW = STDP_WINDOW;
    
    // spatial parameters structure
    spatial_params.SIDE = SIDE;
    
    spatial_params.Gei_mean = Gei_mean;
    spatial_params.Gei_std = Gei_std;
    spatial_params.Gie_mean = Gie_mean;
    spatial_params.Gie_std = Gie_std;

    spatial_params.A_RA2I = A_RA2I;
    spatial_params.SIGMA_RA2I = SIGMA_RA2I;
	spatial_params.B_I2RA = B_I2RA;
	spatial_params.SIGMA_I2RA = SIGMA_I2RA;
    
    spatial_params.MIN_INTERNEURON_DISTANCE = MIN_INTERNEURON_DISTANCE;

    // noise
    noise_params.mean_soma = mean_soma;
	noise_params.std_soma = std_soma;
	noise_params.mean_dend = mean_dend;
	noise_params.std_dend = std_dend;
    
    // maturation parameters structure
    gaba_params.E_GABA_IMMATURE = E_GABA_IMMATURE;
    gaba_params.E_GABA_MATURE = E_GABA_MATURE;

    gaba_params.GABA_RATE_THRESHOLD = GABA_RATE_THRESHOLD;
    gaba_params.MATURATION_RATE_THRESHOLD = MATURATION_RATE_THRESHOLD;
    gaba_params.DEATH_RATE_THRESHOLD = DEATH_RATE_THRESHOLD;
        
    gaba_params.RATE_WINDOW_SHORT = RATE_WINDOW_SHORT;
    gaba_params.RATE_WINDOW_LONG = RATE_WINDOW_LONG;
        
    gaba_params.GABA_DOWN = GABA_DOWN;
    
    // time parameters
    time_params.timestep = timestep;
    time_params.trial_duration = trial_duration;
    time_params.WAITING_TIME = WAITING_TIME;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    
    // parse command line parameters
    //
    if (argc > 1)
    {
        // number of neurons
        N_RA = atoi(argv[1]);
        N_I = atoi(argv[2]);
	
        time_params.network_update_frequency = atof(argv[3]);

        outputFile = argv[4];
        outputDirectory = argv[5];
        
        if (rank == 0)
            printf("Output directory is %s\n", outputDirectory.c_str());
    }
    else
    {
        printf("Command line parameters are needed!\n");
        return -1;
    }

	PoolParallel pool(N_TR, N_RA, N_I, N_ss, outputDirectory);

    pool.set_spatial_parameters(spatial_params);
    pool.set_synaptic_parameters(syn_params);
    pool.set_gaba_parameters(gaba_params);
    pool.set_noise_parameters(noise_params);
    pool.set_time_parameters(time_params);
    
    int num_layers = 10;

    pool.initialize_coordinates();
	pool.initialize_chain_connections(num_layers);

    //pool.print_invariable_connections();

	pool.print_simulation_parameters();
    
   	int save_freq_short = 20;
	int save_freq_long = 50;

    int num_trials = 10;

    double start, end;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    pool.run_trials_no_save(num_trials);
    
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();

    double average_execution_time;

    if (rank == 0)
    {
        average_execution_time = (end - start) / static_cast<double>(num_trials);

        std::cout << "Average execution time is " << average_execution_time << " seconds" << std::endl;
    
        write_execution_time(N_RA, N_I, np, timestep, time_params.network_update_frequency, average_execution_time, (outputDirectory + outputFile).c_str());
    }
	
	MPI_Finalize();


	return 0;

}
