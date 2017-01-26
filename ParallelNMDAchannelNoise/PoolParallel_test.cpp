#include "PoolParallel.h"
#include <mpi.h>
#include "structures.h"

using namespace std;

int main(int argc, char** argv)
{
	double Gei_mean = 0.25;
	double Gei_std = 0.0;
	double Gie_std = 0.0;

    double trial_duration = 1000;
    double timestep = 0.02; // timestep for solving ODE for neurons
    int trials = 100000;

    double network_update_frequency = 1.0; // frequency of network synchronization in ms
    double R = 1; // learning rate
    double SIDE = 100; // length of the square side of HVC

    double E_GABA_MATURE = -80; // GABA reverse potential of mature HVC(RA) neuron
    double E_GABA_IMMATURE = -55; // GABA reverse potential of newborn immature HVC(RA) neuron


    double T_P = 10;
    double T_D = 10;
    double TAU_P = 30;
    double TAU_D = 30;

    double STDP_WINDOW = 100; // window size in which spikes contribute to LTP and LTD
    double MIN_INTERNEURON_DISTANCE = 0.01; // minimum allowed distance between neurons

    double WAITING_TIME = 100; // time in the beginning of the trial before current is injected to training neurons

    int N_RA, N_I, N_ss, N_TR;

    string outputDirectory;
	
    int rank; // MPI process rank

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

    syn_params.STDP_WINDOW = STDP_WINDOW;
    
    // spatial parameters structure
    spatial_params.SIDE = SIDE;
    
    spatial_params.Gei_mean = Gei_mean;
    spatial_params.Gei_std = Gei_std;
    spatial_params.Gie_std = Gie_std;

    spatial_params.MIN_INTERNEURON_DISTANCE = MIN_INTERNEURON_DISTANCE;

    // maturation parameters structure
    gaba_params.E_GABA_IMMATURE = E_GABA_IMMATURE;
    gaba_params.E_GABA_MATURE = E_GABA_MATURE;

    // time parameters
    time_params.timestep = timestep;
    time_params.network_update_frequency = network_update_frequency;
    time_params.trial_duration = trial_duration;
    time_params.WAITING_TIME = WAITING_TIME;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // parse command line parameters
    //
    if (argc > 1)
    {
        // number of neurons
        N_TR = atoi(argv[1]);
        N_RA = atoi(argv[2]);
        N_I = atoi(argv[3]);
        N_ss = atoi(argv[4]);
	
        // spatial params
        spatial_params.Gie_mean = atof(argv[5]);
        
        spatial_params.A_RA2I = atof(argv[6]);
		spatial_params.SIGMA_RA2I = atof(argv[7]);
		spatial_params.B_I2RA = atof(argv[8]);
		spatial_params.SIGMA_I2RA = atof(argv[9]);
        
        // noise
        noise_params.mean_soma = atof(argv[10]);
		noise_params.std_soma = atof(argv[11]);
		noise_params.mean_dend = atof(argv[12]);
		noise_params.std_dend = atof(argv[13]);
		
        // synaptic params
        syn_params.A_P = atof(argv[14]);
        syn_params.A_D = atof(argv[15]);
		syn_params.F_0 = atof(argv[16]);
        syn_params.BETA = atof(argv[17]);
        syn_params.BETA_SUPERSYNAPSE = atof(argv[18]);
        syn_params.ACTIVATION = atof(argv[19]);
        syn_params.SUPERSYNAPSE_THRESHOLD = atof(argv[20]);
        syn_params.G_MAX = atof(argv[21]);
        
        // gaba params
        gaba_params.GABA_RATE_THRESHOLD = atof(argv[22]);
        gaba_params.MATURATION_RATE_THRESHOLD = atof(argv[23]);
        gaba_params.DEATH_RATE_THRESHOLD = atof(argv[24]);
        
        gaba_params.RATE_WINDOW_SHORT = atof(argv[25]);
        gaba_params.RATE_WINDOW_LONG = atof(argv[26]);
        
        gaba_params.GABA_DOWN = atof(argv[27]);
        
        outputDirectory = argv[28];
        
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
    
    pool.initialize_coordinates();
	pool.initialize_connections();

    //pool.print_invariable_connections();

	pool.print_simulation_parameters();
    
   	int save_freq_short = 20;
	int save_freq_long = 20;

    double start_time = MPI_Wtime();
    
    pool.chain_growth(save_freq_short, save_freq_long);
	
    double end_time = MPI_Wtime();

    if (rank == 0)
        printf("Execution time = %f\n", end_time - start_time);

	MPI_Finalize();


	return 0;

}
