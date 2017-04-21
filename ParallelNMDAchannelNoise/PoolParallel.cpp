#include "PoolParallel.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <functional>
#include "training_current.h"

using namespace std::placeholders;

PoolParallel::PoolParallel(const Configuration& cfg)
{

    MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
    
    // initialize internal time
    internal_time = 0.0;
    network_time = 0.0;
    trial_number = 0;

    this->initialize_generator();

    // set network parameters
    struct NetworkParameters network_params = cfg.get_network_parameters();
    
    this->set_network_parameters(network_params);

    // distibute work among processes
	N_RA_sizes.resize(MPI_size); // array with number of RA neurons per process
	N_I_sizes.resize(MPI_size); // array with number of I neurons per process

	for (int i = 0; i < MPI_size; i++)
	{
		N_RA_sizes[i] = N_RA / MPI_size;
		N_I_sizes[i] = N_I / MPI_size;
	}
	int RA_remain = N_RA % MPI_size;
	int I_remain = N_I % MPI_size;
	int j = 0;

	// distribute RA neurons
	while (RA_remain > 0)
	{
		N_RA_sizes[j] += 1;
		RA_remain -= 1;
		j += 1;

		if (j >= MPI_size)
			j -= MPI_size;
	}

	// distribute I neurons
	j = 0;
	while (I_remain > 0)
	{
		
		N_I_sizes[j] += 1;
		I_remain -= 1;
		j += 1;

		if (j >= MPI_size)
			j -= MPI_size;
	}

	N_RA_local = N_RA_sizes[MPI_rank]; // assign number of RA neurons for each process
	N_I_local = N_I_sizes[MPI_rank]; // assign number of I neurons for each process

	printf("My rank is %d; N_RA_local = %d; N_I_local = %d\n", MPI_rank, N_RA_local, N_I_local);

    // resize all data arrays
    this->resize_arrays_for_RA(N_RA_local, N_RA);
    this->resize_arrays_for_I(N_I_local, N_I);
	
    // assign real id to neurons
	for (int i = 0; i < N_RA_local; i++)
	{
		// assign real Id for RA neurons
		int N = 0; // number of neurons in the processes with lower rank
		
		for (int k = 0; k < MPI_rank; k++)
			N += N_RA_sizes[k];

		Id_RA_local[i] = N + i;
	}
	
    // assign real ID for I neurons
    for (int i = 0; i < N_I_local; i++)
	{
		int N = 0; // number of neurons in the processes with lower rank
		
		for (int k = 0; k < MPI_rank; k++)
			N += N_I_sizes[k];

		Id_I_local[i] = N + i;
	}

    // prepare global array with ids of all HVC(RA) neurons
    Id_RA_global.resize(N_RA); 
    

    for (int i = 0; i < N_RA; i++)
        Id_RA_global[i] = i;

    // set all parameters
    struct SynapticParameters synaptic_params = cfg.get_synaptic_parameters();
    struct SpatialParameters spatial_params = cfg.get_spatial_parameters();
    struct GabaParameters gaba_params = cfg.get_gaba_parameters();
    struct TimeParameters time_params = cfg.get_time_parameters();
    struct NoiseParameters noise_params = cfg.get_noise_parameters();
    
    outputDirectory = cfg.get_output_directory();
    
    this->set_synaptic_parameters(synaptic_params);
    this->set_spatial_parameters(spatial_params);
    this->set_gaba_parameters(gaba_params);
    this->set_time_parameters(time_params);
    this->set_noise_parameters(noise_params);

    this->write_global_index_array((outputDirectory + "global_index_array.bin").c_str());
}

void PoolParallel::set_network_parameters(const struct NetworkParameters& network_params)
{
    N_RA = network_params.N_RA;
    N_TR = network_params.N_TR;
    N_I = network_params.N_I;
}

void PoolParallel::set_spatial_parameters(const struct SpatialParameters& spatial_params)
{
    SIDE = spatial_params.SIDE;
    MIN_INTERNEURON_DISTANCE = spatial_params.MIN_INTERNEURON_DISTANCE;

    A_RA2I = spatial_params.A_RA2I;
    SIGMA_RA2I = spatial_params.SIGMA_RA2I;
    B_I2RA = spatial_params.B_I2RA;
    SIGMA_I2RA = spatial_params.SIGMA_I2RA;

    Gei_mean = spatial_params.Gei_mean;
    Gei_std = spatial_params.Gei_std;
    Gie_mean = spatial_params.Gie_mean;
    Gie_std = spatial_params.Gie_std;

}

void PoolParallel::set_synaptic_parameters(const struct SynapticParameters& syn_params)
{
    Nss = syn_params.Nss;
 
    R = syn_params.R;

    A_P = syn_params.A_P;
    A_D = syn_params.A_D;
    F_0 = syn_params.F_0;
    T_P = syn_params.T_P;
    T_D = syn_params.T_D;
    TAU_P = syn_params.TAU_P;
    TAU_D = syn_params.TAU_D;
    
    BETA = syn_params.BETA;
    BETA_SUPERSYNAPSE = syn_params.BETA_SUPERSYNAPSE;
    
    ACTIVATION_THRESHOLD = syn_params.ACTIVATION_THRESHOLD;
    SUPERSYNAPSE_THRESHOLD = syn_params.SUPERSYNAPSE_THRESHOLD;
    WEIGHT_MAX = syn_params.WEIGHT_MAX;
    
    STDP_WINDOW = syn_params.STDP_WINDOW;
}

void PoolParallel::set_noise_parameters(const struct NoiseParameters& noise_params)
{
    white_noise_mean_soma = noise_params.white_noise_mean_soma;
    white_noise_std_soma = noise_params.white_noise_std_soma;
    white_noise_mean_dend = noise_params.white_noise_mean_dend;
    white_noise_std_dend = noise_params.white_noise_std_dend;

	for (int i = 0; i < N_RA_local; i++)
	{	
		HVCRA_local[i].set_noise_generator(&generator);
		HVCRA_local[i].set_white_noise(white_noise_mean_soma, white_noise_std_soma, white_noise_mean_dend, white_noise_std_dend);
	}

    for (int i = 0; i < N_I_local; i++)
	{
        HVCI_local[i].set_noise_generator(&generator);
        HVCI_local[i].set_poisson_noise();
    }
}

void PoolParallel::set_time_parameters(const struct TimeParameters& time_params)
{
    timeStep = time_params.timestep;
    network_update_frequency = time_params.network_update_frequency;
    trial_duration = time_params.trial_duration;

    WAITING_TIME = time_params.WAITING_TIME;

    size = (int) round(trial_duration/timeStep) + 1;

    for (int i = 0; i < N_RA_local; i++)
		HVCRA_local[i].set_dynamics(trial_duration, timeStep);
	
    for (int i = 0; i < N_I_local; i++)
		HVCI_local[i].set_dynamics(trial_duration, timeStep);
}

void PoolParallel::set_gaba_parameters(const struct GabaParameters& gaba_params)
{
    E_GABA_MATURE = gaba_params.E_GABA_MATURE;
    E_GABA_IMMATURE = gaba_params.E_GABA_IMMATURE;

    GABA_RATE_THRESHOLD = gaba_params.GABA_RATE_THRESHOLD;
    MATURATION_RATE_THRESHOLD = gaba_params.MATURATION_RATE_THRESHOLD;
    DEATH_RATE_THRESHOLD = gaba_params.DEATH_RATE_THRESHOLD;
    
    RATE_WINDOW_SHORT = gaba_params.RATE_WINDOW_SHORT;
    RATE_WINDOW_LONG = gaba_params.RATE_WINDOW_LONG;
    
    GABA_DOWN = gaba_params.GABA_DOWN;

	std::fill(num_bursts_in_recent_trials.begin(), num_bursts_in_recent_trials.end(), intBuffer(RATE_WINDOW_LONG));

	for (size_t i = 0; i < num_bursts_in_recent_trials.size(); i++)
		for (int j = 0; j < RATE_WINDOW_LONG; j++)
			num_bursts_in_recent_trials[i].push_back(0);
	
    // set training neurons to be mature, all other neurons are immature
	for (int i = 0; i < N_RA_local; i++)
	{
		if (Id_RA_local[i] < N_TR)
		{
			mature_local[i] = 1;
			gaba_potential_local[i] = E_GABA_MATURE;
		}
		else
		{
			mature_local[i] = 0;
			gaba_potential_local[i] = E_GABA_IMMATURE;
		}
	}
}

void PoolParallel::initialize_generator()
{
    // prepare seeds and initialize each generator with different seed
    srand((unsigned int) time(NULL));

    //unsigned seed = rand();
    //unsigned seed = rand() + MPI_rank;

    unsigned seed = 15000 + MPI_rank*1000;

    generator.set_seed(seed);
}


void PoolParallel::initialize_coordinates_for_added_neurons(int n_total_old)
{
    if (MPI_rank == 0)
    {
        double xx; // temporary x-coordinate
        double yy; // temporary y-coordinate

        bool close; // are neurons too close or not

        // set coordinates for new added HVC(RA) neurons
        for (int i = n_total_old; i < N_RA; i++)
        {
            do
            {
                close = false;
                xx = 0.5*SIDE / (sqrt(N_I)+1) + generator.random(SIDE - SIDE / (sqrt(N_I)+1));
                yy = 0.5*SIDE / (sqrt(N_I)+1) + generator.random(SIDE - SIDE / (sqrt(N_I)+1));

                // check distances to all I neurons

                for (int j = 0; j < N_I; j++)
                {
                    if (distance(xx,yy,xx_I[j],yy_I[j]) < MIN_INTERNEURON_DISTANCE)
                    {
                        close = true;
                        break;
                    }
                }

                // check distances to all previous RA neurons
                if (!close)
                    for (int j = 0; j < i; j++)
                    {
                        if (distance(xx, yy, xx_RA[j], yy_RA[j]) < MIN_INTERNEURON_DISTANCE)
                    {
                        close = true;
                        break;
                    }

                }
            } while(close);

            xx_RA.push_back(xx);
            yy_RA.push_back(yy);
		}
        
        std::string fileRAxy = outputDirectory + "RA_xy_after_addition.bin";

        this->write_coordinates_RA(fileRAxy.c_str());
    }


}

void PoolParallel::initialize_coordinates_for_dispersed_training()
{
	if (MPI_rank == 0)
	{
		this->initialize_coordinates();
		
		// fix training neuron 0 and move other training neurons if necessary
		int neuron_id_to_check = 1; // id of the next neuron to check
        bool found_next_distant_neuron;

		for (int i = 1; i < N_TR; i++)
		{
			
			// calculate distances from next_neuron_to_check to previously confirmed distant training neurons
			do 
			{
			    found_next_distant_neuron = true;
				
                for (int j = 0; j < i; j++)
					if (distance(xx_RA[j], yy_RA[j], xx_RA[neuron_id_to_check], yy_RA[neuron_id_to_check]) < 5 * SIDE / (sqrt(N_I)+1))
					{
                        
                        std::cout << "Neuron " << neuron_id_to_check << std::endl;
                            
                        std::cout << "Distance to neuron " << j << " = " << distance(xx_RA[j], yy_RA[j], xx_RA[neuron_id_to_check], yy_RA[neuron_id_to_check]) 
                                << " is smaller than " << 5 * SIDE / (sqrt(N_I) + 1) << std::endl;

						neuron_id_to_check++;
						found_next_distant_neuron = false;
						break;
					}
			} while(!found_next_distant_neuron);
			
			double tempx = xx_RA[i];
			double tempy = yy_RA[i];
			
            std::cout << "Found neuron " << neuron_id_to_check << std::endl;

			xx_RA[i] = xx_RA[neuron_id_to_check];
			yy_RA[i] = yy_RA[neuron_id_to_check];
			
			xx_RA[neuron_id_to_check] = tempx;
			yy_RA[neuron_id_to_check] = tempy;
			
			neuron_id_to_check++;
			
		}
	}
}

void PoolParallel::generate_default_network(std::string networkDir)
{
	outputNetworkDir = networkDir;
	
    this->initialize_coordinates();
    this->write_all_coordinates();
    this->initialize_connections();    
    
    this->write_global_index_array((networkDir + "global_index_array.bin").c_str());
}

void PoolParallel::generate_network_with_clustered_training(std::string networkDir)
{
	outputNetworkDir = networkDir;
	
	// initialize coordinates with clustered training neurons and connections
    this->initialize_coordinates_for_clustered_training();
    this->write_all_coordinates();
    this->initialize_connections(); 
    this->write_global_index_array((networkDir + "global_index_array.bin").c_str());
}

void PoolParallel::generate_network_with_dispersed_training(std::string networkDir)
{
	outputNetworkDir = networkDir;
	
	// initialize coordinates with dispersed training neurons and connections
    this->initialize_coordinates_for_dispersed_training();
    this->write_all_coordinates();
    this->initialize_connections(); 
    this->write_global_index_array((networkDir + "global_index_array.bin").c_str());
}

void PoolParallel::initialize_coordinates_for_clustered_training()
{
    if (MPI_rank == 0)
    {
        this->initialize_coordinates();

        // now select neuron 0 as the first training neuron and find distances to all other neurons in the pool

        double xx = xx_RA[0];
        double yy = yy_RA[0];

        std::vector<double> distances_to_pool_neurons(xx_RA.size()-1); 

        for (size_t i = 1; i < xx_RA.size(); i++)
            distances_to_pool_neurons[i-1] = distance(xx, yy, xx_RA[i], yy_RA[i]);

        // sort distances
        std::vector<size_t> idx(distances_to_pool_neurons.size()); // vector with sorted indices of distances

        std::iota(idx.begin(), idx.end(), 0);

        std::sort(idx.begin(), idx.end(), [&distances_to_pool_neurons](size_t i1, size_t i2)
                                            {return distances_to_pool_neurons[i1] < distances_to_pool_neurons[i2];});

        // rearrange neuron coordinates: make N_TR - 1 neurons with smallest distances to have id from 1 to N_TR-1
        // temporary coordinates for swapping neurons
        double xx_temp;
        double yy_temp;

        for (int i = 0; i < N_TR-1; i++)
        {
            xx_temp = xx_RA[i+1];
            yy_temp = yy_RA[i+1];

            xx_RA[i+1] = xx_RA[idx[i]+1];
            yy_RA[i+1] = yy_RA[idx[i]+1];

            xx_RA[idx[i]+1] = xx_temp;
            yy_RA[idx[i]+1] = yy_temp;
        }
    }
}

void PoolParallel::initialize_coordinates()
{
    if (MPI_rank == 0)
    {
        double xx; // temporary x-coordinate
        double yy; // temporary y-coordinate

        bool close; // are neurons too close or not

        // set coordinates for HVC(I) neurons
        for (int i = 0; i < (int) sqrt(N_I); i++)
        {
            for (int k = 0; k < (int) sqrt(N_I); k++)
            {
                do
                {
					close = false;
					xx = (double) (i+1) * SIDE / (sqrt(N_I)+1) + generator.random(0.25 * SIDE / (sqrt(N_I) + 1));
					yy = (double) (k+1) * SIDE / (sqrt(N_I)+1) + generator.random(0.25 * SIDE / (sqrt(N_I) + 1));

					// check distances to all previous I neurons

                    for (size_t j = 0; j < xx_I.size(); j++)
                    {
                        if (distance(xx,yy,xx_I[j],yy_I[j]) < MIN_INTERNEURON_DISTANCE)
                        {
                            close = true;
                            break;
                        }
                    }
				} while(close);

			xx_I.push_back(xx);
			yy_I.push_back(yy);

            }
        }


        // set coordinates for HVC(RA) neurons
        for (int i = 0; i < N_RA; i++)
        {
            do
            {
                close = false;
                xx = 0.5*SIDE / (sqrt(N_I)+1) + generator.random(SIDE - SIDE / (sqrt(N_I)+1));
                yy = 0.5*SIDE / (sqrt(N_I)+1) + generator.random(SIDE - SIDE / (sqrt(N_I)+1));

                // check distances to all I neurons

                for (int j = 0; j < N_I; j++)
                {
                    if (distance(xx,yy,xx_I[j],yy_I[j]) < MIN_INTERNEURON_DISTANCE)
                    {
                        close = true;
                        break;
                    }
                }

                // check distances to all previous RA neurons
                if (!close)
                    for (int j = 0; j < i; j++)
                    {
                        if (distance(xx, yy, xx_RA[j], yy_RA[j]) < MIN_INTERNEURON_DISTANCE)
                    {
                        close = true;
                        break;
                    }

                }
            } while(close);

            xx_RA.push_back(xx);
            yy_RA.push_back(yy);
		}
		
		
    }
}

void PoolParallel::print_simulation_parameters()
{
	if (MPI_rank == 0)
	{
        printf("\nSimulation parameters: \n");

		printf("\nNetwork\n");

		printf("\nN_RA = %d\n", N_RA);
		printf("N_I = %d\n", N_I);
		printf("N_TR = %d\n", N_TR);
		
		printf("\nSpatial parameters\n");
		
        printf("\nA_RA2I = %f\n", A_RA2I);
		printf("SIGMA_RA2I = %f\n", SIGMA_RA2I);
		printf("B_I2RA = %f\n", B_I2RA);
		printf("SIGMA_I2RA = %f\n", SIGMA_I2RA);
		
        printf("\nGei_mean = %f\n", Gei_mean);
		printf("Gei_std = %f\n", Gei_std);
		printf("Gie_mean = %f\n", Gie_mean);
		printf("Gie_std = %f\n", Gie_std);

		printf("\nNoise:\n");
		
        printf("\nwhite_noise_mean_soma = %f\n", white_noise_mean_soma);
		printf("white_noise_std_soma = %f\n", white_noise_std_soma);
		printf("white_noise_mean_dend = %f\n", white_noise_mean_dend);
		printf("white_noise_std_dend = %f\n", white_noise_std_dend);

        
		printf("\nTime parameters\n");
		
        printf("\ntrial_duration = %f\n", trial_duration);
		printf("timestep = %f\n", timeStep);
		printf("network_update_frequency = %f\n", network_update_frequency);
		printf("WAITING_TIME = %f\n", WAITING_TIME);

		printf("\nSynaptic parameters\n");
        
        printf("\nNss = %d\n", Nss);
        
        printf("\nR = %f\n", R);
		
		printf("\nA_P = %f\n", A_P);
		printf("A_D = %f\n", A_D);
		printf("F_0 = %f\n", F_0);
		printf("T_P = %f\n", T_P);
		printf("TAU_P = %f\n", TAU_P);
		printf("T_D = %f\n", T_D);
		printf("TAU_D = %f\n", TAU_D);
        
        printf("\nBETA = %f\n", BETA);
		printf("BETA_SUPERSYNAPSE = %f\n", BETA_SUPERSYNAPSE);
		
        printf("\nACTIVATION_THRESHOLD = %f\n", ACTIVATION_THRESHOLD);
		printf("SUPERSYNAPSE_THRESHOLD = %f\n", SUPERSYNAPSE_THRESHOLD);
		printf("WEIGHT_MAX = %f\n", WEIGHT_MAX);
		
        printf("\nSTDP_WINDOW = %f\n", STDP_WINDOW);

		printf("\nGABA developmental switch\n");
		
        printf("\nE_GABA_MATURE = %f\n", E_GABA_MATURE);
        printf("E_GABA_IMMATURE = %f\n", E_GABA_IMMATURE);
        
        printf("\nGABA_DOWN = %f\n", GABA_DOWN);
		
        printf("\nGABA_RATE_THRESHOLD = %f\n", GABA_RATE_THRESHOLD);
        printf("MATURATION_RATE_THRESHOLD = %f\n", MATURATION_RATE_THRESHOLD);
        printf("DEATH_RATE_THRESHOLD = %f\n", DEATH_RATE_THRESHOLD);
		
        printf("\nRATE_WINDOW_SHORT = %d\n", RATE_WINDOW_SHORT);
        printf("RATE_WINDOW_LONG = %d\n", RATE_WINDOW_LONG);
	}
}

void PoolParallel::read_fixed_network(std::string networkDir)
{
	std::string fileRA2I = networkDir + "RA_I_connections.bin";
    std::string fileI2RA = networkDir + "I_RA_connections.bin";
   
    std::string fileGlobalIndexArray = networkDir + "global_index_array.bin";

    this->read_global_index_array(fileGlobalIndexArray.c_str());
    this->read_invariable_synapses(fileRA2I.c_str(), fileI2RA.c_str());
   
    this->send_connections();
    
    
}

void PoolParallel::read_network_state(std::string dirname, int starting_trial)
{
	std::string trial_extension; // extension for file name, which takes into account trial number
	
	if (starting_trial > 0)
		trial_extension = "_" + std::to_string(starting_trial) + "_";
	else
		trial_extension = "";
		
	
    std::string fileRA2I = dirname + "RA_I_connections.bin";
    std::string fileI2RA = dirname + "I_RA_connections.bin";
   
    std::string fileRA2INew = outputDirectory + "RA_I_connections_NEW.bin";
    std::string fileI2RANew = outputDirectory + "I_RA_connections_NEW.bin";

    std::string fileActiveGraph = dirname + "RA_RA_active_connections" + trial_extension + ".bin";
    std::string fileActiveGraphNew = outputDirectory + "RA_RA_active_connections" + trial_extension + "NEW.bin";
    
    std::string fileSuperGraph = dirname + "RA_RA_super_connections" + trial_extension + ".bin";
    std::string fileSuperGraphNew = outputDirectory + "RA_RA_super_connections" + trial_extension + "NEW.bin";
    
    std::string fileMaturationGraph = dirname + "mature" + trial_extension + ".bin";
    std::string fileMaturationGraphNew = outputDirectory + "mature" + trial_extension + "NEW.bin";

    std::string fileWeightsGraph = dirname + "weights"  + trial_extension + ".bin";
    std::string fileWeightsGraphNew = outputDirectory + "weights"  + trial_extension + "NEW.bin";

    std::string fileGlobalIndexArray = outputDirectory + "global_index_array.bin";

    this->read_global_index_array(fileGlobalIndexArray.c_str());
    this->read_super_synapses(fileSuperGraph.c_str());
    this->read_active_synapses(fileActiveGraph.c_str());
    this->read_weights(fileWeightsGraph.c_str());
    this->read_maturation_info(fileMaturationGraph.c_str());
    this->read_invariable_synapses(fileRA2I.c_str(), fileI2RA.c_str());
    //this->print_super();
   
    trial_number++;

    this->send_connections();
    
    this->gather_data();
    
    this->write_supersynapses(fileSuperGraphNew.c_str());
    this->write_active_synapses(fileActiveGraphNew.c_str());
    this->write_weights(fileWeightsGraphNew.c_str());
    this->write_maturation_info(fileMaturationGraphNew.c_str());
    this->write_invariable_synapses(fileRA2INew.c_str(), fileI2RANew.c_str());
    
}

void PoolParallel::read_invariable_synapses(const char* RA_I, const char* I_RA)
{
    // read all connections to master process 0
    if (MPI_rank == 0)
    {
        std::ifstream inp_RA_I, inp_I_RA;
        
        // open files
        inp_RA_I.open(RA_I, std::ios::binary | std::ios::in);
        inp_I_RA.open(I_RA, std::ios::binary | std::ios::in);
    
        int N_ra;
        int N_i;

        inp_RA_I.read(reinterpret_cast<char *>(&N_ra), sizeof(N_ra));
        inp_I_RA.read(reinterpret_cast<char *>(&N_i), sizeof(N_i));

        if (N_ra != N_RA)
            std::cerr << "Number of HVC(RA) neurons read from file with invariable synapses: N_ra = " << N_ra << "is different from N_RA = " << N_RA << std::endl;
        
        if (N_i != N_I)
            std::cerr << "Number of HVC(I) neurons read from file with invariable synapses: N_i = " << N_i << "is different from N_I = " << N_I << std::endl;
        
        // read connections from RA to I
        for (int i = 0; i < N_RA; i++)
        {
            int n_id; // neuronal id
            int size; // number of outgoing connections

            inp_RA_I.read(reinterpret_cast<char *>(&n_id), sizeof(n_id));
            inp_RA_I.read(reinterpret_cast<char *>(&size), sizeof(size)); // write neuron's ID

            syn_ID_RA_I_global[i].resize(size);
            weights_RA_I_global[i].resize(size);

            for (int j = 0; j < size; j++)
            {
                inp_RA_I.read(reinterpret_cast<char *>(&syn_ID_RA_I_global[i][j]), sizeof(syn_ID_RA_I_global[i][j]));
                inp_RA_I.read(reinterpret_cast<char *>(&weights_RA_I_global[i][j]), sizeof(weights_RA_I_global[i][j]));

            }

        }

        // read connections from I to RA
        for (int i = 0; i < N_I; i++)
        {
            int n_id; // neuronal id
            int size; // number of outgoing connections

            inp_I_RA.read(reinterpret_cast<char *>(&n_id), sizeof(n_id));
            inp_I_RA.read(reinterpret_cast<char *>(&size), sizeof(size)); // write neuron's ID
            
            syn_ID_I_RA_global[i].resize(size);
            weights_I_RA_global[i].resize(size);
            
            for (int j = 0; j < size; j++)
            {
                inp_I_RA.read(reinterpret_cast<char *>(&syn_ID_I_RA_global[i][j]), sizeof(syn_ID_I_RA_global[i][j])); // write targets ID
                inp_I_RA.read(reinterpret_cast<char *>(&weights_I_RA_global[i][j]), sizeof(weights_I_RA_global[i][j])); // write targets conductance

            }
        }
        // close files
        inp_I_RA.close();
        inp_RA_I.close();
	}
}

void PoolParallel::read_global_index_array(const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ifstream inp;
      
        inp.open(filename, std::ios::in | std::ios::binary );

        // read number of HVC(RA) neurons
      int N;
           
        inp.read(reinterpret_cast<char *>(&N), sizeof(N));

        for (size_t i = 0; i < Id_RA_global.size(); i++)
            inp.read(reinterpret_cast<char *>(&Id_RA_global[i]), sizeof(Id_RA_global[i]));
		
        inp.close();
	}
}

void PoolParallel::read_maturation_info(const char* filename)
{
    std::ifstream inp;

    inp.open(filename, std::ios::in | std::ios::binary);

    // read number of HVC(RA) neurons
  int N;

    inp.read(reinterpret_cast<char *>(&N), sizeof(N));

    if (N != N_RA)
        std::cerr << "Number of HVC(RA) neurons read from file with active synapses: N = " << N << "is different from N_RA = " << N_RA << std::endl;

    inp.read(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));

    if (MPI_rank == 0)
        std::cout << "Trial number read from maturation file: trial_number = " << trial_number << std::endl;
    
    int counter = 0; // counter of neuron id in the process
   
    for (int i = 0; i < N_RA; i++)
    {
        bool neuron_found = false; // indicator that we found a neuron

        if (counter < N_RA_local)
            if (Id_RA_local[counter] == i)
            {
                neuron_found = true;
            }

        double gaba; // reverse GABA potential
        double firing_rate; // neuronal firing rate
        int remodeled; // indicator of neuronal axon remodeling
        int mature; // indicator of neuronal maturation
        
        inp.read(reinterpret_cast<char *>(&gaba), sizeof(gaba));
        inp.read(reinterpret_cast<char *>(&firing_rate), sizeof(firing_rate));
        inp.read(reinterpret_cast<char *>(&remodeled), sizeof(remodeled));
        inp.read(reinterpret_cast<char *>(&mature), sizeof(mature));

        if (neuron_found)
        {
            gaba_potential_local[counter] = gaba;
            firing_rate_local[counter] = firing_rate;
            remodeled_local[counter] = remodeled;
            mature_local[counter] = mature;
        }

        // increase counter
        if (neuron_found)
            counter++;

    }
    inp.close();
}

void PoolParallel::read_active_synapses(const char *filename)
{
    std::ifstream inp;

    inp.open(filename, std::ios::in | std::ios::binary);

    // read number of HVC(RA) neurons
  int N;

    inp.read(reinterpret_cast<char *>(&N), sizeof(N));

    if (N != N_RA)
        std::cerr << "Number of HVC(RA) neurons read from file with active synapses: N = " << N << "is different from N_RA = " << N_RA << std::endl;

    int counter = 0; // counter of neuron id in the process

    for (int i = 0; i < N_RA; i++)
    {
        // read neuron's ID number
        int n_id;
        
        inp.read(reinterpret_cast<char *>(&n_id), sizeof(n_id));
        
        // read number of connections to RA neurons
      int num_active;
        bool source_neuron_found = false; // indicator that we found a source neuron

        inp.read(reinterpret_cast<char *>(&num_active), sizeof(num_active));
        
        if (counter < N_RA_local)
            if (Id_RA_local[counter] == i)
            {
                active_synapses_local[counter].resize(num_active);
                source_neuron_found = true;
            }

        // read active synapses
        for (int j = 0; j < num_active; j++)
        {
            int target_id; // id of active target
            double synapse_weight; // weight of synapse

            inp.read(reinterpret_cast<char *>(&target_id), sizeof(target_id));
            inp.read(reinterpret_cast<char *>(&synapse_weight), sizeof(synapse_weight));

            if (source_neuron_found)
            {
                active_synapses_local[counter][j] = target_id;
                active_indicators_local[counter][target_id] = 1;
            }

        }
        // increase counter
        if (source_neuron_found)
            counter++;

    }
    inp.close();
}

void PoolParallel::print_active()
{
    for (int i = 0; i < N_RA_local; i++)
        for (int j = 0; j < N_RA; j++)
            if (active_indicators_local[i][j] == 1)
                std::cout << "Active synapse " << Id_RA_local[i] << " -> " << j << std::endl; 
}

void PoolParallel::print_super()
{
    for (int i = 0; i < N_RA_local; i++)
        for (int j = 0; j < N_RA; j++)
            if (supersynapses_indicators_local[i][j] == 1)
                std::cout << "Active supersynapse " << Id_RA_local[i] << " -> " << j << std::endl; 
}

void PoolParallel::read_weights(const char* filename)
{
    std::ifstream inp;

    inp.open(filename, std::ios::in | std::ios::binary);

    // read number of HVC(RA) neurons
  int N;

    inp.read(reinterpret_cast<char *>(&N), sizeof(N));

    if (N != N_RA)
        std::cerr << "Number of HVC(RA) neurons read from file with supersynapses: N = " << N << "is different from N_RA = " << N_RA << std::endl;

    int num_trial;
    
    inp.read(reinterpret_cast<char *>(&num_trial), sizeof(num_trial));
    
    if (MPI_rank == 0)
        std::cout << "Trial number read from file with synaptic weights: " << num_trial << std::endl;

    int counter = 0; // counter of neuron id in the process



    for (int i = 0; i < N_RA; i++)
    {
        bool source_neuron_found = false; // indicator that source neuron is found
        
        if (counter < N_RA_local)
            if (Id_RA_local[counter] == i)
            {
              source_neuron_found = true;
            }
       
        for (int j = 0; j < N_RA; j++)
        {
            
            double synapse_weight; // weight of synapse

            inp.read(reinterpret_cast<char *>(&synapse_weight), sizeof(synapse_weight));

            if (source_neuron_found)
            {
                weights_local[counter][j] = synapse_weight;
            }

        }
        // increase counter
        if (source_neuron_found)
            counter++;

    }
    inp.close();
}

void PoolParallel::read_super_synapses(const char *filename)
{
    std::ifstream inp;

    inp.open(filename, std::ios::in | std::ios::binary);

    // read number of HVC(RA) neurons
  int N;

    inp.read(reinterpret_cast<char *>(&N), sizeof(N));

    if (N != N_RA)
        std::cerr << "Number of HVC(RA) neurons read from file with supersynapses: N = " << N << "is different from N_RA = " << N_RA << std::endl;

    int counter = 0; // counter of neuron id in the process

    for (int i = 0; i < N_RA; i++)
    {
        // read neuron's ID number
        int n_id;
        bool source_neuron_found = false; // indicator that source neuron is found
        
        inp.read(reinterpret_cast<char *>(&n_id), sizeof(n_id));
        
        // read number of connections to RA neurons
      int num_super;
        
        inp.read(reinterpret_cast<char *>(&num_super), sizeof(num_super));
        
        if (counter < N_RA_local)
            if (Id_RA_local[counter] == i)
            {
                supersynapses_local[counter].resize(num_super);
              source_neuron_found = true;
            }
        // read supersynapses
        for (int j = 0; j < num_super; j++)
        {
            int target_id; // id of supersynaptic target
            double synapse_weight; // weight of synapse

            inp.read(reinterpret_cast<char *>(&target_id), sizeof(target_id));
            inp.read(reinterpret_cast<char *>(&synapse_weight), sizeof(synapse_weight));

            if (source_neuron_found)
            {
                supersynapses_local[counter][j] = target_id;
                supersynapses_indicators_local[counter][target_id] = 1;
            }

        }
        // increase counter
        if (source_neuron_found)
            counter++;

    }
    inp.close();
}

void PoolParallel::send_connections()
{
	MPI_Status status;
	int* sendcounts_syn_num_RA;
	int* sendcounts_syn_num_I;
	int* displs_syn_num_RA;
	int* displs_syn_num_I;
	int* syn_num_RA;
	int* syn_num_I;
	int* syn_num_RA_local = new int[N_RA_local];
	int* syn_num_I_local = new int[N_I_local];

	if (MPI_rank == 0)
	{
		sendcounts_syn_num_RA = new int[MPI_size];
		sendcounts_syn_num_I = new int[MPI_size];
		syn_num_RA = new int[N_RA];
		syn_num_I = new int[N_I];

		displs_syn_num_RA = new int[MPI_size];
		displs_syn_num_I = new int[MPI_size];

		// send connections to all processes

		sendcounts_syn_num_RA[0] = N_RA_sizes[0];
		sendcounts_syn_num_I[0] = N_I_sizes[0];
		displs_syn_num_RA[0] = 0;
		displs_syn_num_I[0] = 0;

		for (int i = 1; i < MPI_size; i++)
		{
			sendcounts_syn_num_RA[i] = N_RA_sizes[i];
			sendcounts_syn_num_I[i] = N_I_sizes[i];
			displs_syn_num_RA[i] = displs_syn_num_RA[i-1] + sendcounts_syn_num_RA[i-1];
			displs_syn_num_I[i] = displs_syn_num_I[i-1] + sendcounts_syn_num_I[i-1];
		}

		for (int i = 0; i < N_RA; i++)
		{
			syn_num_RA[i] = static_cast<int>(syn_ID_RA_I_global[i].size());
			//printf("Master process; syn_num_RA[%d] = %d\n", i, syn_num_RA[i]);
		}

		for (int i = 0; i < N_I; i++)
		{
			syn_num_I[i] = static_cast<int>(syn_ID_I_RA_global[i].size());

			//printf("Master process; syn_num_I[%d] = %d\n", i, syn_num_I[i]);
		}


    }

	// send number of connections for each neuron

	MPI_Scatterv(&syn_num_RA[0], sendcounts_syn_num_RA, displs_syn_num_RA, MPI_INT,
		&syn_num_RA_local[0], N_RA_local, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Scatterv(&syn_num_I[0], sendcounts_syn_num_I, displs_syn_num_I, MPI_INT,
		&syn_num_I_local[0], N_I_local, MPI_INT, 0, MPI_COMM_WORLD);

	//for (int i = 0; i < N_RA_local; i++)
	//	printf("My rank = %d; sun_num_RA[%d] = %d\n", MPI_rank, Id_RA_local[i], syn_num_RA_local[i]);
	
    //for (int i = 0; i < N_I_local; i++)
	//	printf("My rank = %d; sun_num_I[%d] = %d\n", MPI_rank, Id_I_local[i], syn_num_I_local[i]);

	for (int i = 0; i < N_RA_local; i++)
	{
		syn_ID_RA_I_local[i].resize(syn_num_RA_local[i]);
		weights_RA_I_local[i].resize(syn_num_RA_local[i]);
	}

	for (int i = 0; i < N_I_local; i++)
	{
		syn_ID_I_RA_local[i].resize(syn_num_I_local[i]);
		weights_I_RA_local[i].resize(syn_num_I_local[i]);
	}

	// send connection ID and weights

	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_RA_local; i++)
		{
			weights_RA_I_local[i] = weights_RA_I_global[i];
			syn_ID_RA_I_local[i] = syn_ID_RA_I_global[i];
		}

		for (int i = 0; i < N_I_local; i++)
		{
			weights_I_RA_local[i] = weights_I_RA_global[i];
			syn_ID_I_RA_local[i] = syn_ID_I_RA_global[i];

		}

		int offset_RA = N_RA_local;
		int offset_I = N_I_local;

		for (int i = 1; i < MPI_size; i++)
		{
            //  send ID and weights of RA targets
			for (int j = 0; j < sendcounts_syn_num_RA[i]; j++)
			{
				//printf("Master. Here!\n");
				MPI_Send(&syn_ID_RA_I_global[offset_RA+j][0], syn_num_RA[offset_RA+j], MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&weights_RA_I_global[offset_RA+j][0], syn_num_RA[offset_RA+j], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				//for (int k = 0; k < syn_num_RA[offset_RA+j]; k++)
				//	printf("Master. RA neuron %d; RA2I %d; w = %f\n", offset_RA+j, syn_ID_RA_I_global[offset_RA+j][k],
                 //       weights_RA_I_global[offset_RA+j][k]);

			}
			offset_RA += sendcounts_syn_num_RA[i];

            // send ID and weights of I targets
			for (int j = 0; j < sendcounts_syn_num_I[i]; j++)
			{
				//printf("Master. Here!\n");
				MPI_Send(&syn_ID_I_RA_global[offset_I+j][0], syn_num_I[offset_I+j], MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&weights_I_RA_global[offset_I+j][0], syn_num_I[offset_I+j], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

				//for (int k = 0; k < syn_num_I[offset_I+j]; k++)
				//	printf("Master. I neuron %d; I2RA %d; w = %f\n", offset_I+j, syn_ID_I_RA_global[offset_I+j][k],
                  //      weights_I_RA_global[offset_I+j][k]);

			}
			offset_I += sendcounts_syn_num_I[i];
			//printf("Master. Number of synapses;\n");

		}

	}
	else
	{
        // receive ID and weights of RA targets
		for (int i = 0; i < N_RA_local; i++)
		{
			//printf("Rank = %d. Here!\n", MPI_rank);
			MPI_Recv(&syn_ID_RA_I_local[i][0], syn_num_RA_local[i], MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&weights_RA_I_local[i][0], syn_num_RA_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

			//for (int j = 0; j < syn_num_RA_local[i]; j++)
			//	printf("My rank = %d; RA neuron %d; RA2I %d; w = %f\n", MPI_rank, Id_RA_local[i], syn_ID_RA_I_local[i][j],
              //      weights_RA_I_local[i][j]);
        }

        // receive ID and weights of I targets
        for (int i = 0; i < N_I_local; i++)
        {
            MPI_Recv(&syn_ID_I_RA_local[i][0], syn_num_I_local[i], MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&weights_I_RA_local[i][0], syn_num_I_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

			//for (int j = 0; j < syn_num_I_local[i]; j++)
			//	printf("My rank = %d; I neuron %d; I2RA %d; w = %f\n", MPI_rank, Id_I_local[i], syn_ID_I_RA_local[i][j],
              //      weights_I_RA_local[i][j]);

		}
	}

	if (MPI_rank == 0)
	{
        delete [] sendcounts_syn_num_RA;
        delete [] sendcounts_syn_num_I;
        delete [] syn_num_RA;
        delete [] syn_num_I;
        delete [] syn_num_RA_local;
        delete [] syn_num_I_local;
        delete [] displs_syn_num_RA;
        delete [] displs_syn_num_I;
    }
}

void PoolParallel::disable_RA2I_connections()
{
	for (int i = 0; i < N_RA_local; i++)
	{
		weights_RA_I_local[i].clear();
		syn_ID_RA_I_local[i].clear();
	}
}

void PoolParallel::initialize_ideal_chain_connections(int num_layers)
{
    if (MPI_rank == 0)
    {
        // connections from HVC(RA) to HVC(I) and vice versa
        
        for (int i = 0; i < num_layers; i++)
        {
            // from HVC(RA) to HVC(I)
            for (int j = 0; j < Nss; j++)
            {
                double G = this->sample_Ge2i();

                weights_RA_I_global[i*Nss + j].push_back(G);
                syn_ID_RA_I_global[i*Nss + j].push_back(i*Nss + j);
            }
        
            // from HVC(I) to HVC(RA)
            for (int j = 0; j < Nss; j++)
            {
                double G = this->sample_Gi2e();

                weights_I_RA_global[i*Nss + j].push_back(G);
                syn_ID_I_RA_global[i*Nss + j].push_back((i+1)*Nss + j);

            }
        }
    }

    this->send_connections();

    std::string fileRA2I = outputDirectory + "RA_I_connections.bin";
    std::string fileI2RA = outputDirectory + "I_RA_connections.bin";
    std::string filePajekFixed = outputDirectory + "fixed.net";

    this->write_invariable_synapses(fileRA2I.c_str(), fileI2RA.c_str());
    this->write_pajek_fixed(filePajekFixed.c_str());
    
    // initialize excitatory connections for all chain layers
    for (int i = 0; i < num_layers; i++)
    {
        for (int j = 0; j < N_RA_local; j++)
        {
            if ( (Id_RA_local[j] >= Nss * i) && (Id_RA_local[j] < Nss * (i + 1)) )
            {
                int diff = Nss - (Id_RA_local[j] - Nss * i);

                for (int k = 0; k < Nss; k ++)
                    weights_local[j][Id_RA_local[j] + diff + k] = WEIGHT_MAX;
            }
        }
    }
    

    this->update_all_synapses();
}

void PoolParallel::initialize_random_chain_connections(int num_layers)
{
    std::vector<std::vector<int>> chain; // neurons in the chain

    chain.resize(num_layers);
    chain[0].resize(N_TR);

    for (int i = 0; i < N_TR; i++)
        chain[0][i] = i;

    for (int i = 1; i < num_layers; i++)
        chain[i].resize(Nss);

    if (MPI_rank == 0)
    {

        // connections for HVC(RA) neurons
		for (int i = 0; i < N_RA; i++)
     	{
	 		for (int j = 0; j < N_I; j++)
	        {
	         	 if (generator.random(1) < p_RA2I(i,j))
	             {
		        	 double G = this->sample_Ge2i();

		             weights_RA_I_global[i].push_back(G);
		             syn_ID_RA_I_global[i].push_back(j);
		         }
			 }

		 }
		// connections for HVC(I) neurons

		for (int i = 0; i < N_I; i++)
     	{
	 		for (int j = 0; j < N_RA; j++)
	        {
	         	 if (generator.random(1) < p_I2RA(i,j))
	             {
		        	 double G = this->sample_Gi2e();

		             weights_I_RA_global[i].push_back(G);
		             syn_ID_I_RA_global[i].push_back(j);
		         }
			 }
		 }
        
        // chain connections
        std::vector<int> neuronsInChain; // id of HVC(RA) neurons in chain

        for (int i = 0; i < N_TR; i++)
            neuronsInChain.push_back(i);

        // recruit layers
        for (int i = 1; i < num_layers; i++)
        {
            // recruit neurons in the layer
            for (int j = 0; j < Nss; j++)
            {
                bool neuronAlreadyInChain = false; // indicator that selected neuron is already in the chain

                int target_id; // id of the recruited target
            
                do
                {
                    // sample target id
                    target_id = generator.sample_integer(N_TR, N_RA);
                    
                    // check if target is already in the chain
                    std::vector<int>::iterator pos = std::find(neuronsInChain.begin(), neuronsInChain.end(), target_id);
                    
                    if (pos != neuronsInChain.end())
                        neuronAlreadyInChain = true;
                    else
                        neuronAlreadyInChain = false;

                } while (neuronAlreadyInChain);

                chain[i][j] = target_id;
                neuronsInChain.push_back(target_id);
            }
        }

        // print chain
        std::cout << "Chain:\n" << std::endl;
        for (int i = 0; i < num_layers; i++)
        {
            for (size_t j = 0; j < chain[i].size(); j++)
                std::cout << chain[i][j] << "\t";
            std::cout << std::endl;
        }
	}

    this->send_connections();

    std::string fileRA2I = outputDirectory + "RA_I_connections.bin";
    std::string fileI2RA = outputDirectory + "I_RA_connections.bin";
    std::string filePajekFixed = outputDirectory + "fixed.net";

    this->write_invariable_synapses(fileRA2I.c_str(), fileI2RA.c_str());
    this->write_pajek_fixed(filePajekFixed.c_str());

    // send neurons in the chain
    for (int i = 1; i < num_layers; i++)
        MPI_Bcast(&chain[i][0], Nss, MPI_INT, 0, MPI_COMM_WORLD);

    // connect neurons according to the chain

    for (int i = 0; i < num_layers-1; i++)
    {
        // all neurons in the source layer
        for (size_t j = 0; j < chain[i].size(); j++)
        {
            // connect to all neurons in the next layer
            for (size_t k = 0; k < chain[i+1].size(); k++)
            {
                // determine the location of the source neuron
                int rank, shift;
                
                this->get_neuronRA_location(chain[i][j], &rank, &shift);
                

                if (MPI_rank == rank)
                {
                    weights_local[shift][chain[i+1][k]] = WEIGHT_MAX;
            
                    std::cout << "Rank = " << MPI_rank << " shift = " << shift << " Id_RA_local = " << Id_RA_local[shift] << " source_id = " << chain[i][j] 
                              << std::endl;
                }
            }
        }
    }

    this->update_all_synapses();
    this->gather_data();
    
    std::string fileActiveGraph = outputDirectory + "RA_RA_active_connections.bin";
    std::string fileSuperGraph = outputDirectory + "RA_RA_super_connections.bin";
	    	
    this->write_supersynapses(fileSuperGraph.c_str());
    this->write_active_synapses(fileActiveGraph.c_str());
}

void PoolParallel::initialize_test_connections(int num_RA_targets, int num_RA_target_groups)
{

    if (MPI_rank == 0)
    {

        if (num_RA_targets > N_RA - N_TR)
            std::cerr << "Number of RA targets exceeds number of pool neurons!" << std::endl;

        // connect first training neuron to a single HVC(I) neuron
        weights_RA_I_global[0].push_back(Gei_mean);
        syn_ID_RA_I_global[0].push_back(0);

		// connect HVC(I) neuron to 
      double inhibition_strength = 0.6;  
            
        for (int j = N_TR; j < num_RA_targets; j++)
        {

             if ((j - N_TR) % num_RA_target_groups == 0)
                inhibition_strength += Gie_mean;
                
             weights_I_RA_global[0].push_back(inhibition_strength);
             syn_ID_I_RA_global[0].push_back(j);
        }
	}

    this->send_connections();

    std::string fileRA2I = outputDirectory + "RA_I_connections.bin";
    std::string fileI2RA = outputDirectory + "I_RA_connections.bin";
    std::string filePajekFixed = outputDirectory + "fixed.net";

    this->write_invariable_synapses(fileRA2I.c_str(), fileI2RA.c_str());
    this->write_pajek_fixed(filePajekFixed.c_str());

}

void PoolParallel::initialize_connections_for_added_neurons(int n_total_old)
{
    if (MPI_rank == 0)
    {

        // connections for new HVC(RA) neurons
		for (int i = n_total_old; i < N_RA; i++)
     	{
	 		for (int j = 0; j < N_I; j++)
	        {
	         	 if (generator.random(1) < p_RA2I(i,j))
	             {
		        	 double G = this->sample_Ge2i();

		             weights_RA_I_global[i].push_back(G);
		             syn_ID_RA_I_global[i].push_back(j);
		         }
			 }

		 }
		// connections for HVC(I) neurons to new HVC(RA) neurons

		for (int i = 0; i < N_I; i++)
     	{
	 		for (int j = n_total_old; j < N_RA; j++)
	        {
	         	 if (generator.random(1) < p_I2RA(i,j))
	             {
		        	 double G = this->sample_Gi2e();

		             weights_I_RA_global[i].push_back(G);
		             syn_ID_I_RA_global[i].push_back(j);
		         }
			 }
		 }
	}

    this->send_connections();

    std::string fileRA2I = outputDirectory + "RA_I_connections_after_addition.bin";
    std::string fileI2RA = outputDirectory + "I_RA_connections_after_addition.bin";
    std::string filePajekFixed = outputDirectory + "fixed_after_addition.net";

    this->write_invariable_synapses(fileRA2I.c_str(), fileI2RA.c_str());
    this->write_pajek_fixed(filePajekFixed.c_str());
}


void PoolParallel::initialize_connections()
{

    if (MPI_rank == 0)
    {

        // connections for HVC(RA) neurons
		for (int i = 0; i < N_RA; i++)
     	{
	 		for (int j = 0; j < N_I; j++)
	        {
	         	 if (generator.random(1) < p_RA2I(i,j))
	             {
		        	 double G = this->sample_Ge2i();

		             weights_RA_I_global[i].push_back(G);
		             syn_ID_RA_I_global[i].push_back(j);
		         }
			 }

		 }
		// connections for HVC(I) neurons

		for (int i = 0; i < N_I; i++)
     	{
	 		for (int j = 0; j < N_RA; j++)
	        {
	         	 if (generator.random(1) < p_I2RA(i,j))
	             {
		        	 double G = this->sample_Gi2e();

		             weights_I_RA_global[i].push_back(G);
		             syn_ID_I_RA_global[i].push_back(j);
		         }
			 }
		 }
	}

    this->send_connections();

    std::string fileRA2I = outputNetworkDir + "RA_I_connections.bin";
    std::string fileI2RA = outputNetworkDir + "I_RA_connections.bin";
    std::string filePajekFixed = outputNetworkDir + "fixed.net";

    this->write_invariable_synapses(fileRA2I.c_str(), fileI2RA.c_str());
    this->write_pajek_fixed(filePajekFixed.c_str());

}

void PoolParallel::print_invariable_connections()
{
    if (MPI_rank == 0)
    {
        for (int i = 0; i < N_RA; i++)
            std::cout << "RA neuron " << i << "has " << weights_RA_I_global[i].size() << " connections to I neurons\n" << std::endl;

        for (int i = 0; i < N_I; i++)
            std::cout << "I neuron " << i << "has " << weights_I_RA_global[i].size() << " connections to RA neurons" << std::endl;
    }
}

void PoolParallel::print_received_invariable_connections()
{
    /*
    for (int i = 0; i < N_RA_local; i++)
    {
        std::cout << "RA neuron " << Id_RA_local[i] << "has " << weights_RA_I_local[i].size() << " connections to I neurons:\n";
        
        for (size_t j = 0; j < weights_RA_I_local[i].size(); j++)
            std::cout << syn_ID_RA_I_local[i][j] << "\t";

        std:: cout << std::endl;
    }
    */
    
    for (int i = 0; i < N_I_local; i++)
    {
        std::cout << "I neuron " << Id_I_local[i] << "has " << weights_I_RA_local[i].size() << " connections to RA neurons:\n";
        
        for (size_t j = 0; j < weights_I_RA_local[i].size(); j++)
            std::cout << syn_ID_I_RA_local[i][j] << "\t";

        std:: cout << std::endl;
    }
}

void PoolParallel::set_all_mature()
{
    for (int i = 0; i < N_RA_local; i++)
        mature_local[i] = 1;

}

void PoolParallel::randomize_after_trial()
{
    std::fill(spike_times_dend_global.begin(), spike_times_dend_global.end(), -200.0);
    
    for (int i = 0; i < N_RA; i++)
    {
        spikes_in_trial_soma_global[i].clear();
        spikes_in_trial_dend_global[i].clear();

    }

    for (int i = 0; i < N_I; i++)
        spikes_in_trial_interneuron_global[i].clear();

    for (int i = 0; i < N_RA_local; i++)
    {
	    HVCRA_local[i].set_to_rest();

        spikes_in_trial_soma_local[i].clear();
        spikes_in_trial_dend_local[i].clear();
    }

    for (int i = 0; i < N_I_local; i++)
    {
        HVCI_local[i].set_to_rest();
        spikes_in_trial_interneuron_local[i].clear();
    }
}


void PoolParallel::set_training_current(double t)
{
    if (MPI_rank == 0)
    {
        std::function<double (double)> f = std::bind(&training_current, t, _1);

        for (int i = 0; i < N_TR; i++)
            HVCRA_local[i].set_dend_current(f);
    }
}



void PoolParallel::set_training_current()
{
    if (MPI_rank == 0)
    {
        current_injection_time = WAITING_TIME + generator.random(trial_duration-2*WAITING_TIME);

        std::function<double (double)> f = std::bind(&training_current, current_injection_time, _1);

        for (int i = 0; i < N_TR; i++)
            HVCRA_local[i].set_dend_current(f);
    }
}

void PoolParallel::test_ideal_chain(int num_layers, int num_trials)
{
    // initialize coordintates and ideal chain connections
    this->initialize_coordinates();
    this->write_all_coordinates();
    this->initialize_ideal_chain_connections(num_layers);
   
    this->gather_data();

    //this->print_received_invariable_connections();

    

    // write active and super synapses to files
    std::string fileActiveGraph = outputDirectory + "RA_RA_active_connections.bin";
    std::string fileSuperGraph = outputDirectory + "RA_RA_super_connections.bin";
	    	
    this->write_supersynapses(fileSuperGraph.c_str());
    this->write_active_synapses(fileActiveGraph.c_str());

    // make all neurons mature
    for (int i = 0; i < N_RA_local; i++)
        gaba_potential_local[i] = E_GABA_MATURE;

    // unused array with average dendritic spike times. Needed for gathering mature data
    std::vector<std::vector<double>> average_dendritic_spike_time; // array with average dendritic spike time in every trial

	average_dendritic_spike_time.resize(N_RA);
    
    for (int i = 0; i < num_trials; i++)
    {
        this->mature_trial();
        this->gather_mature_data(average_dendritic_spike_time);

        for (int n = 1; n <= num_layers; n++)
            for (int j = 0; j < Nss; j++)
                this->write_RA((outputDirectory + "RA/RA" + std::to_string(n*Nss + j) + "_trial" + std::to_string(i+1) + ".bin").c_str(), n*Nss + j);
        
        for (int n = 0; n <= num_layers; n++)
            for (int j = 0; j < Nss; j++)
                this->write_I((outputDirectory + "I/I" + std::to_string(n*Nss + j) + "_trial" + std::to_string(i+1) + ".bin").c_str(), n*Nss + j);

        // write all occured spikes to files
        this->write_soma_spike_times((outputDirectory + "soma_spikes_trial" + std::to_string(i+1) + ".bin").c_str());
        this->write_dend_spike_times((outputDirectory + "dend_spikes_trial" + std::to_string(i+1) + ".bin").c_str());
        
        this->randomize_after_trial();
    }   
}

void PoolParallel::continue_growth(std::string dataDir, int starting_trial, int save_freq_short, int save_freq_long)
{
	this->read_network_state(dataDir, starting_trial); // read data from file
    
    outputDirectory = dataDir;
    
    std::cout << "Trial number : " << trial_number << std::endl;
    
    bool training = true;
    
    this->chain_growth(training, save_freq_short, save_freq_long);
    
}

void PoolParallel::test_grown_chain(int num_trials, std::string dataDir, int starting_trial, std::string outputDir)
{
    this->read_network_state(dataDir, starting_trial); // read data from file
    outputDirectory = outputDir;
    this->test_mature_chain(num_trials); // test networkc
}

void PoolParallel::test_random_chain(int num_layers, int num_trials)
{
    this->initialize_coordinates();
    this->write_all_coordinates();
    this->initialize_random_chain_connections(num_layers);

    // set all neurons to be mature
    for (int i = 0; i < N_RA_local; i++)
    {
        mature_local[i] = 1;
        gaba_potential_local[i] = E_GABA_MATURE;
    }

    this->test_mature_chain(num_trials);
}

void PoolParallel::test_mature_chain(int num_trials)
{
    std::string file_soma_spikes = outputDirectory + "soma_spikes_in_trial.bin"; // file with somatic spikes in trial
    std::string file_dend_spikes = outputDirectory + "dend_spikes_in_trial.bin"; // file with dendritic spikes in trial
    std::string file_chain_test = outputDirectory + "mature_chain_test.bin"; // trial with mature chain test info

	std::vector<std::vector<double>> average_dendritic_spike_time; // array with average dendritic spike time in every trial
	std::vector<std::vector<int>> num_dendritic_spikes_in_trials; // number of dendritic spikes produced in all trials
	std::vector<std::vector<int>> num_somatic_spikes_in_trials; // number of somatic spikes produced in all trials

	average_dendritic_spike_time.resize(N_RA);
    num_dendritic_spikes_in_trials.resize(N_RA);
    num_somatic_spikes_in_trials.resize(N_RA);

    if (MPI_rank == 0)
    {
        for (int j = 0; j < N_RA; j++)
        {
            num_dendritic_spikes_in_trials[j].resize(num_trials);
            num_somatic_spikes_in_trials[j].resize(num_trials);
        }
    }
    // neurons for gabaMaturation 300117	
    /*std::vector<int> RAtoWrite = {71, 186, 187, 84, 44, 175, 219, 238, 224, 70, 288, 117, 99, 276, 23, 24, 165, 
                                  128, 184, 155, 114, 203, 257, 65, 273, 183, 294, 19, 35, 97, 142, 233, 6, 192, 
                                  248, 295, 38, 69, 207, 268, 49, 263, 132, 101, 33, 206, 90, 252, 77, 43, 293, 36, 
                                  5, 180, 282, 65, 34, 267, 208, 66, 146, 179};
    */
    /*
    // neurons for gabaMaturation 310117	
    std::vector<int> RAtoWrite = {201, 209, 124, 275, 40, 87, 66, 282, 222, 285, 115, 58, 183, 123, 244, 96, 226,
                                  110, 15, 20, 178, 215, 192, 128, 280, 38, 7, 235, 273, 258, 227, 132, 169, 172, 
                                  243, 100, 188};
    */
    // neurons for gabaMaturation 010217	
    /*std::vector<int> RAtoWrite = {179, 129, 128, 268, 130, 142, 15, 115, 273, 19, 23, 282, 29, 261, 290, 163, 292, 
                                  37, 167, 168, 169, 199, 172, 51, 182, 60, 68, 69, 256, 201, 207, 208, 209, 82, 85, 
                                  87, 90, 92, 122, 144, 226, 227, 131, 101, 81, 259, 231, 110, 114, 243, 117, 120, 250, 123, 124, 213};
    */
	/*
    // neurons for gabaMaturation270317 huxley
    std::vector<int> RAtoWrite = {179, 66, 11, 123, 173, 129, 148, 287, 199, 174, 285, 298, 144, 20, 161, 165, 205, 89, 17}; 
    */
    // neurons for gabaMaturation270317 huxley
    //std::vector<int> RAtoWrite = {179, 66, 11, 123, 173, 129, 148, 287, 199, 174, 285, 298, 144, 20, 161, 165, 205, 89, 17}; 
    //std::vector<int> ItoWrite;

	/*
    // neurons for gabaMaturation010417 hodgkin
    std::vector<int> RAtoWrite = {281, 156, 84, 52, 16, 92, 238, 75, 47, 10, 283, 171, 115, 194, 225, 78, 268, 221, 289, 104,
                                  185, 285, 287, 21, 58, 55, 229, 222, 145, 239, 123, 173, 295, 179, 240, 134, 280, 42, 228, 178, 
                                  208, 244, 294, 130, 45, 4, 217, 143, 87, 226, 148, 233, 190, 223, 255, 138, 29, 192, 290, 12, 
                                  142, 129, 150, 48, 69, 271, 174, 17, 167, 168, 273, 68, 35, 95, 163, 207, 128, 172, 231, 258, 
                                  99, 30, 100}; 
    */
    /*
    // neurons for gabaMaturation280317 huxley
    std::vector<int> RAtoWrite = {111, 253, 62, 265, 260, 8, 291, 160, 143, 64, 271, 128, 134, 84, 38, 72, 267, 34, 137, 77, 
                                  20, 188, 200, 136, 173, 13, 206, 5, 118};
    */
    /*
    // neurons for gabaMaturation040417 huxley
    std::vector<int> RAtoWrite = {85, 197, 201, 44, 262, 247, 228, 249, 185, 46, 199, 212, 64, 140, 174, 210, 236, 77, 129, 
								  15, 39, 298, 168, 216, 142, 295, 204, 13, 23, 34, 280, 186, 299, 121, 54, 269, 292, 105, 9, 
								  35, 57, 251, 100, 69, 260, 182, 136, 237, 134, 26, 66, 157, 286, 135, 193, 45, 219, 80, 20, 
								  126, 196, 211, 6, 190, 257, 81, 104, 36, 253, 25, 90, 115, 30, 183, 63, 109, 266, 202, 94, 113, 
								  222, 187, 246, 86, 206, 232, 160, 125, 240, 117, 282, 152, 19, 259, 198, 128};
    */
    // neurons for gabaMaturation130417 huxley
    std::vector<int> RAtoWrite = {51, 48, 146, 172, 132, 277, 203, 175, 275, 28, 31, 37, 140, 235, 67, 245, 21, 50, 138, 93, 76,
									228, 46, 225, 187, 231, 156, 210, 246, 148, 7, 49, 195, 74, 124, 255, 169, 152, 269, 206, 260, 
									94, 83, 259, 57, 171, 114, 23, 222, 248, 113, 165, 20, 104, 116, 59, 257, 25, 26, 89, 252, 151, 
									229, 253, 106, 176, 115, 183, 283, 30, 112, 226, 267, 139, 238, 158, 167, 95, 84, 268, 162, 111, 164, 163};
    
    std::vector<int> ItoWrite;

    for (int i = 0; i < N_I; i++)
        ItoWrite.push_back(i);

    for (int i = 0; i < num_trials; i++)
	{
        if (MPI_rank == 0)
            std::cout << "Trial " << i << std::endl;

		this->mature_trial();
		this->gather_mature_data(average_dendritic_spike_time);

        if (MPI_rank == 0)
        {
            for (int j = 0; j < N_RA; j++)
            {
                num_dendritic_spikes_in_trials[j][i] = static_cast<int>(spikes_in_trial_dend_global[j].size());
                num_somatic_spikes_in_trials[j][i] = static_cast<int>(spikes_in_trial_soma_global[j].size());
            }

        }

        for (size_t j = 0; j < RAtoWrite.size(); j++)
            this->write_RA((outputDirectory + "RA/RA" + std::to_string(RAtoWrite[j]) + "_trial" + std::to_string(i+1) + ".bin").c_str(), RAtoWrite[j]);
        
        for (size_t j = 0; j < ItoWrite.size(); j++)
            this->write_I((outputDirectory + "I/I" + std::to_string(ItoWrite[j]) + "_trial" + std::to_string(i+1) + ".bin").c_str(), ItoWrite[j]);
        
        if (i == 0)
        {
		    this->write_soma_spike_times(file_soma_spikes.c_str());
		    this->write_dend_spike_times(file_dend_spikes.c_str());
        }
		this->randomize_after_trial();
	}
	
	// process dendritic spikes

	std::vector<double> mean_burst_time; // average of dendritic spike time
	std::vector<double> std_burst_time; // standard deviation of dendritic spike time
    std::vector<double> average_num_dend_spikes_in_trial; // average number of dendritic spikes in trial
    std::vector<double> average_num_soma_spikes_in_trial; // average number of somatic spikes in trials

    std::vector<int> num_trials_with_dend_spikes; // number of trials in which neuron produced dendritic spikes

	mean_burst_time.resize(N_RA);
	std_burst_time.resize(N_RA);
    average_num_dend_spikes_in_trial.resize(N_RA);
    average_num_soma_spikes_in_trial.resize(N_RA);
    num_trials_with_dend_spikes.resize(N_RA);

	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_RA; i++)
		{
            average_num_dend_spikes_in_trial[i] = std::accumulate(num_dendritic_spikes_in_trials[i].begin(), num_dendritic_spikes_in_trials[i].end(), 0.0) 
                                                / static_cast<double>(num_trials);

            average_num_soma_spikes_in_trial[i] = std::accumulate(num_somatic_spikes_in_trials[i].begin(), num_somatic_spikes_in_trials[i].end(), 0.0) 
                                                / static_cast<double>(num_trials);
            
            num_trials_with_dend_spikes[i] = static_cast<int>(average_dendritic_spike_time[i].size());
			
            /*
            if (i == 288)
            {
                std::cout << "Number of somatic spikes in all trials: " << std::endl;
                for (int j = 0; j < num_trials; j++)
                    std::cout << num_somatic_spikes_in_trials[i][j] << '\t';

                std::cout << '\n' << std::endl;
                
                std::cout << "Number of dendritic spikes in all trials: " << std::endl;
                for (int j = 0; j < num_trials; j++)
                    std::cout << num_dendritic_spikes_in_trials[i][j] << "\t";

                std::cout << std::endl << std::endl;
                
                std::cout << "Number of trials in which neuron produced dendritic bursts: " << num_trials_with_dend_spikes[i] << "\n" << std::endl;
                std::cout << "All dendritic burst times in these trials: " << std::endl;
                
                for (size_t j = 0; j < average_dendritic_spike_time[i].size(); j++)
                    std::cout << average_dendritic_spike_time[i][j] << "\t";

                std::cout << std::endl;

            }
            */

            if (num_trials_with_dend_spikes[i] > 0)
            {

                //for (int j = 0; j < (int) average_dendritic_spike_time[i].size(); j++)
                  //  printf("average_dendritic_spike_time[%d][%d] = %f\n", i, j, average_dendritic_spike_time[i][j]);
				mean_burst_time[i] = std::accumulate(average_dendritic_spike_time[i].begin(), average_dendritic_spike_time[i].end(), 0.0) / 
                                        static_cast<double>(num_trials_with_dend_spikes[i]);

            }

			else
				mean_burst_time[i] = -1;
			// calculate standard deviation of burst times

			double accum = 0;
			double mean = mean_burst_time[i];

			std::for_each(average_dendritic_spike_time[i].begin(), average_dendritic_spike_time[i].end(), [&accum, mean](const double t)
			{
				accum += (t - mean) * (t - mean);
			});

			if (static_cast<int>(average_dendritic_spike_time[i].size() > 1))
				std_burst_time[i] = sqrt(accum / (static_cast<double>(average_dendritic_spike_time[i].size()) - 1));
			else
				std_burst_time[i] = -1;


		}
	}

	this->write_chain_test(num_trials, num_trials_with_dend_spikes, average_num_dend_spikes_in_trial, average_num_soma_spikes_in_trial, 
                           mean_burst_time, std_burst_time, file_chain_test.c_str());
}

void PoolParallel::mature_trial()
{
	int some_RA_neuron_fired_soma_local;

	int some_I_neuron_fired_local;
	int some_RA_neuron_fired_soma_global;

	int some_I_neuron_fired_global;

    trial_number++;
    
    internal_time = 0;
    network_time = internal_time + network_update_frequency;

	// set training current
	double current_injection_time = 100; //ms
    this->set_training_current(current_injection_time);

	// initialize update arrays and fired indicators
	std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);

    some_RA_neuron_fired_soma_local = 0;
	some_RA_neuron_fired_soma_global = 0;

    some_I_neuron_fired_local = 0;
    some_I_neuron_fired_global = 0;
	
    // evolve dynamics
    for (int t = 1; t < size; t++)
	{
		internal_time += timeStep;
		
		for (int i = 0; i < N_RA_local; i++)
		{
            // set GABA potential
            HVCRA_local[i].set_Ei(gaba_potential_local[i]);
            
            // Debraband step
            HVCRA_local[i].Debraband_step_no_target_update();
            
            // if some neuron produced somatic spike
            if (HVCRA_local[i].get_fired_soma())
            {
                spikes_in_trial_soma_local[i].push_back(internal_time);
                
                // if neuron is mature, update conductances of targets
                if (mature_local[i] == 1)
                {
                    some_RA_neuron_fired_soma_local = 1;
                    // loop over all inhibitory targets of fired neurons
                    size_t num_I_targets = syn_ID_RA_I_local[i].size();
                    for (size_t j = 0; j < num_I_targets; j++)
                    {
                        int syn_ID = syn_ID_RA_I_local[i][j];
                        update_Ge_I_local[syn_ID] += weights_RA_I_local[i][j];

                    }
                    
                    // loop over all excitatory targets
                    size_t num_RA_targets = active_synapses_local[i].size();
                    
                    //std::cout << "Neuron fired: " << Id_RA_local[i] << " num_RA_targets: " << num_RA_targets << std::endl;

                    for (size_t j = 0; j < num_RA_targets; j++)
                    {
                        int syn_ID = active_synapses_local[i][j];
                        update_Ge_RA_local[syn_ID] += weights_local[i][syn_ID];
                        //std::cout << "Neuron fired: " << Id_RA_local[i] << " target: " << syn_ID << " synaptic weight: " << weights_local[i][syn_ID] << std::endl;
                    }
                }
            } 

            if (HVCRA_local[i].get_fired_dend())
            {
                spikes_in_trial_dend_local[i].push_back(internal_time);
            }
		}

		for (int i = 0; i < N_I_local; i++)
		{
            HVCI_local[i].DP8_step_no_target_update();
            
            //  if some I neuron spikes, change conductance update array
            if (HVCI_local[i].get_fired())
            {
                //printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
                some_I_neuron_fired_local = 1;
                spikes_in_trial_interneuron_local[i].push_back(internal_time);

                size_t num_RA_targets = syn_ID_I_RA_local[i].size();
                // loop over all targets of fired neurons
                for (size_t j = 0; j < num_RA_targets; j++)
                {
                    int syn_ID = syn_ID_I_RA_local[i][j];
                    update_Gi_RA_local[syn_ID] += weights_I_RA_local[i][j];
                    //printf("Rank = %d; i = %d; update_Gi_RA_local[%d] = %f; weights_I_RA_local[%d][%d] = %f\n", MPI_rank, i, syn_ID,
                     //   update_Gi_RA_local[syn_ID], weights_I_RA_local[fired_ID][j], fired_ID, j);
                }
            }
		}

        // if we need to update network state
        // get if any neurons fired in some process

        if (internal_time > network_time)
        {
            MPI_Allreduce(&some_RA_neuron_fired_soma_local, &some_RA_neuron_fired_soma_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_I_neuron_fired_local, &some_I_neuron_fired_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            if (some_I_neuron_fired_global > 0)
            {
            // sum update array and send to all processes

                MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                for (int i = 0; i < N_RA_local; i++)
                {
                    HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
                }

				// reset fired indicators and arrays
				some_I_neuron_fired_global = 0;
				some_I_neuron_fired_local = 0;


				std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            }

            //if (some_RA_neuron_fired_global == 1)
            //    printf("Rank %d; some_RA_neuron_fired_global: %d\n", MPI_rank, some_RA_neuron_fired_global);

            // if somatic compartment of any neuron in the pool fired, update synaptic conductances
             if (some_RA_neuron_fired_soma_global > 0)
             {
            // sum all update arrays and send to all processes

                MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                // now update excitatory conductances of all neurons
                for (int i = 0; i < N_RA_local; i++)
                {
                    HVCRA_local[i].raiseE(update_Ge_RA_global[Id_RA_local[i]]); // update conductance
				}

                for (int i = 0; i < N_I_local; i++)
                {
                    HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);
                }

				// reset conductance arrays and fired indicators
				some_RA_neuron_fired_soma_global = 0;
				some_RA_neuron_fired_soma_local = 0;

				std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
				std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
            }

            network_time += network_update_frequency;
        } // end network update


        //MPI_Barrier(MPI_COMM_WORLD);
    } // end evolve dynamics
}

void PoolParallel::run_trials_no_save(int num_trials)
{
    bool training = true;

	while (trial_number < num_trials)
    {
        this->trial(training);

        trial_number++;
        this->randomize_after_trial();
    }
}


void PoolParallel::run_trials_with_save(int num_trials)
{
	std::string fileNumSynapses = outputDirectory + "num_synapses.bin";
	
    std::string fileActiveGraph = outputDirectory + "RA_RA_connections.bin";
    std::string fileSuperGraph = outputDirectory + "RA_RA_super_connections.bin";
    std::string fileTimeSoma = outputDirectory + "spike_times_soma.bin";
    std::string fileTimeDend = outputDirectory + "spike_times_dend.bin";
    std::string fileTimeInterneuron = outputDirectory + "spike_times_interneuron.bin";
    std::string fileWeightsGraph = outputDirectory + "weights.bin";
    std::string fileMaturationGraph = outputDirectory + "mature.bin";

    bool training = true;

	while (trial_number < num_trials)
    {
        this->trial(training);

        this->gather_data();

        this->write_num_synapses(fileNumSynapses.c_str());
        this->write_soma_spike_times(fileTimeSoma.c_str());
        this->write_dend_spike_times(fileTimeDend.c_str());
        this->write_interneuron_spike_times(fileTimeInterneuron.c_str());
       
        this->write_supersynapses(fileSuperGraph.c_str());
        this->write_maturation_info(fileMaturationGraph.c_str());
    
        this->write_weights(fileWeightsGraph.c_str());
        this->write_active_synapses(fileActiveGraph.c_str());
        this->write_maturation_info((outputDirectory + "mature" + std::to_string(trial_number) + ".bin").c_str());
        
        trial_number++;
        this->randomize_after_trial();
    }

}

void PoolParallel::chain_growth(bool training, int save_freq_short, int save_freq_long)
{
    std::string fileRA = outputDirectory + "RA.bin";
    std::string fileI = outputDirectory + "I.bin";
    
	std::string fileNumSynapses = outputDirectory + "num_synapses.bin";
	std::string fileWeightStatistics = outputDirectory + "weight_statistics.bin";
	
    std::string fileActiveGraph = outputDirectory + "RA_RA_active_connections.bin";
    std::string fileSuperGraph = outputDirectory + "RA_RA_super_connections.bin";
    std::string fileTimeSoma = outputDirectory + "spike_times_soma.bin";
    std::string fileTimeDend = outputDirectory + "spike_times_dend.bin";
    std::string fileTimeInterneuron = outputDirectory + "spike_times_interneuron.bin";
    std::string fileMaturationGraph = outputDirectory + "mature.bin";
    
    std::string fileMaturationTimeSequence = outputDirectory + "maturation_time_sequence.bin";

    std::string fileWeightsGraph = outputDirectory + "weights.bin";
    std::string fileWeightsSourceToTarget = outputDirectory + "weightsTimeSequence.bin";
    std::string RAdir = outputDirectory + "RAneurons/";
    std::string Idir = outputDirectory + "Ineurons/";

    bool data_gathered; // indicator if data was already gathered

    int trial_to_add_neurons = 2500;
    int N = 100;

    std::vector<int> RAtoWrite{0, 1, 2, 3};

    std::vector<int> source{0, 1, 2, 3};
    std::vector<int> target{};

    for (int i = 0; i < N_RA; i++)
        target.push_back(i);

    // make all neurons able to make output connections
    //this->set_all_mature();

    while (true)
    {
       /* if (trial_number == trial_to_add_neurons)
        {
            if (MPI_rank == 0)
                std::cout << "Adding new neurons!" << std::endl;
            
            this->add_new_neurons(N);
        }
        */
		//break;
		data_gathered = false;

        //this->trial(training);
		this->trial_burst_stdp(training);

        if (MPI_rank == 0)
            std::cout << "Trial " << trial_number + 1 << std::endl;
        //pool.gather_data();

        //pool.statistics();


        if (trial_number % save_freq_short == 0)
        {
			this->gather_data();
			data_gathered = true;

           	this->write_num_synapses(fileNumSynapses.c_str());
		    this->write_soma_spike_times(fileTimeSoma.c_str());
            this->write_dend_spike_times(fileTimeDend.c_str());
            //this->write_dend_spike_times((outputDirectory + "spike_times_dend_" + std::to_string(trial_number) + ".bin").c_str());
            this->write_interneuron_spike_times(fileTimeInterneuron.c_str());
           
	    	this->write_supersynapses(fileSuperGraph.c_str());
	    	this->write_active_synapses(fileActiveGraph.c_str());
			this->write_maturation_info(fileMaturationGraph.c_str());
		
            this->write_weights_time_sequence_from_source_to_target(source, target, fileWeightsSourceToTarget.c_str());
            this->write_maturation_time_sequence(target, fileMaturationTimeSequence.c_str());

            for (int i = 0; i < (int) RAtoWrite.size(); i++)
	    	{
				this->write_RA((RAdir + "RA" + std::to_string(RAtoWrite[i]) + ".bin").c_str(), RAtoWrite[i]);
	    	}

			//fileMature = fileMaturePerm + "mature" + std::to_string(count) + ".bin";
			//pool.write_mature(fileMature.c_str());
			//pool.write_pajek_all(filePajekAll.c_str());
        }
	
			

		if (trial_number % save_freq_long == 0)
		{
			if (!data_gathered)
				this->gather_data();

            this->write_weight_statistics(fileWeightStatistics.c_str());
            this->write_weights(fileWeightsGraph.c_str());
            this->write_active_synapses((outputDirectory + "RA_RA_active_connections_" + std::to_string(trial_number) + "_.bin").c_str());
            this->write_supersynapses((outputDirectory + "RA_RA_super_connections_" + std::to_string(trial_number) + "_.bin").c_str());
			this->write_weights((outputDirectory + "weights_" + std::to_string(trial_number) + "_.bin").c_str());
			this->write_maturation_info((outputDirectory + "mature_" + std::to_string(trial_number) + "_.bin").c_str());
		}
		
		//break;    
		//for (int i = 0; i < N_I; i++)
	    //{
	      //  fileAllIneurons = Idir + "I" + std::to_string(i) + ".bin";
		//	pool.write_I(fileAllIneurons.c_str(), i);
	    //}

        trial_number++;
        this->randomize_after_trial();
		
		//break;
    }

}

void PoolParallel::trial_burst_stdp(int training)
{
	int some_RA_neuron_fired_soma_local;
	int some_RA_neuron_fired_dend_local;

	int some_I_neuron_fired_local;
	int some_RA_neuron_fired_soma_global;
	int some_RA_neuron_fired_dend_global;

	int some_I_neuron_fired_global;

    std::vector<double> spike_times_fired_dend_local;
    
    std::vector<int> RA_neurons_fired_dend_global;
    std::vector<int> RA_neurons_fired_dend_realID;

    internal_time = 0;
    network_time = network_update_frequency;

    // advance internal time of each neuron
    for (int i = 0; i < N_RA_local; i++)
        num_trials_after_replacement_local[i] += 1;

    //printf("network time = %f\n", network_time);

    // if training trial
    if (training > 0)
        this->set_training_current();

	// initialize update arrays and fired indicators
	std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);

    some_RA_neuron_fired_soma_local = 0;
	some_RA_neuron_fired_soma_global = 0;

	some_RA_neuron_fired_dend_local = 0;
    some_RA_neuron_fired_dend_global = 0;

    some_I_neuron_fired_local = 0;
    some_I_neuron_fired_global = 0;
	
    // evolve dynamics
    for (int t = 1; t < size; t++)
	{
		internal_time += timeStep;
		
		for (int i = 0; i < N_RA_local; i++)
		{
            // set GABA potential
            HVCRA_local[i].set_Ei(gaba_potential_local[i]);
            
            // Debraband step
            HVCRA_local[i].Debraband_step_no_target_update();
            
            // if some neuron produced somatic spike, do LTD for all previous dendritic spikes
            if (HVCRA_local[i].get_fired_soma())
            {
                spikes_in_trial_soma_local[i].push_back(internal_time);
                
                // allow neuron output only if neuron is mature
                if (mature_local[i] == 1)
                {

                    some_RA_neuron_fired_soma_local = 1;
                    // loop over all inhibitory targets of fired neurons
                    size_t num_I_targets = syn_ID_RA_I_local[i].size();
                    for (size_t j = 0; j < num_I_targets; j++)
                    {
                        int syn_ID = syn_ID_RA_I_local[i][j];
                        update_Ge_I_local[syn_ID] += weights_RA_I_local[i][j];

                    }
                    
                    // loop over all excitatory targets
                    size_t num_RA_targets = active_synapses_local[i].size();
                    for (size_t j = 0; j < num_RA_targets; j++)
                    {
                        int syn_ID = active_synapses_local[i][j];
                        update_Ge_RA_local[syn_ID] += weights_local[i][syn_ID];
                    }

                }
                //for (int j = 0; j < last_soma_spikes_local[i].size(); j++)
                //{
                  //  printf("My rank = %d; Soma spikes of neuron %d:  spike_time = %f\n", MPI_rank, Id_RA_local[i],
                    //    last_soma_spikes_local[i][j]);

                //}

                //printf("My rank = %d; Soma neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], spike_times_soma_local[i]);

            } // end if get fired soma
            
            // if some neuron produced dendritic spike, store this neuron in array
            if (HVCRA_local[i].get_fired_dend())
            {
                some_RA_neuron_fired_dend_local = 1;
                spikes_in_trial_dend_local[i].push_back(internal_time);
                //printf("My rank = %d; RA neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], internal_time);
                RA_neurons_fired_dend_realID.push_back(Id_RA_local[i]);
                spike_times_fired_dend_local.push_back(internal_time);
                //printf("My rank = %d; Dend neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], spike_times_dend_local[i]);

                // if neuron is saturated apply LTD only to supersynapses
                if (static_cast<int>(supersynapses_local[i].size()) == Nss)
                {
                    for (size_t k = 0; k < supersynapses_local[i].size(); k++)
                    {
                        int supersynapse_id = supersynapses_local[i][k];
                        
                        if ( (spikes_in_trial_dend_global[supersynapse_id].size() > 0) && (Id_RA_local[i] != supersynapse_id) )
                        {
                            double dt = internal_time - spikes_in_trial_dend_global[supersynapse_id].back();

                            //std::cout << "dt in saturated LTD = " << dt << std::endl;

                            if (dt < STDP_WINDOW)
                            {
                                //double w = weights_local[i][supersynapse_id];     
                       
                                LTD_burst(weights_local[i][supersynapse_id], dt);
                                
                                //printf("LTD from saturated neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n",
                                //            Id_RA_local[i], supersynapse_id, spikes_in_trial_soma_local[i].back(), spike_times_dend_global[supersynapse_id],
                                //            dt, weights_local[i][supersynapse_id] - w);
                                
                                update_synapse(i, supersynapse_id);
                                
                            }
                        }
                    }

					// if some supersynapse desaturated, update all synapses
                	if (static_cast<int>(supersynapses_local[i].size()) < Nss)
						for (int j = 0; j < N_RA; j++)
							this->update_synapse(i, j);
                }
                // if not saturated apply LTD rule with the last dendritic spike of all neurons and add glutamate to all neuron except the fired one
                else
                {
                    for (int j = 0; j < N_RA; j++)
                    {
                        if ( (spikes_in_trial_dend_global[j].size() > 0) && (Id_RA_local[i] != j) )
                        {
                            double dt = internal_time - spikes_in_trial_dend_global[j].back();
                            
                            //std::cout << "dt in LTD = " << dt << std::endl;

                            if (dt < STDP_WINDOW)
                            {

                                //double w = weights_local[i][j];

                                LTD_burst(weights_local[i][j], dt);
                     
                                //printf("LTD from neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n", Id_RA_local[i], j,
                                //	spikes_in_trial_soma_local[i].back(), spike_times_dend_global[j], dt, weights_local[i][j] - w);
                                
                                update_synapse(i, j);
                                    
                            }
                        }
                    }
                }

            } // end if get_fired_dend

		} // end for i -> N_RA_local

		for (int i = 0; i < N_I_local; i++)
		{
            HVCI_local[i].DP8_step_no_target_update();
            
            //  if some I neuron spikes, change conductance update array
            if (HVCI_local[i].get_fired())
            {
                //printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
                some_I_neuron_fired_local = 1;
                spikes_in_trial_interneuron_local[i].push_back(internal_time);

                size_t num_RA_targets = syn_ID_I_RA_local[i].size();
                // loop over all targets of fired neurons
                for (size_t j = 0; j < num_RA_targets; j++)
                {
                    int syn_ID = syn_ID_I_RA_local[i][j];
                    update_Gi_RA_local[syn_ID] += weights_I_RA_local[i][j];
                    //printf("Rank = %d; i = %d; update_Gi_RA_local[%d] = %f; weights_I_RA_local[%d][%d] = %f\n", MPI_rank, i, syn_ID,
                     //   update_Gi_RA_local[syn_ID], weights_I_RA_local[fired_ID][j], fired_ID, j);
                }
            }
		} // end if i -> N_I_local

        // if we need to update network state
        // get if any neurons fired in some process

        if (internal_time > network_time)
        {
            //if (MPI_rank == 0)
            //{
            //    std::cout << "internal time = " << internal_time << " > network_time = " << network_time << std::endl; 
            //}

            MPI_Allreduce(&some_RA_neuron_fired_soma_local, &some_RA_neuron_fired_soma_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_RA_neuron_fired_dend_local, &some_RA_neuron_fired_dend_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_I_neuron_fired_local, &some_I_neuron_fired_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            if (some_I_neuron_fired_global > 0)
            {
            // sum update array and send to all processes

                MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                for (int i = 0; i < N_RA_local; i++)
                {
                    HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
                }
            	
				// update conductance arrays and fired indicators
				some_I_neuron_fired_local = 0;
            	some_I_neuron_fired_global = 0;
				
				std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            }


            //if (some_RA_neuron_fired_global == 1)
            //    printf("Rank %d; some_RA_neuron_fired_global: %d\n", MPI_rank, some_RA_neuron_fired_global);

            // if somatic compartment of any neuron in the pool fired, update synaptic conductances
             if (some_RA_neuron_fired_soma_global > 0)
             {
            // sum all update arrays and send to all processes

                MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                // now update excitatory conductances of all neurons
                for (int i = 0; i < N_RA_local; i++)
                {
                    HVCRA_local[i].raiseE(update_Ge_RA_global[Id_RA_local[i]]); // update conductance
				}

                for (int i = 0; i < N_I_local; i++)
                {
                    HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);
                }
				
				// update conductance arrays and fired indicators
            	some_RA_neuron_fired_soma_local = 0;
	        	some_RA_neuron_fired_soma_global = 0;
				
				std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
				
				std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
            }

            // if dendritic compartment of any neuron in the pool fired, update weights
            if (some_RA_neuron_fired_dend_global > 0)
            {
                // gather all bursted neurons
                this->gather_bursts(RA_neurons_fired_dend_realID, RA_neurons_fired_dend_global, spike_times_fired_dend_local);

              	/* 
                if (MPI_rank == 0)
                {
					for (size_t i = 0; i < RA_neurons_fired_dend_global.size(); i++)
                        printf("neuron %d bursted; spike_time_dend = %f\n", RA_neurons_fired_dend_global[i], 
								spikes_in_trial_dend_global[RA_neurons_fired_dend_global[i]].back());
                }
                */

                // apply LTP rule for all recent dendritic bursts and dendritic spike of fired neurons
                // but only if source neuron is mature

                for (int i = 0; i < N_RA_local; i++)
                {
                    // if neuron is saturated apply LTP only if dendritic spike occured in supersynapse
                    if (( static_cast<int>(supersynapses_local[i].size()) == Nss) && (mature_local[i] == 1) )
                    {
                        for (size_t j = 0; j < RA_neurons_fired_dend_global.size(); j++)
                        {
                            int fired_ID = RA_neurons_fired_dend_global[j];

                            std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
                                        supersynapses_local[i].end(), fired_ID);

                            if (pos!=supersynapses_local[i].end())
                            {
                                for  (size_t k = 0; k < spikes_in_trial_dend_local[i].size(); k++)
                                {
                                    double dt = spikes_in_trial_dend_global[fired_ID].back() - spikes_in_trial_dend_local[i][k];

									//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
                                    if (dt < STDP_WINDOW)
                                    {
                                        //double w = weights_local[i][fired_ID];
                                        
                                        LTP_burst(weights_local[i][fired_ID], dt);

                                        //printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
                                         //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
                                         //           dt, weights_local[i][fired_ID] - w);
                                        
                                        update_synapse(i, fired_ID);
									}
                                }
                            }
                        }

                    }
                    // if not saturated and is mature apply LTP for all dendritic spikes
                    else if (mature_local[i] == 1)
                    {
                        for (size_t j = 0; j < RA_neurons_fired_dend_global.size(); j++)
                        {
                            int fired_ID = RA_neurons_fired_dend_global[j];
                            // don't allow self-to-self connections
                            if (fired_ID != Id_RA_local[i])
                            {
                                // loop over last somatic spikes
                                for  (size_t k = 0; k < spikes_in_trial_dend_local[i].size(); k++)
                                {
                                    double dt = spikes_in_trial_dend_global[fired_ID].back() - spikes_in_trial_dend_local[i][k];

									//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
                                    if (dt < STDP_WINDOW)
                                    {
                                        //double w = weights_local[i][fired_ID];
                                        
                                        LTP_burst(weights_local[i][fired_ID], dt);
                                        
                                        //printf("LTP from neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
                                          //          Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
                                            //        dt, weights_local[i][fired_ID] - w);
                                        
                                        update_synapse(i, fired_ID);
                                  		
									}
                                }
                            }
                        }
                   }

                } // end for i -> N_RA_local

               // check if we need axon remodeling

                for (int i = 0; i < N_RA_local; i++)
                {
                    if ( (static_cast<int>(supersynapses_local[i].size()) == Nss) && (remodeled_local[i] == 0) )
                    {
					    this->axon_remodeling(i);

				    }
                }
	        	
				// update fired arrays and indicators
				some_RA_neuron_fired_dend_local = 0;
            	some_RA_neuron_fired_dend_global = 0;
            	
				RA_neurons_fired_dend_global.clear();
            	RA_neurons_fired_dend_realID.clear();

                spike_times_fired_dend_local.clear();
            } // end if some_neuron_dend_fired
            
            network_time += network_update_frequency;
            
            //printf("network_time = %f\n", network_time);

        }


        //MPI_Barrier(MPI_COMM_WORLD);
    }
    this->potentiation_decay();
    //printf("After potentiation decay")
    this->update_all_synapses();
	
	// calculate new gaba reverse potential
	// update rate
	for (int i = 0; i < N_RA_local; i++)
	{
		num_bursts_in_recent_trials[i].push_front(static_cast<int>(spikes_in_trial_dend_local[i].size()));
        
        firing_rate_local[i] = std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].begin() + RATE_WINDOW_SHORT, 0.0)
                                    / static_cast<double>(RATE_WINDOW_SHORT);


		//if (std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].end(), 0.0) >= 1)
		//{
		//	std::cout << "Recent bursts of neuron " << Id_RA_local[i] << " :" << std::endl;
			
		//	for (size_t j = 0; j < num_bursts_in_recent_trials[i].size(); j++)
		//		std::cout << num_bursts_in_recent_trials[i][j] << "\t";
			
		//	std::cout << std::endl;
		//}
	}
    // update GABA potential based on firing rates
	this->update_Ei();
    
    // gather all neurons that are to be replaced
    this->gather_neurons_2replace();
    
    // if some neurons are to be replaced, replace them
    if (replace_real_id_global.size() > 0)
        this->replace_neurons();


}



void PoolParallel::trial(int training)
{
	int some_RA_neuron_fired_soma_local;
	int some_RA_neuron_fired_dend_local;

	int some_I_neuron_fired_local;
	int some_RA_neuron_fired_soma_global;
	int some_RA_neuron_fired_dend_global;

	int some_I_neuron_fired_global;

    std::vector<double> spike_times_fired_dend_global;
    std::vector<double> spike_times_fired_dend_local;
    
    std::vector<int> RA_neurons_fired_dend_global;
    std::vector<int> RA_neurons_fired_dend_realID;

    internal_time = 0;
    network_time = network_update_frequency;

    // advance internal time of each neuron
    for (int i = 0; i < N_RA_local; i++)
        num_trials_after_replacement_local[i] += 1;

    //printf("network time = %f\n", network_time);

    // if training trial
    if (training > 0)
        this->set_training_current();

	// initialize update arrays and fired indicators
	std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);

    some_RA_neuron_fired_soma_local = 0;
	some_RA_neuron_fired_soma_global = 0;

	some_RA_neuron_fired_dend_local = 0;
    some_RA_neuron_fired_dend_global = 0;

    some_I_neuron_fired_local = 0;
    some_I_neuron_fired_global = 0;
	
    // evolve dynamics
    for (int t = 1; t < size; t++)
	{
		internal_time += timeStep;
		
		for (int i = 0; i < N_RA_local; i++)
		{
            // set GABA potential
            HVCRA_local[i].set_Ei(gaba_potential_local[i]);
            
            // Debraband step
            HVCRA_local[i].Debraband_step_no_target_update();
            
            // if some neuron produced somatic spike, do LTD for all previous dendritic spikes
            if (HVCRA_local[i].get_fired_soma())
            {
                some_RA_neuron_fired_soma_local = 1;
                spikes_in_trial_soma_local[i].push_back(internal_time);

                // loop over all inhibitory targets of fired neurons
                size_t num_I_targets = syn_ID_RA_I_local[i].size();
                for (size_t j = 0; j < num_I_targets; j++)
                {
                    int syn_ID = syn_ID_RA_I_local[i][j];
                    update_Ge_I_local[syn_ID] += weights_RA_I_local[i][j];

                }
				
                // loop over all excitatory targets
                size_t num_RA_targets = active_synapses_local[i].size();
                for (size_t j = 0; j < num_RA_targets; j++)
                {
                    int syn_ID = active_synapses_local[i][j];
                    update_Ge_RA_local[syn_ID] += weights_local[i][syn_ID];
                }

                // if neuron is saturated apply LTD only to supersynapses
                if (static_cast<int>(supersynapses_local[i].size()) == Nss)
                {
                    for (size_t k = 0; k < supersynapses_local[i].size(); k++)
                    {
                        int supersynapse_id = supersynapses_local[i][k];

                        double dt = internal_time - spike_times_dend_global[supersynapse_id];

						//std::cout << "dt in saturated LTD = " << dt << std::endl;

                        if ( (Id_RA_local[i] != supersynapse_id) && (dt < STDP_WINDOW) )
                        {
                            //double w = weights_local[i][supersynapse_id];     
                   
                            LTD(weights_local[i][supersynapse_id], dt);
                            
                            //printf("LTD from saturated neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n",
                            //            Id_RA_local[i], supersynapse_id, spikes_in_trial_soma_local[i].back(), spike_times_dend_global[supersynapse_id],
                            //            dt, weights_local[i][supersynapse_id] - w);
                            
                            update_synapse(i, supersynapse_id);
                       		
					   	}
                    }

					// if some supersynapse desaturated, update all synapses
                	if (static_cast<int>(supersynapses_local[i].size()) < Nss)
						for (int j = 0; j < N_RA; j++)
							this->update_synapse(i, j);
                }
                // if not saturated apply LTD rule with the last dendritic spike of all neurons and add glutamate to all neuron except the fired one
                else
                {
                    for (int j = 0; j < N_RA; j++)
                    {
                        double dt = internal_time - spike_times_dend_global[j];
						
						//std::cout << "dt in LTD = " << dt << std::endl;

                        if ( (Id_RA_local[i] != j) && (dt < STDP_WINDOW) )
                        {

                            //double w = weights_local[i][j];

                           	LTD(weights_local[i][j], dt);
                 
                            //printf("LTD from neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n", Id_RA_local[i], j,
			    		    //	spikes_in_trial_soma_local[i].back(), spike_times_dend_global[j], dt, weights_local[i][j] - w);
                           	
                            update_synapse(i, j);
								
                        }
                    }
                }


                //for (int j = 0; j < last_soma_spikes_local[i].size(); j++)
                //{
                  //  printf("My rank = %d; Soma spikes of neuron %d:  spike_time = %f\n", MPI_rank, Id_RA_local[i],
                    //    last_soma_spikes_local[i][j]);

                //}

                //printf("My rank = %d; Soma neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], spike_times_soma_local[i]);

            } // end if get fired soma
            
            // if some neuron produced dendritic spike, store this neuron in array
            if (HVCRA_local[i].get_fired_dend())
            {
                some_RA_neuron_fired_dend_local = 1;
                spikes_in_trial_dend_local[i].push_back(internal_time);
                //printf("My rank = %d; RA neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], internal_time);
                RA_neurons_fired_dend_realID.push_back(Id_RA_local[i]);
                spike_times_fired_dend_local.push_back(internal_time);
                //printf("My rank = %d; Dend neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], spike_times_dend_local[i]);

            }

		} // end for i -> N_RA_local

		for (int i = 0; i < N_I_local; i++)
		{
            HVCI_local[i].DP8_step_no_target_update();
            
            //  if some I neuron spikes, change conductance update array
            if (HVCI_local[i].get_fired())
            {
                //printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
                some_I_neuron_fired_local = 1;
                spikes_in_trial_interneuron_local[i].push_back(internal_time);

                size_t num_RA_targets = syn_ID_I_RA_local[i].size();
                // loop over all targets of fired neurons
                for (size_t j = 0; j < num_RA_targets; j++)
                {
                    int syn_ID = syn_ID_I_RA_local[i][j];
                    update_Gi_RA_local[syn_ID] += weights_I_RA_local[i][j];
                    //printf("Rank = %d; i = %d; update_Gi_RA_local[%d] = %f; weights_I_RA_local[%d][%d] = %f\n", MPI_rank, i, syn_ID,
                     //   update_Gi_RA_local[syn_ID], weights_I_RA_local[fired_ID][j], fired_ID, j);
                }
            }
		} // end if i -> N_I_local

        // if we need to update network state
        // get if any neurons fired in some process

        if (internal_time > network_time)
        {
            //if (MPI_rank == 0)
            //{
            //    std::cout << "internal time = " << internal_time << " > network_time = " << network_time << std::endl; 
            //}

            MPI_Allreduce(&some_RA_neuron_fired_soma_local, &some_RA_neuron_fired_soma_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_RA_neuron_fired_dend_local, &some_RA_neuron_fired_dend_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_I_neuron_fired_local, &some_I_neuron_fired_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            if (some_I_neuron_fired_global > 0)
            {
            // sum update array and send to all processes

                MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                for (int i = 0; i < N_RA_local; i++)
                {
                    HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
                }
            	
				// update conductance arrays and fired indicators
				some_I_neuron_fired_local = 0;
            	some_I_neuron_fired_global = 0;
				
				std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            }


            //if (some_RA_neuron_fired_global == 1)
            //    printf("Rank %d; some_RA_neuron_fired_global: %d\n", MPI_rank, some_RA_neuron_fired_global);

            // if somatic compartment of any neuron in the pool fired, update synaptic conductances
             if (some_RA_neuron_fired_soma_global > 0)
             {
            // sum all update arrays and send to all processes

                MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                // now update excitatory conductances of all neurons
                for (int i = 0; i < N_RA_local; i++)
                {
                    HVCRA_local[i].raiseE(update_Ge_RA_global[Id_RA_local[i]]); // update conductance
				}

                for (int i = 0; i < N_I_local; i++)
                {
                    HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);
                }
				
				// update conductance arrays and fired indicators
            	some_RA_neuron_fired_soma_local = 0;
	        	some_RA_neuron_fired_soma_global = 0;
				
				std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
				
				std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
            }

            // if dendritic compartment of any neuron in the pool fired, update weights
            if (some_RA_neuron_fired_dend_global > 0)
            {
                // send all fired neurons
                int num_RA_fired_dend_local = RA_neurons_fired_dend_realID.size();
                int num_RA_fired_dend_global;
                //printf("Rank %d; num_RA_fired_local: %d\n", MPI_rank, num_RA_fired_local);

                int* recvcounts = new int[MPI_size];
                int* displs = new int[MPI_size];

                // get total number of fired neurons
                MPI_Allreduce(&num_RA_fired_dend_local, &num_RA_fired_dend_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                RA_neurons_fired_dend_global.resize(num_RA_fired_dend_global);
                spike_times_fired_dend_global.resize(num_RA_fired_dend_global);

                //printf("Rank %d; fired RA num: %d\n", MPI_rank, num_RA_fired_local);

                //if (MPI_rank == 0)
                //{
                //   printf("Master; fired RA num global: %d\n", num_RA_fired_global);
                //   printf("Master; RA_neurons_fired_global.size(): %d\n", RA_neurons_fired_global.size());
                //}

                // get array with number of fired neurons in each process
                MPI_Allgather(&num_RA_fired_dend_local, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, MPI_COMM_WORLD);

                displs[0] = 0;
                for (int i = 1; i < MPI_size; i++)
                {
                    displs[i] = displs[i-1] + recvcounts[i-1];
                }


                //for (int i = 0; i < RA_neurons_fired_realID.size(); i++)
                //        printf("Rank %d; fired RA neuron: %d\n", MPI_rank, RA_neurons_fired_realID[i]);

                // get fired neurons
                MPI_Allgatherv(&RA_neurons_fired_dend_realID[0], num_RA_fired_dend_local, MPI_INT,
                    &RA_neurons_fired_dend_global[0], recvcounts, displs, MPI_INT, MPI_COMM_WORLD);
                
                MPI_Allgatherv(&spike_times_fired_dend_local[0], num_RA_fired_dend_local, MPI_DOUBLE,
                    &spike_times_fired_dend_global[0], recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
                /*
                if (MPI_rank == 0)
                 {
                    printf("Master; RA_neurons_fired_dend_global.size(): %d\n", RA_neurons_fired_dend_global.size());
                    for (int i = 0; i < RA_neurons_fired_dend_global.size(); i++)
                        printf("Rank %d; Dend fired RA neuron: %d; spike_time = %f\n", MPI_rank, RA_neurons_fired_dend_global[i],
                                internal_time);
                }
                 */
                // change spike times
                for (size_t i = 0; i < RA_neurons_fired_dend_global.size(); i++)
                {
                    spike_times_dend_global[RA_neurons_fired_dend_global[i]] = internal_time;
                    spike_times_dend_global[RA_neurons_fired_dend_global[i]] = spike_times_fired_dend_global[i];
                    
                }

                delete [] recvcounts;
                delete [] displs;

                /*
                if (MPI_rank == 0)
                {   for (int i = 0; i < N_RA; i++)
                        printf("neuron %d; spike_time_dend = %f\n", i ,spike_times_dend_global[i]);
                }
                */

                // apply LTP rule for last NUM_SOMA_SPIKES and dendritic spike of fired neurons
                for (int i = 0; i < N_RA_local; i++)
                {
                    // if neuron is saturated apply LTP only if dendritic spike occured in supersynapse
                    if (static_cast<int>(supersynapses_local[i].size()) == Nss)
                    {
                        for (size_t j = 0; j < RA_neurons_fired_dend_global.size(); j++)
                        {
                            int fired_ID = RA_neurons_fired_dend_global[j];

                            std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
                                        supersynapses_local[i].end(), fired_ID);

                            if (pos!=supersynapses_local[i].end())
                            {
                                for  (size_t k = 0; k < spikes_in_trial_soma_local[i].size(); k++)
                                {
                                    double dt = internal_time - spikes_in_trial_soma_local[i][k];

                                    if (dt < STDP_WINDOW)
                                    {
                                        //double w = weights_local[i][fired_ID];
                                        
                                        LTP(weights_local[i][fired_ID], dt);

                                        //printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
                                         //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
                                         //           dt, weights_local[i][fired_ID] - w);
                                        
                                        update_synapse(i, fired_ID);
									}
                                }
                            }
                        }

                    }
                    // if not saturated apply LTP for all dendritic spikes
                    else
                    {
                        for (size_t j = 0; j < RA_neurons_fired_dend_global.size(); j++)
                        {
                            int fired_ID = RA_neurons_fired_dend_global[j];
                            // don't allow self-to-self connections
                            if (fired_ID != Id_RA_local[i])
                            {
                                // loop over last somatic spikes
                                for  (size_t k = 0; k < spikes_in_trial_soma_local[i].size(); k++)
                                {
                                    double dt = internal_time - spikes_in_trial_soma_local[i][k];

                                    if (dt < STDP_WINDOW)
                                    {
                                        //double w = weights_local[i][fired_ID];
                                        
                                        LTP(weights_local[i][fired_ID], dt);
                                        
                                        //printf("LTP from neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
                                          //          Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
                                            //        dt, weights_local[i][fired_ID] - w);
                                        
                                        update_synapse(i, fired_ID);
                                  		
									}
                                }
                            }
                        }
                   }

                } // end for i -> N_RA_local

               // check if we need axon remodeling

                for (int i = 0; i < N_RA_local; i++)
                {
                    if ( (static_cast<int>(supersynapses_local[i].size()) == Nss) && (remodeled_local[i] == 0) )
                    {
					    this->axon_remodeling(i);

				    }
                }
	        	
				// update fired arrays and indicators
				some_RA_neuron_fired_dend_local = 0;
            	some_RA_neuron_fired_dend_global = 0;
            	
				RA_neurons_fired_dend_global.clear();
            	RA_neurons_fired_dend_realID.clear();

                spike_times_fired_dend_local.clear();
                spike_times_fired_dend_global.clear();
            } // end if some_neuron_dend_fired
            
            network_time += network_update_frequency;
            
            //printf("network_time = %f\n", network_time);

        }


        //MPI_Barrier(MPI_COMM_WORLD);
    }
    this->potentiation_decay();
    //printf("After potentiation decay")
    this->update_all_synapses();
	
	// calculate new gaba reverse potential
	// update rate
	for (int i = 0; i < N_RA_local; i++)
	{
		num_bursts_in_recent_trials[i].push_front(static_cast<int>(spikes_in_trial_dend_local[i].size()));
        
        firing_rate_local[i] = std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].begin() + RATE_WINDOW_SHORT, 0.0)
                                    / static_cast<double>(RATE_WINDOW_SHORT);


		//if (std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].end(), 0.0) >= 1)
		//{
		//	std::cout << "Recent bursts of neuron " << Id_RA_local[i] << " :" << std::endl;
			
		//	for (size_t j = 0; j < num_bursts_in_recent_trials[i].size(); j++)
		//		std::cout << num_bursts_in_recent_trials[i][j] << "\t";
			
		//	std::cout << std::endl;
		//}
	}
    // update GABA potential based on firing rates
	this->update_Ei();
    
    // gather all neurons that are to be replaced
    this->gather_neurons_2replace();
    
    // if some neurons are to be replaced, replace them
    if (replace_real_id_global.size() > 0)
        this->replace_neurons();


}

void PoolParallel::gather_bursts(std::vector<int>& RA_neurons_fired_dend_realID, std::vector<int>& RA_neurons_fired_dend_global, 
                                 std::vector<double>& spike_times_fired_dend_local)
{
        std::vector<double> spike_times_fired_dend_global;
        int num_RA_fired_dend_local = RA_neurons_fired_dend_realID.size();
        int num_RA_fired_dend_global;
        //printf("Rank %d; num_RA_fired_local: %d\n", MPI_rank, num_RA_fired_local);

        int* recvcounts = new int[MPI_size];
        int* displs = new int[MPI_size];

        // get total number of fired neurons
        MPI_Allreduce(&num_RA_fired_dend_local, &num_RA_fired_dend_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        RA_neurons_fired_dend_global.resize(num_RA_fired_dend_global);
        spike_times_fired_dend_global.resize(num_RA_fired_dend_global);

        //printf("Rank %d; fired RA num: %d\n", MPI_rank, num_RA_fired_local);

        //if (MPI_rank == 0)
        //{
        //   printf("Master; fired RA num global: %d\n", num_RA_fired_global);
        //   printf("Master; RA_neurons_fired_global.size(): %d\n", RA_neurons_fired_global.size());
        //}

        // get array with number of fired neurons in each process
        MPI_Allgather(&num_RA_fired_dend_local, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, MPI_COMM_WORLD);

        displs[0] = 0;
        for (int i = 1; i < MPI_size; i++)
        {
            displs[i] = displs[i-1] + recvcounts[i-1];
        }


        //for (int i = 0; i < RA_neurons_fired_realID.size(); i++)
        //        printf("Rank %d; fired RA neuron: %d\n", MPI_rank, RA_neurons_fired_realID[i]);

        // get fired neurons
        MPI_Allgatherv(&RA_neurons_fired_dend_realID[0], num_RA_fired_dend_local, MPI_INT,
            &RA_neurons_fired_dend_global[0], recvcounts, displs, MPI_INT, MPI_COMM_WORLD);
        
        MPI_Allgatherv(&spike_times_fired_dend_local[0], num_RA_fired_dend_local, MPI_DOUBLE,
            &spike_times_fired_dend_global[0], recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        /*
        if (MPI_rank == 0)
         {
            printf("Master; RA_neurons_fired_dend_global.size(): %d\n", RA_neurons_fired_dend_global.size());
            for (int i = 0; i < RA_neurons_fired_dend_global.size(); i++)
                printf("Rank %d; Dend fired RA neuron: %d; spike_time = %f\n", MPI_rank, RA_neurons_fired_dend_global[i],
                        internal_time);
        }
         */
        // change spike times
        for (size_t i = 0; i < RA_neurons_fired_dend_global.size(); i++)
        {
			//if (MPI_rank == 0)
			//	std::cout << "neuron " << RA_neurons_fired_dend_global[i] << " bursted at " << spike_times_fired_dend_global[i] << std::endl;
            
			spikes_in_trial_dend_global[RA_neurons_fired_dend_global[i]].push_back(spike_times_fired_dend_global[i]);   
        }

        delete [] recvcounts;
        delete [] displs;
}

void PoolParallel::gather_neurons_2replace()
{
    // gather all neurons that are to be replaced
    // get number of neurons to replace
    int num_2replace_local = static_cast<int>(replace_local_id_local.size()); // number of neurons to replace on the process

    int num_2replace_total; // total number of neurons to replace

    MPI_Allreduce(&num_2replace_local, &num_2replace_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    replace_local_id_global.resize(num_2replace_total);
    replace_real_id_global.resize(num_2replace_total);
    replace_process_rank_global.resize(num_2replace_total);

    // gather ids and process ranks of neurons to replace
    int* recvcounts = new int[MPI_size];
    int* displs = new int[MPI_size];

    // get array with number of neurons to replace in each process
    MPI_Allgather(&num_2replace_local, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, MPI_COMM_WORLD);

    displs[0] = 0;
    for (int i = 1; i < MPI_size; i++)
    {
        displs[i] = displs[i-1] + recvcounts[i-1];
    }

    // get neurons to replace
    MPI_Allgatherv(&replace_local_id_local[0], num_2replace_local, MPI_INT,
        &replace_local_id_global[0], recvcounts, displs, MPI_INT, MPI_COMM_WORLD);
    
    MPI_Allgatherv(&replace_real_id_local[0], num_2replace_local, MPI_INT,
        &replace_real_id_global[0], recvcounts, displs, MPI_INT, MPI_COMM_WORLD);
    
    MPI_Allgatherv(&replace_process_rank_local[0], num_2replace_local, MPI_INT,
        &replace_process_rank_global[0], recvcounts, displs, MPI_INT, MPI_COMM_WORLD);

    delete[] recvcounts;
    delete[] displs;
}

void PoolParallel::kill_neuron(int local_id, int global_id, int process_rank)
{
    // if process contains neuron to replace
    if (MPI_rank == process_rank)
    {
        // zero all local weight matrices and delete all existing synapses from and out of the neuron
        // zero weights from neuron
        for (int i = 0; i < N_RA; i++)
            std::fill(weights_local[local_id].begin(), weights_local[local_id].end(), 0.0);

        // delete active and super synapses
        supersynapses_local[local_id].clear();
        active_synapses_local[local_id].clear();
        
        // set active and supersynapse indicators to zero
        active_indicators_local[local_id].reset();
        supersynapses_indicators_local[local_id].reset();

        // reset gaba_potential, remodeled and mature states
        gaba_potential_local[local_id] = E_GABA_IMMATURE;
        remodeled_local[local_id] = 0;
        mature_local[local_id] = 0;

        gaba_reached_mature_local[local_id] = 0;
	
        firing_rate_local[local_id] = 0.0;
    
        num_trials_after_replacement_local[local_id] = 0;

        // clear all recent bursts
        for (int j = 0; j < RATE_WINDOW_LONG; j++)
            num_bursts_in_recent_trials[local_id].push_back(0);
    }

    // do for all processes
    // erase all super and active synapses to a neuron
    for (int i = 0; i < N_RA_local; i++)
    {
        std::vector<int>::iterator pos = std::find(active_synapses_local[i].begin(), active_synapses_local[i].end(), global_id);

        if (pos != active_synapses_local[i].end())
        {
            active_synapses_local[i].erase(pos);
            
            // set active indicator to false
            active_indicators_local[i][global_id] = false;
        }

        pos = std::find(supersynapses_local[i].begin(), supersynapses_local[i].end(), global_id);

        if (pos != supersynapses_local[i].end())
        {
            // if neuron was saturated, set remodeled to false
            if (static_cast<int>(supersynapses_local[i].size()) == Nss)
            {
                if (remodeled_local[i] != 1)
                    std::cerr << "In replace_neuron. Neuron " << Id_RA_local[i] << " has maximum supersynapses but is nor remodeled!" << std::endl;
                
            }
            else
            {
                if (remodeled_local[i] == 1)
                    std::cerr << "In replace_neuron. Neuron " << Id_RA_local[i] << " is remodeled but its number of output supersynapses is "
                              <<  supersynapses_local[i].size() << " which is smaller than maximum " << Nss << std::endl;
                
            }

            supersynapses_local[i].erase(pos);
            remodeled_local[i] = 0;   
            
            // set supersynapses indicator to false
            supersynapses_indicators_local[i][global_id] = false;
        }

    }

    // zero weights to neuron
    for (int i = 0; i < N_RA_local; i++)
        weights_local[i][global_id] = 0.0;
}

void PoolParallel::initialize_coordinates_for_replaced_neuron(int global_id)
{
    if (MPI_rank == 0)
    {
        double xx; // temporary x-coordinate
        double yy; // temporary y-coordinate

        bool close; // are neurons too close or not

        do
        {
            close = false;
            xx = 0.5*SIDE / (sqrt(N_I)+1) + generator.random(SIDE - SIDE / (sqrt(N_I)+1));
            yy = 0.5*SIDE / (sqrt(N_I)+1) + generator.random(SIDE - SIDE / (sqrt(N_I)+1));

            // check distances to all I neurons

            for (int j = 0; j < N_I; j++)
            {
                if (distance(xx,yy,xx_I[j],yy_I[j]) < MIN_INTERNEURON_DISTANCE)
                {
                    close = true;
                    break;
                }
            }

            // check distances to all previous RA neurons
            if (!close)
            {
                for (int j = 0; j < N_RA; j++)
                {
                    if (distance(xx, yy, xx_RA[j], yy_RA[j]) < MIN_INTERNEURON_DISTANCE)
                    {
                        close = true;
                        break;
                    }

                }
            }
        } while(close);

        xx_RA[global_id] = xx;
        yy_RA[global_id] = yy;
        
    }
}

void PoolParallel::initialize_connections_for_replaced_neuron(int global_id)
{
    if (MPI_rank == 0)
    {
        // erase all previous connections from replaced neuron to interneurons
        weights_RA_I_global[global_id].clear();
        syn_ID_RA_I_global[global_id].clear();

        // erase all previous connections from interneurons to replaced neuron
        for (int i = 0; i < N_I; i++)
        {
            std::vector<int>::iterator pos = std::find(syn_ID_I_RA_global[i].begin(), syn_ID_I_RA_global[i].end(), global_id);

            // if interneuron targets replaced neuron, erase connection
            if (pos != syn_ID_I_RA_global[i].end())
            {
                int index_in_array = std::distance(syn_ID_I_RA_global[i].begin(), pos);

                syn_ID_I_RA_global[i].erase(pos);
                weights_I_RA_global[i].erase(weights_I_RA_global[i].begin() + index_in_array);
            }

        }

      // make new connections to interneurons
        for (int j = 0; j < N_I; j++)
        {
             if (generator.random(1) < p_RA2I(global_id,j))
             {
                 double G = this->sample_Ge2i();

                 weights_RA_I_global[global_id].push_back(G);
                 syn_ID_RA_I_global[global_id].push_back(j);
             }
         }


		// make new connections from interneurons to replaced neuron
		for (int i = 0; i < N_I; i++)
     	{
             if (generator.random(1) < p_I2RA(i,global_id))
             {
                 double G = this->sample_Gi2e();

                 weights_I_RA_global[i].push_back(G);
                 syn_ID_I_RA_global[i].push_back(global_id);
             }
		 }
    }
}

void PoolParallel::replace_neurons()
{
    // write active, supersynapses and weights before replacing neurons
    this->gather_data();

    this->write_weights((outputDirectory + "weights_before_replacement_trial" + std::to_string(trial_number) + ".bin").c_str());
    this->write_active_synapses((outputDirectory + "RA_RA_active_connections_before_replacement_trial" + std::to_string(trial_number) + ".bin").c_str());
    this->write_supersynapses((outputDirectory + "RA_RA_super_connections_before_replacement_trial" + std::to_string(trial_number) + ".bin").c_str());
    this->write_maturation_info((outputDirectory + "mature_before_replacement_trial" + std::to_string(trial_number) + ".bin").c_str());

    // loop through all neurons that are to be replaced
    for (size_t i = 0; i < replace_local_id_global.size(); i++)
    {
        // zero all weight matrices and delete all existing synapses from and out of the neurons
        this->kill_neuron(replace_local_id_global[i], replace_real_id_global[i], replace_process_rank_global[i]);

        // initialize new coordinates and create new connections
        this->initialize_coordinates_for_replaced_neuron(replace_real_id_global[i]);
        this->initialize_connections_for_replaced_neuron(replace_real_id_global[i]);  
    }

    this->send_connections();

    // write all changes to files
    if (MPI_rank == 0)
    {
        // write new coordinates
        std::string fileRAxy = outputDirectory + "RA_xy_after_replacement_trial" + std::to_string(trial_number) + ".bin";

        this->write_coordinates_RA(fileRAxy.c_str());
        
        // write invariable connections
        std::string fileRA2I = outputDirectory + "RA_I_connections_after_replacement_trial" + std::to_string(trial_number) + ".bin";
        std::string fileI2RA = outputDirectory + "I_RA_connections_after_replacement_trial" + std::to_string(trial_number) + ".bin";
        std::string filePajekFixed = outputDirectory + "fixed_after_replacement_trial" + std::to_string(trial_number) + ".net";

        this->write_invariable_synapses(fileRA2I.c_str(), fileI2RA.c_str());
        this->write_pajek_fixed(filePajekFixed.c_str());

        // write replaced neurons to a file
        std::string fileReplaced = outputDirectory + "replaced_neurons_trial" + std::to_string(trial_number) + ".bin";
        
        this->write_replaced_neurons(replace_real_id_global, fileReplaced.c_str());
    }

    // clear all arrays with neurons to be replaced
    replace_local_id_local.clear();
    replace_real_id_local.clear();
    replace_process_rank_local.clear();
    replace_local_id_global.clear();
    replace_real_id_global.clear();
    replace_process_rank_global.clear();
    
    // write active, supersynapses and weights after replacing neurons
    this->gather_data();
    
    this->write_weights((outputDirectory + "weights_after_replacement_trial" + std::to_string(trial_number) + ".bin").c_str());
    this->write_active_synapses((outputDirectory + "RA_RA_active_connections_after_replacement_trial" + std::to_string(trial_number) + ".bin").c_str());
    this->write_supersynapses((outputDirectory + "RA_RA_super_connections_after_replacement_trial" + std::to_string(trial_number) + ".bin").c_str());
    this->write_maturation_info((outputDirectory + "mature_after_replacement_trial" + std::to_string(trial_number) + ".bin").c_str());
}

void PoolParallel::add_new_neurons(int N)
{
    this->gather_data();

    this->write_weights((outputDirectory + "weights_before_neuron_addition.bin").c_str());
    this->write_active_synapses((outputDirectory + "RA_RA_active_connections_before_neuron_addition.bin").c_str());
    this->write_supersynapses((outputDirectory + "RA_RA_super_connections_before_neuron_addition.bin").c_str());
    this->write_maturation_info((outputDirectory + "mature_before_neuron_addition.bin").c_str());
	
    int new_remain = N; // new neurons remaining to be distributed
    // distribute new neurons among processes

    std::vector<int> num_new_on_each_process; // number of new neurons added to each process

    num_new_on_each_process.resize(MPI_size); // array with number of RA neurons per process

    // first make number of HVC(RA) neurons equal on each process if it differs
    // find process where number of neurons is smaller
    int first_process_with_smaller_neuron_num = -1; // id of the process starting from which all processes have number of neurons less by one

    for (int i = 1; i < MPI_size; i++)
    {
        if (N_RA_sizes[i] < N_RA_sizes[i-1])
            first_process_with_smaller_neuron_num = i;       
    }
    
    // if such process is found, add 1 neuron to all processes with id >= first_process_with_smaller_neuron_num
    if (first_process_with_smaller_neuron_num > 0)
    {
        for (int i = first_process_with_smaller_neuron_num; i < MPI_size; i++)
        {
            N_RA_sizes[i]++;
            num_new_on_each_process[i]++; 
        }
        new_remain -= MPI_size - first_process_with_smaller_neuron_num; // fix number of neurons to be distributed
    }
    
    // now number of neurons is already equal. distribute all remaining new neurons
	for (int i = 0; i < MPI_size; i++)
	{
		num_new_on_each_process[i] += new_remain / MPI_size;
        N_RA_sizes[i] += new_remain / MPI_size;
	}
	new_remain = new_remain % MPI_size;
	
    int j = 0;

	// distribute RA neurons
	while (new_remain > 0)
	{
		N_RA_sizes[j] += 1;
        num_new_on_each_process[j] += 1;
		new_remain -= 1;
		j += 1;

		if (j >= MPI_size)
			j -= MPI_size;
	}
   
    int N_RA_local_old = N_RA_local; // number of HVC(RA) neurons on each process before new neurons were added
    int N_RA_old = N_RA; // total number of HVC(RA) neurons before new neurons were added
    N_RA_local = N_RA_sizes[MPI_rank]; // new number of neurons on each process after new neurons were added
    N_RA = N_RA + N; // new number of HVC(RA) neurons in the network after new neurons were added

	printf("My rank = %d; neurons added = %d; N_RA_local_old = %d; N_RA_local = %d; N_RA_old = %d; N_RA = %d;\n", MPI_rank, num_new_on_each_process[MPI_rank], N_RA_local_old, N_RA_local, N_RA_old, N_RA);

    // resize all arrays
    this->resize_arrays_for_RA(N_RA_local, N_RA);

    
    // assign id to added neurons
    Id_RA_local.resize(N_RA_local);

    for (int i = 0; i < num_new_on_each_process[MPI_rank]; i++)
    {
        int start_id_on_each_process = N_RA_old; // real neuron id which is assigned to the first new neuron on each process   

        for (int j = 0; j < MPI_rank; j++)
            start_id_on_each_process += num_new_on_each_process[j];

        Id_RA_local[N_RA_local_old + i] = start_id_on_each_process + i;
    }

    // resize global array containing all ids of HVC(RA) neurons
    Id_RA_global.resize(N_RA);

	int* recvcounts_id_RA = new int[MPI_size];
	int* displs_id_RA = new int[MPI_size];

	if (MPI_rank == 0)
	{
		// send connections to all processes

		recvcounts_id_RA[0] = N_RA_sizes[0];
		displs_id_RA[0] = 0;

		for (int i = 1; i < MPI_size; i++)
		{
			recvcounts_id_RA[i] = N_RA_sizes[i];
			displs_id_RA[i] = displs_id_RA[i-1] + recvcounts_id_RA[i-1];
		}
    }

	// get global array with neuronal ids
	MPI_Gatherv(&Id_RA_local[0], N_RA_local, MPI_INT, &Id_RA_global[0], recvcounts_id_RA, displs_id_RA, MPI_INT, 0, MPI_COMM_WORLD);

	delete[] recvcounts_id_RA;
	delete[] displs_id_RA;

    this->write_global_index_array((outputDirectory + "global_index_array_after_neuron_addition.bin").c_str());
    /*
    std::cout << "My rank = " << MPI_rank << "\nId_RA_local = ";
    for (size_t i = 0; i < Id_RA_local.size(); i++)
        std::cout << Id_RA_local[i] << "\t";
    std::cout << std::endl;
    */
    if (MPI_rank == 0)
    {
        std::cout << "\nId_RA_global = ";
        for (size_t i = 0; i < Id_RA_global.size(); i++)
            std::cout << Id_RA_global[i] << "\t";
        std::cout << std::endl;
	}
    
    std::fill(num_bursts_in_recent_trials.begin() + N_RA_local_old, num_bursts_in_recent_trials.end(), intBuffer(RATE_WINDOW_LONG));
	
    for (size_t i = N_RA_local_old; i < num_bursts_in_recent_trials.size(); i++)
		for (int j = 0; j < RATE_WINDOW_LONG; j++)
			num_bursts_in_recent_trials[i].push_back(0);

    // initialize coordinates for new neurons
    this->initialize_coordinates_for_added_neurons(N_RA_old);

    // initialize_connections for new neurons
    this->initialize_connections_for_added_neurons(N_RA_old);

    // initialize new neurons

    for (int i = N_RA_local_old; i < N_RA_local; i++)
	{	
		HVCRA_local[i].set_noise_generator(&generator);
		HVCRA_local[i].set_white_noise(white_noise_mean_soma, white_noise_std_soma, white_noise_mean_dend, white_noise_std_dend);
		HVCRA_local[i].set_dynamics(trial_duration, timeStep);
	
        gaba_potential_local[i] = E_GABA_IMMATURE;      
    }

    this->gather_data();

    this->write_weights((outputDirectory + "weights_after_neuron_addition.bin").c_str());
    this->write_active_synapses((outputDirectory + "RA_RA_active_connections_after_neuron_addition.bin").c_str());
    this->write_supersynapses((outputDirectory + "RA_RA_super_connections_after_neuron_addition.bin").c_str());
    this->write_maturation_info((outputDirectory + "mature_after_neuron_addition.bin").c_str());
}

void PoolParallel::resize_arrays_for_I(int n_local, int n_total)
{
    HVCI_local.resize(n_local);

    Id_I_local.resize(n_local);

    // connections and their ids
	weights_I_RA_local.resize(n_local);
    weights_I_RA_global.resize(n_total);
	
	syn_ID_I_RA_local.resize(n_local);
    syn_ID_I_RA_global.resize(n_total);

    // update arrays for conductances
    update_Ge_I_local.resize(n_total);
    update_Ge_I_global.resize(n_total);
	
    // spikes in trials
    spikes_in_trial_interneuron_global.resize(n_total);
    spikes_in_trial_interneuron_local.resize(n_local);
}

void PoolParallel::resize_arrays_for_RA(int n_local, int n_total)
{
	HVCRA_local.resize(n_local);

    Id_RA_local.resize(n_local);
    
    // connections and their ids
    weights_global.resize(n_total);
    weights_local.resize(n_local);
    
    weights_RA_I_local.resize(n_local);
    weights_RA_I_global.resize(n_total);
	
    syn_ID_RA_I_local.resize(n_local);
    syn_ID_RA_I_global.resize(n_total);

    // active and super synapses
	supersynapses_global.resize(n_total);
	active_synapses_global.resize(n_total);

    supersynapses_local.resize(n_local);
	active_synapses_local.resize(n_local);
    
    active_indicators_local.resize(n_local);
    supersynapses_indicators_local.resize(n_local);

    // update arrays for conductances
	update_Ge_RA_local.resize(n_total);
    update_Gi_RA_local.resize(n_total);

    update_Ge_RA_global.resize(n_total);
    update_Gi_RA_global.resize(n_total);
	
    // spikes in trials
    spike_times_dend_global.resize(n_total);
    
    spikes_in_trial_soma_global.resize(n_total);
    spikes_in_trial_dend_global.resize(n_total);

    spikes_in_trial_soma_local.resize(n_local);
    spikes_in_trial_dend_local.resize(n_local);
  
    // axon remodeling
	remodeled_local.resize(n_local);
	remodeled_global.resize(n_total);

    // mature state
	mature_local.resize(n_local);
	mature_global.resize(n_total);

	gaba_potential_local.resize(n_local);
	gaba_potential_global.resize(n_total);

    gaba_reached_mature_local.resize(n_local);
	
    firing_rate_local.resize(n_local);
	firing_rate_global.resize(n_total);

    num_trials_after_replacement_local.resize(n_local);

    // initialize circular buffers for somatic spikes
	num_bursts_in_recent_trials.resize(n_local);

    // initialize last dendritic spike times
    std::fill(spike_times_dend_global.begin(), spike_times_dend_global.end(), -200.0);

    // initialize arrays
	for (int i = 0; i < n_local; i++)
	{
        weights_local[i].resize(n_total);
        active_indicators_local[i].resize(n_total);
        supersynapses_indicators_local[i].resize(n_total);
	}

	for (int i = 0; i < n_total; i++)
		weights_global[i].resize(n_total);
}

void PoolParallel::update_Ei()
{
	for (int i = 0; i < N_RA_local; i++)
	{
        // calculate firing rate in large window:
        double firing_rate_long = std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].end(), 0.0)
                                                                                            / static_cast<double>(RATE_WINDOW_LONG);

	// if not a training neuron, GABA didnt hit the bottom and firing rate exceeds gaba rate threshold, decrease gaba
		if ( (Id_RA_local[i] >= N_TR) && (!gaba_reached_mature_local[i]) && (firing_rate_local[i] >= GABA_RATE_THRESHOLD))
		{
			gaba_potential_local[i] = gaba_potential_local[i] - GABA_DOWN * static_cast<int>(spikes_in_trial_dend_local[i].size()); // decrease gaba potential
					
			if (gaba_potential_local[i] <= E_GABA_MATURE) // if gaba potential hit the bottom, make neuron mature
			{
				gaba_potential_local[i] = E_GABA_MATURE;
				gaba_reached_mature_local[i] = true;
			}

		}
		// if not a training neuron and if firing rate in wide window exceeds maturation threshold, make mature
		if ( (Id_RA_local[i] >= N_TR) && (firing_rate_long >= MATURATION_RATE_THRESHOLD) && (mature_local[i] == 0) )
            mature_local[i] = 1;

		// if not a training neuron and if firing rate in wide window is smaller than death threshold, replace neuron with a new one
		else if ( (Id_RA_local[i] >= N_TR) && (firing_rate_long <= DEATH_RATE_THRESHOLD) && (num_trials_after_replacement_local[i] >= RATE_WINDOW_LONG))
		{
            replace_local_id_local.push_back(i);
            replace_real_id_local.push_back(Id_RA_local[i]);
            replace_process_rank_local.push_back(MPI_rank);
		}
	}
}

void PoolParallel::potentiation_decay()
{
    for (int i = 0; i < N_RA_local; i++)
    {
        for (int j = 0; j < N_RA; j++)
        {
            if (weights_local[i][j] >= SUPERSYNAPSE_THRESHOLD)
                weights_local[i][j] *= BETA_SUPERSYNAPSE;
            else
                weights_local[i][j] *= BETA;
        }
    }
}

void PoolParallel::axon_remodeling(int i)
{
	// erase all active synapses if they are not among supersynapses

	std::vector<int> active_to_erase; // active synapses to erase

	// find all active synapses needed to be erased
    for (size_t j = 0; j < active_synapses_local[i].size(); j++)
    {
        int syn_ID = active_synapses_local[i][j];

        // check if active synapse is among supersynapses
        std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
                                supersynapses_local[i].end(), syn_ID);
        // if not erase it
        if (pos == supersynapses_local[i].end())
        {
            active_to_erase.push_back(active_synapses_local[i][j]);
            active_indicators_local[i][syn_ID] = 0;
        }
    }

	// erase synapses
	for (size_t j = 0; j < active_to_erase.size(); j++)
	{
		std::vector<int>::iterator pos = std::find(active_synapses_local[i].begin(), active_synapses_local[i].end(), active_to_erase[j]);

		if (pos != active_synapses_local[i].end())
			active_synapses_local[i].erase(pos);
		else
			std::cout << "Active synapse " << Id_RA_local[i] << " -> " << active_to_erase[j]
					  << " which should be removed in axon remodeling is not found!" << std::endl;
	}

    remodeled_local[i] = 1;
}

void PoolParallel::update_all_synapses()
{
    //printf("Process %d; Updating all synapses after potentiation decay\n", MPI_rank);
    for (int i = 0; i < N_RA_local; i++)
    {
        for (int j = 0; j < N_RA; j++)
        {
            this->update_synapse(i, j);
			
        }
    }

}

void PoolParallel::update_synapse(int i, int j)
{
	double w = weights_local[i][j];

    if ( (w >= ACTIVATION_THRESHOLD) && (active_indicators_local[i][j] == 0) && (remodeled_local[i] == 0) )
    {
       // printf("Activated synapse from %d onto %d; weight = %f\n", Id_RA_local[i], j, w);
        active_indicators_local[i][j] = 1;
        active_synapses_local[i].push_back(j);
    }

    if ( (w >= SUPERSYNAPSE_THRESHOLD) && (supersynapses_indicators_local[i][j] == 0) && (static_cast<int>(supersynapses_local[i].size()) < Nss) )
    {
       // printf("Activated supersynapse from %d onto %d; weight = %f\n", Id_RA_local[i], j, w);
        supersynapses_indicators_local[i][j] = 1;
        supersynapses_local[i].push_back(j);
    }

    if ( (w < SUPERSYNAPSE_THRESHOLD) && (supersynapses_indicators_local[i][j] == 1) )
    {
       // printf("Deactivated supersynapse from %d onto %d; weight = %f\n", Id_RA_local[i], j, w);
        supersynapses_indicators_local[i][j] = 0;
        remodeled_local[i] = 0;
        std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
                                                    supersynapses_local[i].end(), j);

        if (pos!= supersynapses_local[i].end())
            supersynapses_local[i].erase(pos);
        else
            std::cout << "Supersynapse " << Id_RA_local[i] << " -> " << j << " to be erased is not found!" << std::endl;
    }

    if ( (w < ACTIVATION_THRESHOLD) && (active_indicators_local[i][j] == 1) )
    {
      //  printf("Deactivated synapse from %d onto %d; w = %f\n", Id_RA_local[i], j, w);
        active_indicators_local[i][j] = 0;

        std::vector<int>::iterator pos = std::find(active_synapses_local[i].begin(),
                                                    active_synapses_local[i].end(), j);

        if (pos!= active_synapses_local[i].end())
            active_synapses_local[i].erase(pos);
        else
            std::cout << "Active synapse " << Id_RA_local[i] << " -> " << j << " to be erased is not found!" << std::endl;


    }

}

double PoolParallel::sample_Ge2i()
{
    return Gei_mean + generator.normal_distribution() * Gei_std;
}

double PoolParallel::sample_Gi2e()
{
    return Gie_mean + generator.normal_distribution() * Gie_std;
}

double PoolParallel::p_RA2I(int i_RA, int j_I)
{
	double prob;
    double d;
    d = distance(xx_RA[i_RA], yy_RA[i_RA], xx_I[j_I], yy_I[j_I]);

	if (d < MIN_INTERNEURON_DISTANCE)
		return 0;
	else
	{
		//prob = A_RA2I * exp(-d / LAMBDA_RA2I) + B_RA2I * exp(-(d - MEAN_RA2I)*(d - MEAN_RA2I)/(2*SIGMA_RA2I*SIGMA_RA2I));
		
		prob = A_RA2I * exp(-d * d / (2 * SIGMA_RA2I * SIGMA_RA2I));
		
		//printf("p = %f\n", prob);
 	}
	return prob;
}

double PoolParallel::p_I2RA(int i_I, int j_RA)
{
	double prob;
    double d;
    d = distance(xx_RA[j_RA], yy_RA[j_RA], xx_I[i_I], yy_I[i_I]);

	if (d < MIN_INTERNEURON_DISTANCE)
		return 0;
	else
	{
		prob = B_I2RA * exp(-d*d / (2 * SIGMA_I2RA * SIGMA_I2RA));
		//printf("p = %f\n", prob);
 	}
	return prob;
}

void PoolParallel::LTP_burst(double &w, double t)
{
	if (t <= T_P)
    {
		w = w + R * A_P * t / T_P;
    }
	else
    {
		w = w + R * A_P * exp(-(t - T_P) / TAU_P);
    }

    if (w < 0)
        w = 0;

    if (w > WEIGHT_MAX)
        w = WEIGHT_MAX;
}


void PoolParallel::LTD_burst(double &w, double t)
{
	if (t <= T_D)
	{
		w = w - w * R * A_D * t / T_D;
		

       // std::cout << "w = " << w << std::endl;
	}
	else
	{
		w = w - w * R * A_D * exp(-(t - T_D) / TAU_D);

	}
	if (w < 0)
    	w = 0;
}

void PoolParallel::LTP(double &w, double t)
{
    //printf("LTP\n");
	if (t <= T_P)
    {
		//std::cout << "Weight before LTP = " << w << std::endl;
        //if ((w + R * A_P * G_P * ((1 + F_0) * t / T_P - F_0))<0)
          //  printf("LTP. Linear interval. Weight became negative. w = %f\n", w);
        //double temp = R * A_P * G_P * ((1 + F_0) * t / T_P - F_0);
        //if (temp<0)
         //   printf("Negative LTP!!! t = %f; temp = %f\n", t, temp);
        //else
         //   printf("Positive LTP!!! t = %f; temp = %f\n", t, temp);
		//w = w + temp;
			w = w + R * A_P * ((1 + F_0) * t / T_P - F_0);
	
		//std::cout << "Weight after LTP = " << w << std::endl;

        //std::cout << "w = " << w << std::endl;
    }
	else
    {
        //if ((w + R * A_P * G_P * exp(-(t - T_P) / TAU_P))<0)
         //   printf("LTP. Exponential interval. Weight became negative. w = %f\n", w);
        //double temp = R * A_P * G_P * exp(-(t - T_P) / TAU_P);

    	//if (temp < 0)
         //   printf("Negative LTP!!! t = %f; temp = %f\n", t, temp);
        //else
         //   printf("Positive LTP!!! t = %f; temp = %f\n", t, temp);

        //w = w + temp;
		//std::cout << "Weight before exp LTP = " << w << std::endl;
   		
		w = w + R * A_P * exp(-(t - T_P) / TAU_P);

		//std::cout << "Weight after exp LTP = " << w << std::endl;
        //std::cout << "w = " << w << std::endl;
    }

    if (w < 0)
        w = 0;

    if (w > WEIGHT_MAX)
        w = WEIGHT_MAX;
}

void PoolParallel::LTD(double &w, double t)
{
	if (t <= T_D)
	{

        //if ((w - R * A_D * w * ((1 - F_0) * t / T_D + F_0))<0)
            //printf("LTD. Linear interval. Weight became negative. w = %f\n", w);
		//w = w - R * A_D * w * ((1 - F_0) * t / T_D + F_0);
		//w = w - R * A_D * ((1 - F_0) * t / T_D + F_0);
		//std::cout << "Weight before LTD = " << w << std::endl;
			w = w - w * R * A_D * ( F_0 + (1 -  F_0) * t / T_D);
		
		//std::cout << "Weight after LTD = " << w << std::endl;
		if (w < 0)
            w = 0;

       // std::cout << "w = " << w << std::endl;
	}
	else
	{
        //if ((w - R * A_D * w * exp(-(t - T_D) / TAU_D))<0)
          //  printf("LTD. Exponential interval. Weight became negative. w = %f\n", w);

		//w = w - R * A_D * w * exp(-(t - T_D) / TAU_D);
		//std::cout << "Weight before exp LTD = " << w << std::endl;
			w = w - w * R * A_D * exp(-(t - T_D) / TAU_D);

		//std::cout << "Weight after exp LTD = " << w << std::endl;
        if (w < 0)
            w = 0;

        //std::cout << "w = " << w << std::endl;
	}
}

//double PoolParallel::E_GABA(double t)
//{
//	return E_GABA_MATURE + (E_GABA_IMMATURE - E_GABA_MATURE) * exp(-t/T_GABA);
//}
void PoolParallel::gather_mature_data(std::vector<std::vector<double>>& average_dendritic_spike_time) 
{

    MPI_Status status;
    // Gather all data to master process
    int *spike_num_soma_local = new int[N_RA_local];
    int *spike_num_soma_global = new int[N_RA];
    int *spike_num_dend_local = new int[N_RA_local];
    int *spike_num_dend_global = new int[N_RA];
    int *spike_num_interneuron_local = new int[N_I_local];
    int *spike_num_interneuron_global = new int[N_I];

    int *recvcounts_RA = new int[MPI_size];
    int *displs_RA = new int[MPI_size];
	int *recvcounts_I = new int[MPI_size];
    int *displs_I = new int[MPI_size];


    if (MPI_rank == 0)
    {
        recvcounts_RA[0] = N_RA_sizes[0];
        displs_RA[0] = 0;

        for (int i = 1; i < MPI_size; i++)
        {
            recvcounts_RA[i] = N_RA_sizes[i];
            displs_RA[i] = displs_RA[i-1] + recvcounts_RA[i-1];
        }
		
		recvcounts_I[0] = N_I_sizes[0];
        displs_I[0] = 0;

        for (int i = 1; i < MPI_size; i++)
        {
            recvcounts_I[i] = N_I_sizes[i];
            displs_I[i] = displs_I[i-1] + recvcounts_I[i-1];
        }

    }

    for (int i = 0; i < N_RA_local; i++)
    {
        spike_num_soma_local[i] = spikes_in_trial_soma_local[i].size();
        spike_num_dend_local[i] = spikes_in_trial_dend_local[i].size();

        //printf("Rank = %d, supersyn_sizes_local[%d] = %d\n", MPI_rank, Id_RA_local[i], supersyn_sizes_local[i]);
        //printf("Rank = %d, syn_sizes_local[%d] = %d\n", MPI_rank, Id_RA_local[i], syn_sizes_local[i]);
        //printf("Rank = %d, spike_num_local[%d] = %d\n", MPI_rank, Id_RA_local[i], spike_num_local[i]);
    }

	for (int i = 0; i < N_I_local; i++)
		spike_num_interneuron_local[i] = (int) spikes_in_trial_interneuron_local[i].size();

	MPI_Gatherv(&spike_num_interneuron_local[0], N_I_local, MPI_INT,
        &spike_num_interneuron_global[0], recvcounts_I, displs_I, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gatherv(&spike_num_soma_local[0], N_RA_local, MPI_INT,
        &spike_num_soma_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gatherv(&spike_num_dend_local[0], N_RA_local, MPI_INT,
        &spike_num_dend_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);

    
	if (MPI_rank == 0)
    {
        for (int i = 0; i < N_RA; i++)
        {
            spikes_in_trial_soma_global[Id_RA_global[i]].resize(spike_num_soma_global[i]);
            spikes_in_trial_dend_global[Id_RA_global[i]].resize(spike_num_dend_global[i]);

            //printf("Master, spike_num_global[%d] = %d\n", i, spike_num_global[i]);
        }


        // Copy from master's local arrays
        for (int i = 0; i < N_RA_local; i++)
        {
            spikes_in_trial_soma_global[Id_RA_global[i]] = spikes_in_trial_soma_local[i];
            spikes_in_trial_dend_global[Id_RA_global[i]] = spikes_in_trial_dend_local[i];
        }

			
    // Gather from others
		int N = N_RA_sizes[0]; // number of RA neurons in the processes with lower rank

        for (int i = 1; i < MPI_size; i++)
        {

            for (int j = 0; j < N_RA_sizes[i]; j++)
            {
                int count;
				int index_in_global_array = N + j;
                int receive_index = Id_RA_global[index_in_global_array];
                
                if (spike_num_soma_global[index_in_global_array] != 0)
                {
                    MPI_Recv(&spikes_in_trial_soma_global[receive_index][0],
                        spike_num_soma_global[index_in_global_array], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv spikes in trial; from i = %d  count = %d\n", i, count);
                }

                if (spike_num_dend_global[index_in_global_array] != 0)
                {
                    MPI_Recv(&spikes_in_trial_dend_global[receive_index][0],
                        spike_num_dend_global[index_in_global_array], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv spikes in trial; from i = %d  count = %d\n", i, count);
                }

               // for (int k = 0; k < spikes_in_trial_global[N_RA_local + (i-1)*offset + j].size(); k++)
                   // printf("Master; spikes_in_trial_global[%d][%d] = %f\n", N_RA_local + (i-1)*offset + j, k,
                      //  spikes_in_trial_global[N_RA_local + (i-1)*offset + j][k]);
            }

			N += N_RA_sizes[i];

        }

    }

    else
    {
        for (int i = 0; i < N_RA_local; i++)
        {

            if (spike_num_soma_local[i] != 0)
                MPI_Send(&spikes_in_trial_soma_local[i][0],
                        spike_num_soma_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

            if (spike_num_dend_local[i] != 0)
                MPI_Send(&spikes_in_trial_dend_local[i][0],
                        spike_num_dend_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }

	// gather spikes of interneurons
	if (MPI_rank == 0)
    {
    
        for (int i = 0; i < N_I; i++)
			spikes_in_trial_interneuron_global[i].resize(spike_num_interneuron_global[i]);
        
	    for (int i = 0; i < N_I_local; i++)
			spikes_in_trial_interneuron_global[i] = spikes_in_trial_interneuron_local[i];
	
        int N = N_I_sizes[0]; // number of I neurons in the processes with lower rank

        for (int i = 1; i < MPI_size; i++)
        {
            for (int j = 0; j < N_I_sizes[i]; j++)
            {
                int count;
				int receive_index = N + j;

                if (spike_num_interneuron_global[receive_index] != 0)
                {
                    MPI_Recv(&spikes_in_trial_interneuron_global[receive_index][0],
                        spike_num_interneuron_global[receive_index], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv spikes in trial; from i = %d  count = %d\n", i, count);
                }

               // for (int k = 0; k < spikes_in_trial_global[N_RA_local + (i-1)*offset + j].size(); k++)
                   // printf("Master; spikes_in_trial_global[%d][%d] = %f\n", N_RA_local + (i-1)*offset + j, k,
                      //  spikes_in_trial_global[N_RA_local + (i-1)*offset + j][k]);
            }

			N += N_I_sizes[i];
        }
    }

    
    else
    {
        for (int i = 0; i < N_I_local; i++)
        {
            if (spike_num_interneuron_local[i] != 0)
                MPI_Send(&spikes_in_trial_interneuron_local[i][0],
                    spike_num_interneuron_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);


        }
    }
    
	// process dendritic spikes
	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_RA; i++)
		{
			if (spike_num_dend_global[i] > 0)
            {
                double average_spike_time = std::accumulate(spikes_in_trial_dend_global[i].begin(), spikes_in_trial_dend_global[i].end(), 0.0) / static_cast<double>(spike_num_dend_global[i]);

                //printf("Average dendritic spike time = %f\n", average_spike_time);
				average_dendritic_spike_time[i].push_back(average_spike_time);
            }
		}


	}


    delete [] recvcounts_RA;
    delete [] displs_RA;
    delete [] recvcounts_I;
    delete [] displs_I;
    delete [] spike_num_soma_local;
    delete [] spike_num_soma_global;
    delete [] spike_num_dend_local;
    delete [] spike_num_dend_global;
	delete [] spike_num_interneuron_local;
	delete [] spike_num_interneuron_global;


}

void PoolParallel::gather_data()
{
    MPI_Status status;
    // Gather all data to master process
    int *supersyn_sizes_global = new int[N_RA];
    int *supersyn_sizes_local = new int[N_RA_local];
    int *syn_sizes_global = new int[N_RA];
    int *syn_sizes_local = new int[N_RA_local];
    int *spike_num_soma_local = new int[N_RA_local];
    int *spike_num_soma_global = new int[N_RA];
    int *spike_num_dend_local = new int[N_RA_local];
    int *spike_num_dend_global = new int[N_RA];
    int *spike_num_interneuron_local = new int[N_I_local];
    int *spike_num_interneuron_global = new int[N_I];

    int *recvcounts_RA = new int[MPI_size];
    int *displs_RA = new int[MPI_size];
	int *recvcounts_I = new int[MPI_size];
    int *displs_I = new int[MPI_size];


    if (MPI_rank == 0)
    {
        recvcounts_RA[0] = N_RA_sizes[0];
        displs_RA[0] = 0;

        for (int i = 1; i < MPI_size; i++)
        {
            recvcounts_RA[i] = N_RA_sizes[i];
            displs_RA[i] = displs_RA[i-1] + recvcounts_RA[i-1];
        }
		
		recvcounts_I[0] = N_I_sizes[0];
        displs_I[0] = 0;

        for (int i = 1; i < MPI_size; i++)
        {
            recvcounts_I[i] = N_I_sizes[i];
            displs_I[i] = displs_I[i-1] + recvcounts_I[i-1];
        }

    }

    for (int i = 0; i < N_RA_local; i++)
    {
        supersyn_sizes_local[i] = supersynapses_local[i].size();
        syn_sizes_local[i] = active_synapses_local[i].size();
        spike_num_soma_local[i] = spikes_in_trial_soma_local[i].size();
        spike_num_dend_local[i] = spikes_in_trial_dend_local[i].size();

        //printf("Rank = %d, supersyn_sizes_local[%d] = %d\n", MPI_rank, Id_RA_local[i], supersyn_sizes_local[i]);
        //printf("Rank = %d, syn_sizes_local[%d] = %d\n", MPI_rank, Id_RA_local[i], syn_sizes_local[i]);
        //printf("Rank = %d, spike_num_local[%d] = %d\n", MPI_rank, Id_RA_local[i], spike_num_local[i]);
    }

    // all gatherv functions collect data from local processes. Data is a concatenated version of local 
    // data on the processes. Therefore neuronal ids may not be monotonic!!! Data is rearranged only 
    // writing to files

	for (int i = 0; i < N_I_local; i++)
		spike_num_interneuron_local[i] = (int) spikes_in_trial_interneuron_local[i].size();

    MPI_Gatherv(&supersyn_sizes_local[0], N_RA_local, MPI_INT,
        &supersyn_sizes_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gatherv(&syn_sizes_local[0], N_RA_local, MPI_INT,
        &syn_sizes_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);
	
	MPI_Gatherv(&spike_num_interneuron_local[0], N_I_local, MPI_INT,
        &spike_num_interneuron_global[0], recvcounts_I, displs_I, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gatherv(&spike_num_soma_local[0], N_RA_local, MPI_INT,
        &spike_num_soma_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gatherv(&spike_num_dend_local[0], N_RA_local, MPI_INT,
        &spike_num_dend_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);

	// gather maturation info
    MPI_Gatherv(&mature_local[0], N_RA_local, MPI_INT,
        &mature_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);
    
    MPI_Gatherv(&remodeled_local[0], N_RA_local, MPI_INT,
        &remodeled_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);
    
	MPI_Gatherv(&gaba_potential_local[0], N_RA_local, MPI_DOUBLE,
        &gaba_potential_global[0], recvcounts_RA, displs_RA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	MPI_Gatherv(&firing_rate_local[0], N_RA_local, MPI_DOUBLE,
        &firing_rate_global[0], recvcounts_RA, displs_RA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
    // Receive functions on the other hand rearrange data immediately. Thus, there is no need
    // to take special care while writing to files

    if (MPI_rank == 0)
    {
        for (int i = 0; i < N_RA; i++)
        {
            //printf("Master; supersyn_sizes_global[%d] = %d\n", i, supersyn_sizes_global[i]);
            supersynapses_global[Id_RA_global[i]].resize(supersyn_sizes_global[i]);
            active_synapses_global[Id_RA_global[i]].resize(syn_sizes_global[i]);
            spikes_in_trial_soma_global[Id_RA_global[i]].resize(spike_num_soma_global[i]);
            spikes_in_trial_dend_global[Id_RA_global[i]].resize(spike_num_dend_global[i]);

            //printf("Master, supersyn_sizes_global[%d] = %d\n", i, supersyn_sizes_global[i]);
            //printf("Master, syn_sizes_global[%d] = %d\n", i, syn_sizes_global[i]);
            //printf("Master, spike_num_global[%d] = %d\n", i, spike_num_global[i]);
        }


        // Copy from master's local arrays to master's global arrays
        // Note that data is copied according to neuronal ids stored in Id_RA_global
        // thus it is already ordered according to neuronal ids.
        for (int i = 0; i < N_RA_local; i++)
        {
            supersynapses_global[Id_RA_global[i]] = supersynapses_local[i];
            active_synapses_global[Id_RA_global[i]] = active_synapses_local[i];
            spikes_in_trial_soma_global[Id_RA_global[i]] = spikes_in_trial_soma_local[i];
            spikes_in_trial_dend_global[Id_RA_global[i]] = spikes_in_trial_dend_local[i];

            //for (int k = 0; k < active_supersynapses_global[i].size(); k++)
              //      printf("Master; active_supersynapses_global[%d][%d] = %u\n", i, k,
                //       active_supersynapses_global[i][k]);

            //for (int k = 0; k < active_synapses_global[i].size(); k++)
              //      printf("Master; active_synapses_global[%d][%d] = %u\n", i, k,
                //        active_synapses_global[i][k]);
        }

			
		for (int i = 0; i < N_RA_local; i++)
        {
            for (int j = 0; j < N_RA; j++)
            {
               weights_global[Id_RA_global[i]][j] = weights_local[i][j];

               //printf("Master; weights_global[%d][%d] = %f\n", i, j,
                 //   weights_global[i][j]);
            }

        }
    // Receive from others
		int N = N_RA_sizes[0]; // number of RA neurons in the processes with lower rank

        for (int i = 1; i < MPI_size; i++)
        {

            for (int j = 0; j < N_RA_sizes[i]; j++)
            {
                int count;
				int index_in_global_array = N + j;
                int receive_index = Id_RA_global[index_in_global_array];

                MPI_Recv(&weights_global[receive_index][0],
                                        N_RA, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_DOUBLE, &count);
                //printf("Recv weights; from i = %d  count = %d\n", i, count);
                if (supersyn_sizes_global[index_in_global_array] != 0)
                {
                    MPI_Recv(&supersynapses_global[receive_index][0],
                        supersyn_sizes_global[index_in_global_array], MPI_INT, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv supersynapses; from i = %d  count = %d\n", i, count);
                }

                if (syn_sizes_global[index_in_global_array] != 0)
                {
                    MPI_Recv(&active_synapses_global[receive_index][0],
                        syn_sizes_global[index_in_global_array], MPI_INT, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv synapses; from i = %d  count = %d\n", i, count);
                }

                if (spike_num_soma_global[index_in_global_array] != 0)
                {
                    MPI_Recv(&spikes_in_trial_soma_global[receive_index][0],
                        spike_num_soma_global[index_in_global_array], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv spikes in trial; from i = %d  count = %d\n", i, count);
                }

                if (spike_num_dend_global[index_in_global_array] != 0)
                {
                    MPI_Recv(&spikes_in_trial_dend_global[receive_index][0],
                        spike_num_dend_global[index_in_global_array], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv spikes in trial; from i = %d  count = %d\n", i, count);
                }
                //printf("Master;  N_RA_local + (i-1)*offset + j = %d\n", N_RA_local + (i-1)*offset + j);
                //printf("Master;  supersyn_sizes_global[N_RA_local + (i-1)*offset + j] = %d\n",
                 //   supersyn_sizes_global[N_RA_local + (i-1)*offset + j]);

                //for (int k = 0; k < N_RA; k++)
                  // printf("Master; weights_global[%d][%d] = %f\n", N_RA_local + (i-1)*offset + j, k,
                    //    weights_global[N_RA_local + (i-1)*offset + j][k]);

                //for (int k = 0; k < active_supersynapses_global[N_RA_local + (i-1)*offset + j].size(); k++)
                  //  printf("Master; active_supersynapses_global[%d][%d] = %u\n", N_RA_local + (i-1)*offset + j, k,
                    //   active_supersynapses_global[N_RA_local + (i-1)*offset + j][k]);

                //for (int k = 0; k < active_synapses_global[N_RA_local + (i-1)*offset + j].size(); k++)
                  //  printf("Master; active_synapses_global[%d][%d] = %u\n", N_RA_local + (i-1)*offset + j, k,
                    //    active_synapses_global[N_RA_local + (i-1)*offset + j][k]);

               // for (int k = 0; k < spikes_in_trial_global[N_RA_local + (i-1)*offset + j].size(); k++)
                   // printf("Master; spikes_in_trial_global[%d][%d] = %f\n", N_RA_local + (i-1)*offset + j, k,
                      //  spikes_in_trial_global[N_RA_local + (i-1)*offset + j][k]);
            }

			N += N_RA_sizes[i];

        }

    }

    else
    {
        for (int i = 0; i < N_RA_local; i++)
        {
            //for (int j = 0; j < N_RA; j++)
              //  printf("Rank %d; weights_local[%d][%d] = %f\n", MPI_rank, i, j, weights_local[i][j]);
            MPI_Send(&weights_local[i][0], N_RA, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

            if (supersyn_sizes_local[i] != 0)
                MPI_Send(&supersynapses_local[i][0],
                        supersyn_sizes_local[i], MPI_INT, 0, 0, MPI_COMM_WORLD);

            if (syn_sizes_local[i] != 0)
                MPI_Send(&active_synapses_local[i][0],
                        syn_sizes_local[i], MPI_INT, 0, 0, MPI_COMM_WORLD);

            if (spike_num_soma_local[i] != 0)
                MPI_Send(&spikes_in_trial_soma_local[i][0],
                        spike_num_soma_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

            if (spike_num_dend_local[i] != 0)
                MPI_Send(&spikes_in_trial_dend_local[i][0],
                        spike_num_dend_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            //printf("Rank %d;  supersyn_sizes_local[i] = %d\n", MPI_rank, supersyn_sizes_local[i]);
            //for (int k = 0; k < active_synapses_local[i].size(); k++)
             //   printf("Rank %d; active_synapses_local[%d][%d] = %u\n", MPI_rank,
              //      Id_RA_local[i], k, active_synapses_local[i][k]);

           // for (int k = 0; k < active_supersynapses_local[i].size(); k++)
               // printf("Rank %d; active_supersynapses_local[%d][%d] = %u\n", MPI_rank,
                 //   Id_RA_local[i], k, active_supersynapses_local[i][k]);


        }
    }

	// gather spikes of interneurons
	if (MPI_rank == 0)
    {
    
        for (int i = 0; i < N_I; i++)
			spikes_in_trial_interneuron_global[i].resize(spike_num_interneuron_global[i]);
        
	    for (int i = 0; i < N_I_local; i++)
			spikes_in_trial_interneuron_global[i] = spikes_in_trial_interneuron_local[i];
	
        int N = N_I_sizes[0]; // number of I neurons in the processes with lower rank

        for (int i = 1; i < MPI_size; i++)
        {
            for (int j = 0; j < N_I_sizes[i]; j++)
            {
                int count;
				int receive_index = N + j;

                if (spike_num_interneuron_global[receive_index] != 0)
                {
                    MPI_Recv(&spikes_in_trial_interneuron_global[receive_index][0],
                        spike_num_interneuron_global[receive_index], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv spikes in trial; from i = %d  count = %d\n", i, count);
                }

               // for (int k = 0; k < spikes_in_trial_global[N_RA_local + (i-1)*offset + j].size(); k++)
                   // printf("Master; spikes_in_trial_global[%d][%d] = %f\n", N_RA_local + (i-1)*offset + j, k,
                      //  spikes_in_trial_global[N_RA_local + (i-1)*offset + j][k]);
            }

			N += N_I_sizes[i];
        }
    }

    
    else
    {
        for (int i = 0; i < N_I_local; i++)
        {
            if (spike_num_interneuron_local[i] != 0)
                MPI_Send(&spikes_in_trial_interneuron_local[i][0],
                    spike_num_interneuron_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

           // for (int k = 0; k < active_supersynapses_local[i].size(); k++)
               // printf("Rank %d; active_supersynapses_local[%d][%d] = %u\n", MPI_rank,
                 //   Id_RA_local[i], k, active_supersynapses_local[i][k]);


        }
    }
    



    delete [] recvcounts_RA;
    delete [] displs_RA;
    delete [] recvcounts_I;
    delete [] displs_I;
    delete [] supersyn_sizes_global;
    delete [] supersyn_sizes_local;
    delete [] syn_sizes_global;
    delete [] syn_sizes_local;
    delete [] spike_num_soma_local;
    delete [] spike_num_soma_global;
    delete [] spike_num_dend_local;
    delete [] spike_num_dend_global;
	delete [] spike_num_interneuron_local;
	delete [] spike_num_interneuron_global;
}

void PoolParallel::get_neuronRA_location(int n, int* rank, int* shift)
{
	int i = 0;
	int N = 0;

	while (n >= N)
	{
		if (n < N + N_RA_sizes[i])
		{
			*rank = i;
			*shift = n - N;
		}
		
		N += N_RA_sizes[i];
		i++;
	}
}

void PoolParallel::get_neuronI_location(int n, int* rank, int* shift)
{
	int i = 0;
	int N = 0;

	while (n >= N)
	{
		if (n < N + N_I_sizes[i])
		{
			*rank = i;
			*shift = n - N;
		}
		
		N += N_I_sizes[i];
		i++;
	}
}

void PoolParallel::write_chain_test(int num_trials, std::vector<int>& num_trials_with_dend_spikes, std::vector<double>& average_num_dend_spikes_in_trials, 
                                    std::vector<double>& average_num_soma_spikes_in_trials, std::vector<double>& mean_burst_time, 
                                    std::vector<double>& std_burst_time, const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;

        out.open(filename, std::ios::out | std::ios::binary);

        // write number of neurons to each file

        out.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));
        out.write(reinterpret_cast<char *>(&num_trials), sizeof(num_trials));

        for (int i = 0; i < N_RA; i++)
        {
            double firing_robustness = num_trials_with_dend_spikes[i] / static_cast<double>(num_trials);

            out.write(reinterpret_cast<char *>(&firing_robustness), sizeof(firing_robustness));
            out.write(reinterpret_cast<char *>(&average_num_dend_spikes_in_trials[i]), sizeof(average_num_dend_spikes_in_trials[i]));
            out.write(reinterpret_cast<char *>(&average_num_soma_spikes_in_trials[i]), sizeof(average_num_soma_spikes_in_trials[i]));
            out.write(reinterpret_cast<char *>(&mean_burst_time[i]), sizeof(mean_burst_time[i]));
            out.write(reinterpret_cast<char *>(&std_burst_time[i]), sizeof(std_burst_time[i]));
        }
        out.close();
    }
}

void PoolParallel::write_active_synapses(const char* RA_RA)
{
    if (MPI_rank == 0)
    {
        std::ofstream out_RA_RA;

        int size;

        out_RA_RA.open(RA_RA, std::ios::out | std::ios::binary);

        // write number of neurons to each file

        out_RA_RA.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));

        for (int i = 0; i < N_RA; i++)
        {
            // write neuron's ID number
            out_RA_RA.write(reinterpret_cast<char *>(&i), sizeof(i));

            // write number of connections to RA neurons

            size = active_synapses_global[i].size();
            out_RA_RA.write(reinterpret_cast<char *>(&size), sizeof(size));

            for (int j = 0; j < size; j++)
            {
                int k = active_synapses_global[i][j];
                out_RA_RA.write(reinterpret_cast<char *>(&k), sizeof(k));
                out_RA_RA.write(reinterpret_cast<char *>(&weights_global[i][k]), sizeof(weights_global[i][k]));

            }

        }
        out_RA_RA.close();
    }
}

void PoolParallel::write_supersynapses(const char* RA_RA)
{
    if (MPI_rank == 0)
    {
        std::ofstream out_RA_RA;

        int size;

        out_RA_RA.open(RA_RA, std::ios::out | std::ios::binary);

        // write number of neurons to each file

        out_RA_RA.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));

        for (int i = 0; i < N_RA; i++)
        {
            // write neuron's ID number
            out_RA_RA.write(reinterpret_cast<char *>(&i), sizeof(i));

            // write number of connections to RA neurons

            size = supersynapses_global[i].size();
            out_RA_RA.write(reinterpret_cast<char *>(&size), sizeof(size));

            for (int j = 0; j < size; j++)
            {
                int k = supersynapses_global[i][j];
                out_RA_RA.write(reinterpret_cast<char *>(&k), sizeof(k));
                out_RA_RA.write(reinterpret_cast<char *>(&weights_global[i][k]), sizeof(weights_global[i][k]));

            }

        }
        out_RA_RA.close();
    }
}

void PoolParallel::write_weight_statistics(const char* filename)
{
	if (MPI_rank == 0)
	{
        std::ofstream out;

        if (trial_number == 0)
	    	out.open(filename, std::ios::binary | std::ios::out);
        else
	    	out.open(filename, std::ios::binary | std::ios::out | std::ios::app);

        out.write(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));

		double mean_synaptic_weight = 0.0;
        
        // calculate sum of all synaptic weights
        for (size_t i = 0; i < weights_global.size(); i++)
            mean_synaptic_weight += std::accumulate(weights_global[i].begin(), weights_global[i].end(), 0.0);

        // normalize sum to get mean synaptic weight
        if (static_cast<int>(weights_global.size() > 0))
        {
            if (static_cast<int>(weights_global[0].size() > 0))
                mean_synaptic_weight /= (static_cast<double>(weights_global.size()) * static_cast<double>(weights_global[0].size()));
            else
                std::cerr << "Matrix weights_global has no columns!" << std::endl;
        }
        else
            std::cerr << "Matrix weights_global has no rows!" << std::endl;
        
        double std_synaptic_weight = 0; // standard deviation of synaptic weights

        for (size_t i = 0; i < weights_global.size(); i++)
        {
            std::for_each(weights_global[i].begin(), weights_global[i].end(), [&std_synaptic_weight, mean_synaptic_weight](const double w)
            {
                std_synaptic_weight += (w - mean_synaptic_weight) * (w - mean_synaptic_weight);
            });
        }
        
        // if there are more than two weights in the weight matrix
        if (static_cast<int>(weights_global.size()) * static_cast<int>(weights_global[0].size()) > 1)
            std_synaptic_weight = sqrt(std_synaptic_weight / (static_cast<double>(weights_global.size()) * static_cast<double>(weights_global[0].size()) - 1));
        else
            std::cerr << "Matrix weights_global has less than two elements!" << std::endl;

        out.write(reinterpret_cast<char *>(&mean_synaptic_weight), sizeof(mean_synaptic_weight));
        out.write(reinterpret_cast<char *>(&std_synaptic_weight), sizeof(std_synaptic_weight));

        out.close();
	}
}

void PoolParallel::write_weights(const char* filename)
{
	if (MPI_rank == 0)
	{
        std::ofstream output;

        output.open(filename, std::ios::out | std::ios::binary);

        output.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));
        output.write(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));

        for (int i = 0; i < N_RA; i++)
        {
            for (int j = 0; j < N_RA; j++)
            {
                //printf("weigths[%d][%d] = %1.10f\n", i, j, weights[i][j]);
                output.write(reinterpret_cast<char *>(&weights_global[i][j]), sizeof(weights_global[i][j]));

            }
        }
        output.close();
	}
}

void PoolParallel::write_num_synapses(const char* fileSynapses)
{
	if (MPI_rank == 0)
	{
		std::ofstream out_synapses;
   
        if (trial_number == 0)
	    	out_synapses.open(fileSynapses, std::ios::binary | std::ios::out);
        else
	    	out_synapses.open(fileSynapses, std::ios::binary | std::ios::out | std::ios::app);

		// count active synapses and supersynapses
		int active = 0;
		int super = 0;

		for (int i = 0; i < N_RA; i++)
		{
			active += (int) active_synapses_global[i].size(); 
			super += (int) supersynapses_global[i].size(); 

		}
		
        out_synapses.write(reinterpret_cast<char*>(&trial_number), sizeof(trial_number)); // write current trial numbe
		out_synapses.write(reinterpret_cast<char*>(&active), sizeof(active)); // write number of active synapses
		out_synapses.write(reinterpret_cast<char*>(&super), sizeof(super)); // write number of supersynapses
	
        out_synapses.close();
	}

}

void PoolParallel::write_all_coordinates()
{
    if (MPI_rank == 0)
    {
        std::string fileRAxy = outputNetworkDir + "RA_xy.bin";
        std::string fileIxy = outputNetworkDir + "I_xy.bin";

        this->write_coordinates_RA(fileRAxy.c_str());
        this->write_coordinates_I(fileIxy.c_str());
    }
}

void PoolParallel::write_coordinates_RA(const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;

        // open files
        out.open(filename, std::ios::binary | std::ios::out);

        // write number of neurons
        out.write(reinterpret_cast<char*>(&N_RA), sizeof(N_RA));

        // write coordinates
        for (int i = 0; i < N_RA; i++)
        {
            out.write(reinterpret_cast<char*>(&xx_RA[i]), sizeof(xx_RA[i]));
            out.write(reinterpret_cast<char*>(&yy_RA[i]), sizeof(yy_RA[i]));
        }

        // close files
        out.close();
	}
}

void PoolParallel::write_coordinates_I(const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;

        // open files
        out.open(filename, std::ios::binary | std::ios::out);

        // write number of neurons
        out.write(reinterpret_cast<char*>(&N_I), sizeof(N_I));

        // write coordinates
        for (int i = 0; i < N_I; i++)
        {
            out.write(reinterpret_cast<char*>(&xx_I[i]), sizeof(xx_I[i]));
            out.write(reinterpret_cast<char*>(&yy_I[i]), sizeof(yy_I[i]));
        }

        // close files
        out.close();
	}
}

void PoolParallel::write_soma_spike_times(const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;

        out.open(filename, std::ios::out | std::ios::binary );
        out.write(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));
        out.write(reinterpret_cast<char *>(&internal_time), sizeof(internal_time));
        out.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));

        //for (int i = 0; i < N_RA; i++)
        //    printf("spike_times[%d] = %f\n", i, spike_times[i]);
        // write spike times
        for (int i = 0; i < N_RA; i++)
        {
            int spike_array_size = spikes_in_trial_soma_global[i].size();
            //printf("Neuron %d; number of somatic spikes in trial: %d\n", i, spike_array_size);
            out.write(reinterpret_cast<char *>(&spike_array_size), sizeof(int));

            for (int j = 0; j < spike_array_size; j++)
	        {
       		//printf("relative_spike_time_soma = %f\n", relative_spike_time);
                out.write(reinterpret_cast<char *>(&spikes_in_trial_soma_global[i][j]), sizeof(spikes_in_trial_soma_global[i][j]));
	        }
        }
        //out.write(reinterpret_cast<char *>(spike_times), N_RA*sizeof(double));

        out.close();
    }
}

void PoolParallel::write_dend_spike_times(const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;

        out.open(filename, std::ios::out | std::ios::binary );
        out.write(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));
        out.write(reinterpret_cast<char *>(&internal_time), sizeof(internal_time));
        out.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));

        //for (int i = 0; i < N_RA; i++)
        //    printf("spike_times[%d] = %f\n", i, spike_times[i]);
        // write spike times
        for (int i = 0; i < N_RA; i++)
        {
            int spike_array_size = spikes_in_trial_dend_global[i].size();
            //printf("Neuron %d; number of dendritic spikes in trial: %d\n", i, spike_array_size);
            out.write(reinterpret_cast<char *>(&spike_array_size), sizeof(int));
	    
            for (int j = 0; j < spike_array_size; j++)
	    {
                //out.write(reinterpret_cast<char *>(&spikes_in_trial_dend_global[i][j]), sizeof(double));
        	out.write(reinterpret_cast<char *>(&spikes_in_trial_dend_global[i][j]), sizeof(spikes_in_trial_dend_global[i][j]));
			//printf("Neuron %d; relative spike time = %f\n", i, relative_spike_time);
	    }
	}
        //out.write(reinterpret_cast<char *>(spike_times), N_RA*sizeof(double));

        out.close();
    }
}

void PoolParallel::write_interneuron_spike_times(const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;

        out.open(filename, std::ios::out | std::ios::binary );
        out.write(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));
        out.write(reinterpret_cast<char *>(&internal_time), sizeof(internal_time));
        out.write(reinterpret_cast<char *>(&N_I), sizeof(N_I));

        //for (int i = 0; i < N_RA; i++)
        //    printf("spike_times[%d] = %f\n", i, spike_times[i]);
        // write spike times
        for (int i = 0; i < N_I; i++)
        {
            int spike_array_size = spikes_in_trial_interneuron_global[i].size();
            //printf("Neuron %d; number of dendritic spikes in trial: %d\n", i, spike_array_size);
            out.write(reinterpret_cast<char *>(&spike_array_size), sizeof(int));
	    
            for (int j = 0; j < spike_array_size; j++)
	    {
                //out.write(reinterpret_cast<char *>(&spikes_in_trial_dend_global[i][j]), sizeof(double));
        	out.write(reinterpret_cast<char *>(&spikes_in_trial_interneuron_global[i][j]), sizeof(spikes_in_trial_interneuron_global[i][j]));
			//printf("Neuron %d; relative spike time = %f\n", i, relative_spike_time);
	    }
	}
        //out.write(reinterpret_cast<char *>(spike_times), N_RA*sizeof(double));

        out.close();
    }
}

void PoolParallel::write_invariable_synapses(const char* RA_I, const char* I_RA)
{
    if (MPI_rank == 0)
    {
        std::ofstream out_RA_I, out_I_RA;
        int size;
        // open files
        out_RA_I.open(RA_I, std::ios::binary | std::ios::out);
        out_I_RA.open(I_RA, std::ios::binary | std::ios::out);

        // write number of neurons
        out_RA_I.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));
        out_I_RA.write(reinterpret_cast<char *>(&N_I), sizeof(N_I));

        // write connections from RA to I
        for (int i = 0; i < N_RA; i++)
        {
            size = syn_ID_RA_I_global[i].size();

            out_RA_I.write(reinterpret_cast<char *>(&i), sizeof(i));
            out_RA_I.write(reinterpret_cast<char *>(&size), sizeof(size)); // write neuron's ID

            for (int j = 0; j < size; j++)
            {
                int k = syn_ID_RA_I_global[i][j];
                double G = weights_RA_I_global[i][j];

                out_RA_I.write(reinterpret_cast<char *>(&k), sizeof(k));
                out_RA_I.write(reinterpret_cast<char *>(&G), sizeof(G));

            }

        }

        // write connections from I to RA
        for (int i = 0; i < N_I; i++)
        {
            out_I_RA.write(reinterpret_cast<char *>(&i), sizeof(i)); // write neuron's ID number

            size = syn_ID_I_RA_global[i].size();
            out_I_RA.write(reinterpret_cast<char *>(&size), sizeof(size)); // write number of targets a neuron has
            for (int j = 0; j < size; j++)
            {
                    int k = syn_ID_I_RA_global[i][j];
                    double G = weights_I_RA_global[i][j];
                    out_I_RA.write(reinterpret_cast<char *>(&k), sizeof(k)); // write targets ID
                    out_I_RA.write(reinterpret_cast<char *>(&G), sizeof(G)); // write targets conductance

            }
        }
        // close files
        out_I_RA.close();
        out_RA_I.close();
	}
}

void PoolParallel::write_pajek_fixed(const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;
        out.open(filename, std::ios::out);

        out << "*Vertices " << N_RA + N_I << "\n";

        for (int i = 0; i < N_RA; i++)
        {

            out << i + 1 << " \"" << i << "\" " << xx_RA[i]/SIDE << " " << yy_RA[i]/SIDE << " ic Green\n";
        }

        for (int i = 0; i < N_I; i++)
        {

            out << i + N_RA + 1 << " \"" << i + N_RA << "\" " << xx_I[i]/SIDE << " " << yy_I[i]/SIDE << " ic Red\n";
        }

        out << "*Arcs\n";
        for (int i = 0; i < N_RA; i++)
        {
            for (int j = 0; j < syn_ID_RA_I_global[i].size(); j++)
                {
                    int syn_ID = syn_ID_RA_I_global[i][j] + N_RA;
                    out << i + 1 << " " << syn_ID + 1 << " " << weights_RA_I_global[i][j] << " c Green\n";
                }
        }
        for (int i = 0; i < N_I; i++)
        {
            for (int j = 0; j < syn_ID_I_RA_global[i].size(); j++)
                {
                    int syn_ID = syn_ID_I_RA_global[i][j];
                    out << i + N_RA + 1 << " " << syn_ID + 1 << " " << weights_I_RA_global[i][j] << " c Red\n";
                }
        }
    }
}

void PoolParallel::write_RA(const char* filename, int n)
{
	if (n >= N_RA)
    {

        if (MPI_rank == 0)
        {
			printf("Selected neuron ID doesn't exist in the pool! Instead wrote neuron 0 to the file.\n");
            HVCRA_local[0].writeToFile(filename);
    	}
	}
    else
    {
		int N = 0; // number of neurons in all previous processes
        int i = 0;

		while (n >= N)
		{
			if (n < N + N_RA_sizes[i])
			{
				if (MPI_rank == i)
				{
					HVCRA_local[n-N].writeToFile(filename);
					//printf("My rank is %d; Writing neuron %d for n-N = %d\n", MPI_rank, n, n-N);
				}
			}

			N += N_RA_sizes[i];
			i++;
		}
		
    }
}


void PoolParallel::write_I(const char* filename, int n)
{
    if (n >= N_I)
    {
        if (MPI_rank == 0)
        {
			printf("Selected neuron %d doesn't exist in the pool! Instead wrote neuron 0 to the file.\n",n);
            HVCI_local[0].writeToFile(filename);
    	}
	}
    else
    {
		int N = 0; // number of neurons in all previous processes
        int i = 0;

		while (n >= N)
		{
			if (n < N + N_I_sizes[i])
			{
				//printf("My rank is %d; n = %d; N = %d\n", MPI_rank, n, N);
				if (MPI_rank == i)
				{
					//printf("My rank is %d; Writing neuron %d\n", MPI_rank, n);
					HVCI_local[n-N].writeToFile(filename);

				}
			}

			N += N_I_sizes[i];
			i++;
		}
		
    }
}

void PoolParallel::write_maturation_info(const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;
      
        out.open(filename, std::ios::out | std::ios::binary );
        out.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));

      
        out.write(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));

        for (int i = 0; i < N_RA; i++)
        {
            // find position in global id array
            std::vector<int>::iterator pos = std::find(Id_RA_global.begin(), Id_RA_global.end(), i);

            if (pos == Id_RA_global.end())
                std::cerr << "Neuronal id " << i << " is not found in Id_RA_global in write_maturation_info" << std::endl; 
            else
            {
                int ind = std::distance(Id_RA_global.begin(), pos);
 
                //std::cout << "i = " << i << " ; ind = " << ind << std::endl;

                out.write(reinterpret_cast<char *>(&gaba_potential_global[ind]), sizeof(gaba_potential_global[ind]));
                out.write(reinterpret_cast<char *>(&firing_rate_global[ind]), sizeof(firing_rate_global[ind]));
                out.write(reinterpret_cast<char *>(&remodeled_global[ind]), sizeof(remodeled_global[ind]));
                out.write(reinterpret_cast<char *>(&mature_global[ind]), sizeof(mature_global[ind]));
	        }
	    }
		out.close();
	}

}

void PoolParallel::write_maturation_time_sequence(const std::vector<int>& neurons, const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;
      
        if (trial_number == 0)
        {
            out.open(filename, std::ios::out | std::ios::binary );
            
            int size = static_cast<int>(neurons.size());
            
            out.write(reinterpret_cast<char *>(&size), sizeof(size));

            for (size_t i = 0; i < neurons.size(); i++)
                out.write(reinterpret_cast<const char *>(&neurons[i]), sizeof(neurons[i]));

        }
        else
            out.open(filename, std::ios::out | std::ios::binary | std::ios::app);

        out.write(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));

        for (size_t i = 0; i < neurons.size(); i++)
        {
            // find position in global id array
            std::vector<int>::iterator pos = std::find(Id_RA_global.begin(), Id_RA_global.end(), neurons[i]);

            if (pos == Id_RA_global.end())
                std::cerr << "Neuronal id " << neurons[i] << " is not found in Id_RA_global in write_maturation_time_sequence" << std::endl; 
            else
            {
                int ind = std::distance(Id_RA_global.begin(), pos);
 
                out.write(reinterpret_cast<char *>(&gaba_potential_global[ind]), sizeof(gaba_potential_global[ind]));
                out.write(reinterpret_cast<char *>(&firing_rate_global[ind]), sizeof(firing_rate_global[ind]));
                out.write(reinterpret_cast<char *>(&remodeled_global[ind]), sizeof(remodeled_global[ind]));
                out.write(reinterpret_cast<char *>(&mature_global[ind]), sizeof(mature_global[ind]));
	        }
	    
	    }
		out.close();
	}
}

void PoolParallel::write_weights_time_sequence_from_source_to_target(const std::vector<int>& source, const std::vector<int>& target, const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;
      
        if (trial_number == 0)
        {
            out.open(filename, std::ios::out | std::ios::binary );
            
            int size_source = static_cast<int>(source.size());
            int size_target = static_cast<int>(target.size());
            
            out.write(reinterpret_cast<char *>(&size_source), sizeof(size_source));
            out.write(reinterpret_cast<char *>(&size_target), sizeof(size_target));

            for (size_t i = 0; i < source.size(); i++)
                out.write(reinterpret_cast<const char *>(&source[i]), sizeof(source[i]));

            for (size_t i = 0; i < target.size(); i++)
                out.write(reinterpret_cast<const char *>(&target[i]), sizeof(target[i]));
        }
        else
            out.open(filename, std::ios::out | std::ios::binary | std::ios::app);

        out.write(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));

        for (size_t i = 0; i < source.size(); i++)
        {
            for (size_t j = 0; j < target.size(); j++)
            {
                out.write(reinterpret_cast<char *>(&weights_global[source[i]][target[j]]), sizeof(double));
	        }
	    }
		out.close();
	}

}

void PoolParallel::write_global_index_array(const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;
      
        out.open(filename, std::ios::out | std::ios::binary );

        out.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));

        for (size_t i = 0; i < Id_RA_global.size(); i++)
            out.write(reinterpret_cast<char *>(&Id_RA_global[i]), sizeof(Id_RA_global[i]));
		
        out.close();
	}
}

void PoolParallel::write_replaced_neurons(std::vector<int>& global_id_2replace, const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;
      
        out.open(filename, std::ios::out | std::ios::binary);

        int num_replaced = static_cast<int>(global_id_2replace.size()); // number of replaced neurons

        // write number of replaced neurons
        out.write(reinterpret_cast<char *>(&num_replaced), sizeof(num_replaced));

        for (int i = 0; i < num_replaced; i++)
            out.write(reinterpret_cast<char *>(&global_id_2replace[i]), sizeof(global_id_2replace[i]));
		
        out.close();
	}
}
