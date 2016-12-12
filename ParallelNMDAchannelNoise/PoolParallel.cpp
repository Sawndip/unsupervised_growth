#include "PoolParallel.h"
#include "HHI_final_pool.h"
#include "HH2_final_pool.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <functional>
#include "training_current.h"

using namespace std::placeholders;

PoolParallel::PoolParallel(double a, double s_rai, double b, double s_ira, double network_update,
						double Ei, double beta, double beta_s, double Tp, double Td, double tauP, double tauD, double Ap, double Ad, double Ap_super, 
						double Ad_super, double f0, double activation, double super_threshold, 
                        double Gmax, double gaba_down, int N_ra, int Ni, int N_ss, int N_tr) : A_RA2I(a), 
						SIGMA_RA2I(s_rai), B_I2RA(b), SIGMA_I2RA(s_ira), network_update_frequency(network_update),
						E_GABA_IMMATURE(Ei), BETA(beta), BETA_SUPERSYNAPSE(beta_s), A_P(Ap), A_D(Ad), T_P(Tp), T_D(Td), TAU_P(tauP), 
						TAU_D(tauD), A_P_SUPER(Ap_super),
						A_D_SUPER(Ad_super), F_0(f0), ACTIVATION(activation), SUPERSYNAPSE_THRESHOLD(super_threshold), 
						GABA_DOWN(gaba_down), G_MAX(Gmax),
				        N_RA(N_ra), N_I(Ni),  Nss(N_ss), N_TR(N_tr)
{

    MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);

    // white noise
    double mu_soma = 30; // dc component of white noise to somatic compartment
    double sigma_soma = 100; // variance of white noise to somatic compartment
    double mu_dend = 50; // dc component of white noise to dendritic compartment
    double sigma_dend = 150; // variance of white noise to dendritic compartment

	double delay = 0;
	int num_related_spikes = 0;

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

	HVCRA_local = new HH2_final_pool[N_RA_local];
	HVCI_local = new HHI_final_pool[N_I_local];

    Id_RA_local.resize(N_RA_local);
    Id_I_local.resize(N_I_local);

    	//active_local = new bool*[N_RA];
	//active_supersynapses_local = new bool*[N_RA_local];
	
	weights_RA_I_local = new std::vector<double>[N_RA_local];
	weights_I_RA_local = new std::vector<double>[N_I_local];
	syn_ID_RA_I_local = new std::vector<int>[N_RA_local];
	syn_ID_I_RA_local = new std::vector<int>[N_I_local];

    weights_global = new double*[N_RA];
    weights_RA_I_global = new std::vector<double>[N_RA];
    weights_I_RA_global = new std::vector<double>[N_I];
    syn_ID_RA_I_global = new std::vector<int>[N_RA];
    syn_ID_I_RA_global = new std::vector<int>[N_I];

    active_global = new bool*[N_RA];
    supersynapses_global = new bool*[N_RA];
	active_supersynapses_global = new std::vector<int>[N_RA];
	active_synapses_global = new std::vector<int>[N_RA];

    weights_local = new double*[N_RA_local];
    active_local = new bool*[N_RA_local];
    supersynapses_local = new bool*[N_RA_local];
    active_supersynapses_local = new std::vector<int>[N_RA_local];
	active_synapses_local = new std::vector<int>[N_RA_local];

	update_Ge_RA_local.resize(N_RA);
    update_Gi_RA_local.resize(N_RA);
    update_Ge_I_local.resize(N_I);

    update_Ge_RA_global.resize(N_RA);
    update_Gi_RA_global.resize(N_RA);
    update_Ge_I_global.resize(N_I);
	
	last_soma_spikes_local = new std::deque<double>[N_RA];
	last_soma_spikes_global = new std::deque<double>[N_RA];

    spikes_in_trial_soma_global = new std::vector<double>[N_RA];
    spike_times_dend_global = new double[N_RA];
    spikes_in_trial_dend_global = new std::vector<double>[N_RA];
    spikes_in_trial_interneuron_global = new std::vector<double>[N_I];

    spikes_in_trial_soma_local = new std::vector<double>[N_RA_local];
    spikes_in_trial_dend_local = new std::vector<double>[N_RA_local];
    spikes_in_trial_interneuron_local = new std::vector<double>[N_I_local];

	remodeled_local.resize(N_RA_local);

	mature_local.resize(N_RA_local);
	mature_global.resize(N_RA);

	maturation_triggered_local.resize(N_RA_local);
	maturation_triggered_global.resize(N_RA);

	gaba_potential_local.resize(N_RA_local);
	gaba_potential_global.resize(N_RA);


	// initialize circular buffers for somatic spikes
	num_bursts_in_recent_trials.resize(N_RA_local);
	std::fill(num_bursts_in_recent_trials.begin(), num_bursts_in_recent_trials.end(), intBuffer(RATE_WINDOW));

	for (size_t i = 0; i < num_bursts_in_recent_trials.size(); i++)
		for (int j = 0; j < RATE_WINDOW; j++)
			num_bursts_in_recent_trials[i].push_back(0);

	for (int i = 0; i < N_RA_local; i++)
	{
        weights_local[i] = new double[N_RA];
        active_local[i] = new bool[N_RA];
        supersynapses_local[i] = new bool[N_RA];
        remodeled_local[i] = false;
	
	    for (int j = 0; j < N_RA; j++)
        {
            weights_local[i][j] = 0.0;
            active_local[i][j] = false;
            supersynapses_local[i][j] = false;
        }

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


	for (int i = 0; i < N_RA; i++)
	{
		active_global[i] = new bool[N_RA];
		supersynapses_global[i] = new bool[N_RA];
		weights_global[i] = new double[N_RA];

	}
	
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
		maturation_triggered_local[i] = 0;
	}

	for (int i = 0; i < N_RA; i++)
	{
		for (int j = 0; j < N_RA; j++)
		{
			active_global[i][j] = false;
			supersynapses_global[i][j] = false;
			weights_global[i][j] = 0.0;

		}
        spike_times_dend_global[i] = -100.0;

		
	}
}

PoolParallel::~PoolParallel()
{
	delete[] HVCRA_local;
	delete[] HVCI_local;

	for (int i = 0; i < N_RA; i++)
	{
		delete[] active_global[i];
		delete[] supersynapses_global[i];
		delete[] weights_global[i];
	}

	delete[] spike_times_dend_global;

	delete[] active_global;
	delete[] supersynapses_global;
	delete[] weights_global;

	delete[] weights_RA_I_local;
	delete[] weights_I_RA_local;
	delete[] syn_ID_RA_I_local;
	delete[] syn_ID_I_RA_local;

	delete[] weights_RA_I_global;
	delete[] weights_I_RA_global;
	delete[] syn_ID_RA_I_global;
	delete[] syn_ID_I_RA_global;

//	MPI_Finalize();
}


const double PoolParallel::DEMOTIVATION_WINDOW = 200; // demotivation window in ms
// coordinates
const double PoolParallel::MIN_INTERNEURON_DISTANCE = 0.01; // minimum distance between neurons

const double PoolParallel::SIDE = 100; // length of HVC side


// developmental GABA switch
const double PoolParallel::E_GABA_MATURE = -80;

const double PoolParallel::BURST_RATE_THRESHOLD = 0.75; // threshold fraction of trials. if neuron exceeds the threshold, its maturation program is triggered
const int PoolParallel::RATE_WINDOW = 10; // window in which firing rate is averaged (num of trials)

// constants for STDP-rules
const double PoolParallel::G_P = 0.1;

const double PoolParallel::R = 1;

const double PoolParallel::STDP_WINDOW = 100;

void PoolParallel::initialize_generator()
{
    // prepare seeds and initialize each generator with different seed
    srand((unsigned int) time(NULL));

    //unsigned seed = rand();
    //unsigned seed = rand() + MPI_rank;

    unsigned seed = 15000 + MPI_rank*1000;

    generator.set_seed(seed);
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
					xx = (double) (i+1) * SIDE / (sqrt(N_I)+1)+ generator.random(SIDE / (5*sqrt(N_I)));
					yy = (double) (k+1) * SIDE / (sqrt(N_I)+1) + generator.random(SIDE / (5*sqrt(N_I)));

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
                xx = 0.5*SIDE / (sqrt(N_I)+1) + generator.random(SIDE - 0.5*SIDE / (sqrt(N_I)+1));
                yy = 0.5*SIDE / (sqrt(N_I)+1) + generator.random(SIDE - 0.5*SIDE / (sqrt(N_I)+1));

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

void PoolParallel::read_from_file(const char* RA_xy, const char* I_xy, const char* RA_RA_all, const char* RA_RA_active, const char* RA_RA_super, const char* RA_I, const char * I_RA, const char* maturation, const char* timeInfo)
{
	if (MPI_rank == 0)
	{
	
		std::ifstream inp_RA_xy, inp_I_xy, inp_RA_RA, inp_RA_RA_super, inp_RA_I, inp_I_RA, inp_maturation, inp_timeInfo;
	
		// input files with coordinates of neurons in the pool
		inp_RA_xy.open(RA_xy, std::ios::binary | std::ios::in);
		inp_I_xy.open(I_xy, std::ios::binary | std::ios::in);

		// input files with fixed connections
		inp_RA_I.open(RA_I, std::ios::binary | std::ios::in);
		inp_I_RA.open(I_RA, std::ios::binary | std::ios::in);
	
		inp_RA_xy.read(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));
		inp_I_xy.read(reinterpret_cast<char *>(&N_I), sizeof(N_I));

		inp_RA_I.read(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));
		inp_I_RA.read(reinterpret_cast<char *>(&N_I), sizeof(N_I));

		printf("Number of RA neurons in the pool: %d\n", N_RA);
		printf("Number of I neurons in the pool: %d\n", N_I);
	

		// read coordinates of neurons in the pool
		// read RA coordinates
		for (int i = 0; i < N_RA; i++)
		{
			double x, y;
			
			inp_RA_xy.read(reinterpret_cast<char *>(&x), sizeof(x));
			inp_RA_xy.read(reinterpret_cast<char *>(&y), sizeof(y));

			//printf("RA neuron %d has coordinates (%f, %f)\n", i, x, y);
			
			xx_RA.push_back(x);
			yy_RA.push_back(y);
		}

		inp_RA_xy.close();
		
		// read I coordinates
		for (int i = 0; i < N_I; i++)
		{
			double x, y;
			
			inp_I_xy.read(reinterpret_cast<char *>(&x), sizeof(x));
			inp_I_xy.read(reinterpret_cast<char *>(&y), sizeof(y));
			
			//printf("I neuron %d has coordinates (%f, %f)\n", i, x, y);
			
			xx_I.push_back(x);
			yy_I.push_back(y);
		}

		inp_I_xy.close();

		// read RA to RA connections
		//num_supersynapses = new int[N_RA];
		//num_active = new int[N_RA];

		std::vector<double>* supersynapses_G = new std::vector<double>[N_RA];
		std::vector<double>* active_G = new std::vector<double>[N_RA];

		// read all RA to RA connections
	
		//this->read_all_connections_from_net(RA_RA_all, &weights_global);	
	    this->read_all_connections(RA_RA_all, weights_global);

        //for (int i = 0; i < N_RA; i++)
            //printf("weights[0][%d] = %f\n", i, weights_global[0][i]);

        

		// read RA to RA supersynapses
		this->read_connections_from_net(RA_RA_super, &active_supersynapses_global, &supersynapses_G);	
		this->read_connections_from_net(RA_RA_active, &active_synapses_global, &active_G);
	
		//for (int i = 0; i < N_RA; i++)
		//{
		//	printf("Neuron %d has %d supersynapses\n", i, (int) active_supersynapses_global[i].size());

		//	for (int j = 0; j < (int) active_supersynapses_global[i].size(); j++)
		//	{
		//		printf("target id: %d\tsynaptic weight: %f\n", active_supersynapses_global[i][j], supersynapses_G[i][j]);
		//	}

		//}
	
		//for (int i = 0; i < N_RA; i++)
		//{
		//	for (int j = 0; j < N_RA; j++)
		//	{
		//		printf("%d to %d %e\n", i, j, weights_global[i][j]);

		//	}

		//}
		//for (int i = 0; i < N_RA; i++)
		//{
		//	num_supersynapses[i] = super_targets_id[i].size();
		//}


	// read RA to I connections

		for (int i = 0; i < N_RA; i++)
		{
			int id;
			int num_targets;
	
			inp_RA_I.read(reinterpret_cast<char *>(&id), sizeof(id)); // neuron's id
			inp_RA_I.read(reinterpret_cast<char *>(&num_targets), sizeof(num_targets));
			//printf("Neuron %d has %d targets:\n", id-1, num_targets);
	
			//printf("Neuron %d has %d inhibitory targets:\n", id-1, num_targets);
			
			for (int j = 0; j < num_targets; j++)
			{
				int target; // id of target
				double G; // synaptic conductance
	
				inp_RA_I.read(reinterpret_cast<char *>(&target), sizeof(target));
				
				inp_RA_I.read(reinterpret_cast<char *>(&G), sizeof(G));
				
				weights_RA_I_global[i].push_back(G);
				syn_ID_RA_I_global[i].push_back(target);
				//printf("target %d with weight %f\n", target, G);
			}
		}
	
		// read I to RA connections
		for (int i = 0; i < N_I; i++)
		{
			int id;
			int num_targets;
	
			inp_I_RA.read(reinterpret_cast<char *>(&id), sizeof(id)); // neuron's id
			inp_I_RA.read(reinterpret_cast<char *>(&num_targets), sizeof(num_targets));
			//printf("Neuron %d has %d targets:\n", id-1, num_targets);
	
			//printf("I neuron %d has %d excitatory targets:\n", id-1, num_targets);
			
			for (int j = 0; j < num_targets; j++)
			{
				int target; // id of target
				double G; // synaptic conductance
	
				inp_I_RA.read(reinterpret_cast<char *>(&target), sizeof(target));
				
				inp_I_RA.read(reinterpret_cast<char *>(&G), sizeof(G));
				
				weights_I_RA_global[i].push_back(G);
				syn_ID_I_RA_global[i].push_back(target);
				//printf("target %d with weight %f\n", target, G);
			}
		}
		
		// read maturation info

		inp_maturation.open(maturation, std::ios::binary | std::ios::in);
		inp_maturation.read(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));

		for (int i = 0; i < N_RA; i++)
		{	
			inp_maturation.read(reinterpret_cast<char *>(&gaba_potential_global[i]), sizeof(gaba_potential_global[i]));
			inp_maturation.read(reinterpret_cast<char *>(&maturation_triggered_global[i]), sizeof(maturation_triggered_global[i]));
			inp_maturation.read(reinterpret_cast<char *>(&mature_global[i]), sizeof(mature_global[i]));
			
			//if (mature_global[i] > 0)
				//printf("Neuron %d is mature\n", i);
			//else
				//printf("Neuron %d is immature\n", i);
		}

		// read time information

        inp_timeInfo.open(timeInfo, std::ios::binary | std::ios::in);

		inp_timeInfo.read(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));
		inp_timeInfo.read(reinterpret_cast<char *>(&internal_time), sizeof(internal_time));
		inp_timeInfo.read(reinterpret_cast<char *>(&network_time), sizeof(network_time));

        printf("trial number = %d\n", trial_number);
        printf("internal_time = %f\n", internal_time);
        printf("network_time = %f\n", network_time);

		inp_RA_I.close();
		inp_I_RA.close();
		inp_maturation.close();
        inp_timeInfo.close();
	}
}

void PoolParallel::print_simulation_parameters()
{
	if (MPI_rank == 0)
	{
		printf("\nSpatial distribution\n");
		printf("A_RA2I = %f\n", A_RA2I);
		printf("SIGMA_RA2I = %f\n", SIGMA_RA2I);
		printf("B_I2RA = %f\n", B_I2RA);
		printf("SIGMA_I2RA = %f\n", SIGMA_I2RA);

		printf("\nPool\n");
		printf("N_RA = %d\n", N_RA);
		printf("N_I = %d\n", N_I);

		printf("\nChain width\n");		
		printf("N_TR = %d\n", N_TR);
		printf("Nss = %d\n", Nss);

		printf("\nNoise:\n");
		printf("mu_soma = %f\n", mu_soma);
		printf("sigma_soma = %f\n", sigma_soma);
		printf("mu_dend = %f\n", mu_dend);
		printf("sigma_dend = %f\n", sigma_dend);

		printf("\nSTDP constants\n");
		printf("BETA = %f\n", BETA);
		printf("BETA_SUPERSYNAPSE = %f\n", BETA_SUPERSYNAPSE);
		printf("T_P = %f\n", T_P);
		printf("TAU_P = %f\n", TAU_P);
		printf("T_D = %f\n", T_D);
		printf("TAU_D = %f\n", TAU_D);
		printf("A_P = %f\n", A_P);
		printf("A_D = %f\n", A_D);
		printf("A_P_SUPER = %f\n", A_P_SUPER);
		printf("A_D_SUPER = %f\n", A_D_SUPER);
		printf("ACTIVATION = %f\n", ACTIVATION);
		printf("SUPERSYNAPSE_THRESHOLD = %f\n", SUPERSYNAPSE_THRESHOLD);
		printf("Gmax = %f\n", G_MAX);

		printf("\nGABA developmental switch\n");
		printf("GABA_DOWN = %f\n", GABA_DOWN);
		printf("BURST_RATE_THRESHOLD = %f\n", BURST_RATE_THRESHOLD);
		printf("RATE_WINDOW = %d\n", RATE_WINDOW);

	}


}

void PoolParallel::read_all_connections(const char* filename, double** weights)
{
    std::ifstream ifile;

    ifile.open(filename, std::ios::binary | std::ios::in);

    int N; // number of RA neurons in the network

    ifile.read(reinterpret_cast<char *>(&N), sizeof(N));
    
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            ifile.read(reinterpret_cast<char *>(&weights[i][j]), sizeof(weights[i][j]));
        }

    }

	ifile.close();
}

void PoolParallel::read_all_connections_from_net(const char* filename, double*** weights)
{

    std::string line;
    std::ifstream ifile; // create ifstream object as input stream

    ifile.open(filename, std::ios::in); // relate ifstream object to input file

    if (!ifile.is_open())
	std::cout << "Error while opening the file" << std::endl;
    
    //while(getline(ifile, line))
    //{
    std::string firstLine;
    
    getline(ifile, firstLine); // read first line
    std::istringstream ss(firstLine); // cretae istringstream object
    std::vector<std::string> items; // vector with "words" in a line
    std::string item; // word in line
	
    while (std::getline(ss, item, ' ')) // split first line into words
    {
         items.push_back(item);
    }

	//for (int i = 0; i < items.size(); i++)
    int N = stoi(items[1], nullptr); // convert word to integer; get number of vertices

    //std::vector<int>* target_id = new vector<int>[N]; // array of vectors with targets
    //std::vector<double>* target_G = new vector<double>[N]; // array of vectors with conductances of synapses

    items.clear();

    //std::cout << "Number of vertices = " << N << std::endl;

    int count = 0;

    while(std::getline(ifile, line)) // read arcs
    {
        count++;

	if (count > N + 1)
        {
	    std::istringstream ss(line); // put line into content of istringstream object
	    //cout << "line = " << line << endl;
	    //cout << "ss.str = " << ss.str() << endl;
	    while (std::getline(ss, item, ' ')) // split line into words
	    {
	        //cout << "item = " << item << endl;
		items.push_back(item);

	    }
	    //cout << "items[0] = " << items[0] << endl;
	    //cout << "items[1] = " << items[1] << endl;
	    //cout << "items[2] = " << items[2] << endl;

	    int source = std::stoi(items[0], nullptr) - 1;
	    int target = std::stoi(items[1], nullptr) - 1;
	    double G = std::stod(items[2]);

	    (*weights)[source][target] = G;

	    items.clear();

	}

    }

    ifile.close();

}

void PoolParallel::read_connections_from_net(const char* filename, std::vector<int>** target_id, std::vector<double>** target_G)
{

    std::string line;
    std::ifstream ifile; // create ifstream object as input stream

    ifile.open(filename, std::ios::in); // relate ifstream object to input file

    if (!ifile.is_open())
	std::cout << "Error while opening the file" << std::endl;
    
    //while(getline(ifile, line))
    //{
    std::string firstLine;
    
    getline(ifile, firstLine); // read first line
    std::istringstream ss(firstLine); // cretae istringstream object
    std::vector<std::string> items; // vector with "words" in a line
    std::string item; // word in line
	
    while (std::getline(ss, item, ' ')) // split first line into words
    {
         items.push_back(item);
    }

	//for (int i = 0; i < items.size(); i++)
    int N = stoi(items[1], nullptr); // convert word to integer; get number of vertices

    //std::vector<int>* target_id = new vector<int>[N]; // array of vectors with targets
    //std::vector<double>* target_G = new vector<double>[N]; // array of vectors with conductances of synapses

    items.clear();

    //std::cout << "Number of vertices = " << N << std::endl;

    int count = 0;

    while(std::getline(ifile, line)) // read arcs
    {
        count++;

	if (count > N + 1)
        {
	    std::istringstream ss(line); // put line into content of istringstream object
	    //cout << "line = " << line << endl;
	    //cout << "ss.str = " << ss.str() << endl;
	    while (std::getline(ss, item, ' ')) // split line into words
	    {
	        //cout << "item = " << item << endl;
		items.push_back(item);

	    }
	    //cout << "items[0] = " << items[0] << endl;
	    //cout << "items[1] = " << items[1] << endl;
	    //cout << "items[2] = " << items[2] << endl;

	    int source = std::stoi(items[0], nullptr) - 1;
	    int target = std::stoi(items[1], nullptr) - 1;
	    double G = std::stod(items[2]);

	    (*target_id)[source].push_back(target);
	    (*target_G)[source].push_back(G);

	    items.clear();

	}

    }

    ifile.close();

}

void PoolParallel::send_simulation_parameters()
{
	MPI_Status status;
	int* syn_num_RA_active_global;
	int* syn_num_RA_active_local = new int[N_RA_local];
	
	int* sendcounts_syn_num_RA;
	int* displs_syn_num_RA;
	int* syn_num_RA_super_global;
	int* syn_num_RA_super_local = new int[N_RA_local];
	
	if (MPI_rank == 0)
	{
		syn_num_RA_active_global = new int[N_RA];

		sendcounts_syn_num_RA = new int[MPI_size];
		syn_num_RA_super_global = new int[N_RA];
		displs_syn_num_RA = new int[MPI_size];
		
		// prepare number of active and supersynaptic connections for each neuron

		sendcounts_syn_num_RA[0] = N_RA_sizes[0];
		displs_syn_num_RA[0] = 0;

		for (int i = 1; i < MPI_size; i++)
		{
			sendcounts_syn_num_RA[i] = N_RA_sizes[i];
			displs_syn_num_RA[i] = displs_syn_num_RA[i-1] + sendcounts_syn_num_RA[i-1];
		}

		for (int i = 0; i < N_RA; i++)
		{
			syn_num_RA_active_global[i] = (int) active_synapses_global[i].size();
			//printf("Master process; syn_num_RA[%d] = %d\n", i, syn_num_RA[i]);
			syn_num_RA_super_global[i] = (int) active_supersynapses_global[i].size();
		}



    }

	// send number of connections for each neuron

	MPI_Scatterv(&syn_num_RA_active_global[0], sendcounts_syn_num_RA, displs_syn_num_RA, MPI_INT,
		&syn_num_RA_active_local[0], N_RA_local, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Scatterv(&syn_num_RA_super_global[0], sendcounts_syn_num_RA, displs_syn_num_RA, MPI_INT,
		&syn_num_RA_super_local[0], N_RA_local, MPI_INT, 0, MPI_COMM_WORLD);


	// send gaba potentials
	MPI_Scatterv(&gaba_potential_global[0], sendcounts_syn_num_RA, displs_syn_num_RA, MPI_DOUBLE,
		&gaba_potential_local[0], N_RA_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
	// send maturation indicators
	MPI_Scatterv(&mature_global[0], sendcounts_syn_num_RA, displs_syn_num_RA, MPI_INT,
		&mature_local[0], N_RA_local, MPI_INT, 0, MPI_COMM_WORLD);
	
	// send maturation trigger indicators
	MPI_Scatterv(&maturation_triggered_global[0], sendcounts_syn_num_RA, displs_syn_num_RA, MPI_INT,
		&maturation_triggered_local[0], N_RA_local, MPI_INT, 0, MPI_COMM_WORLD);
	
	// send internal times
    MPI_Bcast(&internal_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&network_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&trial_number, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	//for (int i = 0; i < N_RA_local; i++)
	//	printf("My rank = %d; mature_local[%d] = %d\n", MPI_rank, Id_RA_local[i], mature_local[i]);

	//for (int i = 0; i < N_RA_local; i++)
		//printf("My rank = %d; sun_num_RA[%d] = %d\n", MPI_rank, i, syn_num_RA_local[i]);

	//	resize vectors that store active and supersynapses
	for (int i = 0; i < N_RA_local; i++)
	{
		active_synapses_local[i].resize(syn_num_RA_active_local[i]);
		active_supersynapses_local[i].resize(syn_num_RA_super_local[i]);
	}

	// send connection ID

	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_RA_local; i++)
		{
			active_synapses_local[i] = active_synapses_global[i];
			active_supersynapses_local[i] = active_supersynapses_global[i];
			
			for (int j = 0; j < N_RA; j++)
				weights_local[i][j] = weights_global[i][j];

		}

		int offset_RA = N_RA_local;
		int offset_I = N_I_local;

		for (int i = 1; i < MPI_size; i++)
		{
            //  send ID RA targets
			for (int j = 0; j < sendcounts_syn_num_RA[i]; j++)
			{
				//printf("Master. Here!\n");
				MPI_Send(&active_synapses_global[offset_RA+j][0], syn_num_RA_active_global[offset_RA+j], MPI_INT, i, 0, MPI_COMM_WORLD);
				MPI_Send(&active_supersynapses_global[offset_RA+j][0], syn_num_RA_super_global[offset_RA+j], MPI_INT, i, 0, MPI_COMM_WORLD);
				MPI_Send(&weights_global[offset_RA+j][0], N_RA, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
			offset_RA += sendcounts_syn_num_RA[i];


		}

	}
	else
	{
        // receive ID and weights of RA targets
		for (int i = 0; i < N_RA_local; i++)
		{
			//printf("Rank = %d. Here!\n", MPI_rank);
			MPI_Recv(&active_synapses_local[i][0], syn_num_RA_active_local[i], MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&active_supersynapses_local[i][0], syn_num_RA_super_local[i], MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&weights_local[i][0], N_RA, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        }

	}

	
	//for (int i = 0; i < N_RA_local; i++)
	//{
		//for (int j = 0; j < syn_num_RA_active_local[i]; j++)
		//	printf("My rank = %d; Neuron %d has active connection on %d\n", MPI_rank, Id_RA_local[i], active_synapses_local[i][j]);
		
		//for (int j = 0; j < syn_num_RA_super_local[i]; j++)
		//	printf("My rank = %d; Neuron %d has super connection on %d\n", MPI_rank, Id_RA_local[i], active_supersynapses_local[i][j]);
	
	//	printf("My rank = %d; Neuron %d on %d: %e\n", MPI_rank, Id_RA_local[i], N_RA - 1, weights_local[i][N_RA - 1]);
	//}
	
	if (MPI_rank == 0)
	{
        delete [] sendcounts_syn_num_RA;
        delete [] syn_num_RA_active_local;
        delete [] syn_num_RA_active_global;
        delete [] syn_num_RA_super_local;
        delete [] syn_num_RA_super_global;
        delete [] displs_syn_num_RA;
    }
	else
	{
        delete [] syn_num_RA_active_local;
        delete [] syn_num_RA_super_local;
	}
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
			syn_num_RA[i] = syn_ID_RA_I_global[i].size();
			//printf("Master process; syn_num_RA[%d] = %d\n", i, syn_num_RA[i]);
		}

		for (int i = 0; i < N_I; i++)
		{
			syn_num_I[i] = syn_ID_I_RA_global[i].size();

		}


    }

	// send number of connections for each neuron

	MPI_Scatterv(&syn_num_RA[0], sendcounts_syn_num_RA, displs_syn_num_RA, MPI_INT,
		&syn_num_RA_local[0], N_RA_local, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Scatterv(&syn_num_I[0], sendcounts_syn_num_I, displs_syn_num_I, MPI_INT,
		&syn_num_I_local[0], N_I_local, MPI_INT, 0, MPI_COMM_WORLD);

	//for (int i = 0; i < N_RA_local; i++)
		//printf("My rank = %d; sun_num_RA[%d] = %d\n", MPI_rank, i, syn_num_RA_local[i]);

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

void PoolParallel::update_synaptic_info()
{
	for (int i = 0; i < N_RA_local; i++)
	{
		if (((int) active_supersynapses_local[i].size() == Nss)&&((int) active_synapses_local[i].size() == Nss))
			remodeled_local[i] = true;

		for (int j = 0; j < (int) active_synapses_local[i].size(); j++)
		{
			int target = active_synapses_local[i][j];
			active_local[i][target] = true;
		}

		for (int j = 0; j < (int) active_supersynapses_local[i].size(); j++)
		{
			int target = active_supersynapses_local[i][j];
			supersynapses_local[i][target] = true;
		}
		//printf("Neuron %d; remodeled_local[%d] = %d\n", Id_RA_local[i], Id_RA_local[i], remodeled_local[i]);
	}

}

void PoolParallel::initialize_connections(double Gei_mean, double Gei_var, double Gie_mean, double Gie_var)
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
		        	 double G = Gei_mean + generator.normal_distribution() * Gei_var;

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
		        	 double G = Gie_mean + generator.normal_distribution() * Gie_var;

		             weights_I_RA_global[i].push_back(G);
		             syn_ID_I_RA_global[i].push_back(j);
		         }
			 }
		 }
	}

}

void PoolParallel::print_invariable_connections()
{
    if (MPI_rank == 0)
    {
        for (int i = 0; i < N_RA; i++)
            printf("RA neuron %d has %d connections to I neurons\n", i, weights_RA_I_global[i].size());

        for (int i = 0; i < N_I; i++)
            printf("I neuron %d has %d connections to RA neurons\n", i, weights_I_RA_global[i].size());
    }
}

void PoolParallel::reset_after_trial()
{

    for (int i = 0; i < N_RA; i++)
    {
        spikes_in_trial_soma_global[i].clear();
        spikes_in_trial_dend_global[i].clear();
    }

	for (int i = 0; i < N_I; i++)	
		spikes_in_trial_interneuron_global[i].clear();

    for (int i = 0; i < N_RA_local; i++)
    {
        HVCRA_local[i].reset();
        spikes_in_trial_soma_local[i].clear();
        spikes_in_trial_dend_local[i].clear();

    }
    for (int i = 0; i < N_I_local; i++)
	{
        HVCI_local[i].reset();
		spikes_in_trial_interneuron_local[i].clear();
	}
}

void PoolParallel::randomize_after_trial()
{
    for (int i = 0; i < N_RA; i++)
    {
        spikes_in_trial_soma_global[i].clear();
        spikes_in_trial_dend_global[i].clear();

        spike_times_dend_global[i] = -100;
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
        current_injection_time = DEMOTIVATION_WINDOW + generator.random(trial_duration-2*DEMOTIVATION_WINDOW);

        std::function<double (double)> f = std::bind(&training_current, current_injection_time, _1);

        for (int i = 0; i < N_TR; i++)
            HVCRA_local[i].set_dend_current(f);
    }
}

void PoolParallel::set_simulation_parameters(double interval, double tS, double mu_s, double sigma_s, double mu_d, double sigma_d)
{
	
	mu_soma = mu_s;
	sigma_soma = sigma_s;
	mu_dend = mu_d;
	sigma_dend = sigma_d;

	size = (int) round(interval/tS) + 1;
	timeStep = tS;
	trial_duration = interval;

    //printf("size = %d\n", size);

	for (int i = 0; i < N_RA_local; i++)
	{	
		HVCRA_local[i].set_noise_generator(&generator);
		HVCRA_local[i].set_white_noise(mu_soma, sigma_soma, mu_dend, sigma_dend);
		HVCRA_local[i].set_dynamics(interval, tS);
	}
	for (int i = 0; i < N_I_local; i++)
	{
		HVCRA_local[i].set_noise_generator(&generator);
		HVCI_local[i].set_dynamics(interval, tS);
	}
}

int PoolParallel::get_trial_number()
{
    return trial_number;
}


void PoolParallel::mature_chain_test(const int num_trials, const char* file_soma_spikes, const char* file_dend_spikes, const char* file_chain_test)
{

	std::vector<std::vector<double>> average_dendritic_spike_time; // array with average dendritic spike time in every trial

	average_dendritic_spike_time.resize(N_RA);

	for (int i = 0; i < num_trials; i++)
	{
		this->mature_trial();
		this->gather_mature_data(average_dendritic_spike_time);
        //this->write_RA("/home/eugene/Output/RA0.bin", 0);
		this->write_soma_time_info(file_soma_spikes);
		this->write_dend_time_info(file_dend_spikes);
		this->randomize_after_trial();
		
	}
	
	// process dendritic spikes

	std::vector<double> mean_burst_time; // average of dendritic spike time
	std::vector<double> std_burst_time; // standard deviation of dendritic spike time
    std::vector<int> num_dend_spikes; // number of dendritic spikes in all trials


	mean_burst_time.resize(N_RA);
	std_burst_time.resize(N_RA);
    num_dend_spikes.resize(N_RA);

	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_RA; i++)
		{
            num_dend_spikes[i] = static_cast<int>(average_dendritic_spike_time[i].size());
			if (num_dend_spikes[i] > 0)
            {

                //for (int j = 0; j < (int) average_dendritic_spike_time[i].size(); j++)
                  //  printf("average_dendritic_spike_time[%d][%d] = %f\n", i, j, average_dendritic_spike_time[i][j]);
				mean_burst_time[i] = std::accumulate(average_dendritic_spike_time[i].begin(), average_dendritic_spike_time[i].end(), 0.0) / static_cast<double>(num_dend_spikes[i]);

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

	this->write_chain_test(num_trials, num_dend_spikes, mean_burst_time, std_burst_time, file_chain_test);
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
	double current_injection_time = 10; //ms
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
            
            // RK4 step
            HVCRA_local[i].R4_step_no_target_update();
            
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
            } 

            if (HVCRA_local[i].get_fired_dend())
            {
                spikes_in_trial_dend_local[i].push_back(internal_time);
            }
		}

		for (int i = 0; i < N_I_local; i++)
		{
            HVCI_local[i].R4_step_no_target_update();
            
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
                    HVCRA_local[i].raise_AMPA(update_Ge_RA_global[Id_RA_local[i]]); // update conductance
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

void PoolParallel::trial(int training)
{
	int some_RA_neuron_fired_soma_local;
	int some_RA_neuron_fired_dend_local;

	int some_I_neuron_fired_local;
	int some_RA_neuron_fired_soma_global;
	int some_RA_neuron_fired_dend_global;

	int some_I_neuron_fired_global;

    std::vector<int> RA_neurons_fired_dend_global;
    std::vector<int> RA_neurons_fired_dend_realID;

    trial_number++;
    
    internal_time = (trial_number - 1) * trial_duration;
    network_time = internal_time + network_update_frequency;

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
            
            // RK4 step
            HVCRA_local[i].R4_step_no_target_update();
            
            // if some neuron produced somatic spike, do LTD for all previous dendritic spikes
            if (HVCRA_local[i].get_fired_soma())
            {
                some_RA_neuron_fired_soma_local = 1;
                spikes_in_trial_soma_local[i].push_back(internal_time);
                last_soma_spikes_local[i].push_back(internal_time);

                while ((!last_soma_spikes_local[i].empty())&&
                    (last_soma_spikes_local[i].back() - last_soma_spikes_local[i].front() > STDP_WINDOW))
                    last_soma_spikes_local[i].pop_front();
            
                
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

                // if neuron is saturated apply LTD only to supersynapses and supply glutamate only to supersynapses
                if (static_cast<int>(active_supersynapses_local[i].size()) == Nss)
                {
                    for (size_t k = 0; k < active_supersynapses_local[i].size(); k++)
                    {
                        int supersynapse_id = active_supersynapses_local[i][k];

						//update_g_local[target_ID] += g_KICK; // update glutamate of supersynaptic target

                        double dt = internal_time - spike_times_dend_global[supersynapse_id];

                        if ( (Id_RA_local[i] != supersynapse_id) && (dt < STDP_WINDOW) )
                        {
                        
         //                   printf("LTD from neuron %d onto %d; somatic spike: %f; dendritic spike: %f\n", Id_RA_local[fired_ID], target_ID,
			    	//spike_times_soma_local[fired_ID], spike_times_dend_global[target_ID]);
                            LTD(weights_local[i][supersynapse_id], dt);
                            update_synapse(i, supersynapse_id);
                       		
					   	}
                    }

					// if some supersynapse desaturated, update all synapses
                	if (static_cast<int>(active_supersynapses_local[i].size()) < Nss)
						for (int j = 0; j < N_RA; j++)
							this->update_synapse(i, j);
                }
                // if not saturated apply LTD rule with the last dendritic spike of all neurons and add glutamate to all neuron except the fired one
                else
                {
                    for (int j = 0; j < N_RA; j++)
                    {
                        double dt = internal_time - spike_times_dend_global[j];

                        if ( (Id_RA_local[i] != j) && (dt < STDP_WINDOW) )
                        {
                           // printf("LTD from neuron %d onto %d; somatic spike: %f; dendritic spike: %f\n", Id_RA_local[fired_ID], j,
		//	    	spike_times_soma_local[fired_ID], spike_times_dend_global[j]);
                           	LTD(weights_local[i][j], dt);
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

                //printf("My rank = %d; Dend neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], spike_times_dend_local[i]);

            }

		} // end for i -> N_RA_local

		for (int i = 0; i < N_I_local; i++)
		{
            HVCI_local[i].R4_step_no_target_update();
            
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
                    HVCRA_local[i].raise_AMPA(update_Ge_RA_global[Id_RA_local[i]]); // update conductance
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
                MPI_Status status;

                // get total number of fired neurons
                MPI_Allreduce(&num_RA_fired_dend_local, &num_RA_fired_dend_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                RA_neurons_fired_dend_global.resize(num_RA_fired_dend_global);


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
                for (int i = 0; i < RA_neurons_fired_dend_global.size(); i++)
                {
                    spike_times_dend_global[RA_neurons_fired_dend_global[i]] = internal_time;

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
                    if (static_cast<int>(active_supersynapses_local[i].size()) == Nss)
                    {
                        for (size_t j = 0; j < RA_neurons_fired_dend_global.size(); j++)
                        {
                            int fired_ID = RA_neurons_fired_dend_global[j];

                            std::vector<int>::iterator pos = std::find(active_supersynapses_local[i].begin(),
                                        active_supersynapses_local[i].end(), fired_ID);

                            if (pos!=active_supersynapses_local[i].end())
                            {
                                for  (size_t k = 0; k < last_soma_spikes_local[i].size(); k++)
                                {
                                    double dt = internal_time - last_soma_spikes_local[i][k];

                                    if (dt < STDP_WINDOW)
                                    {
                                        LTP(weights_local[i][fired_ID], dt);
         //                             printf("LTP from neuron %d onto %d; somatic spike at %f; dendritic spike at %f\n", Id_RA_local[i], fired_ID,
             //                                               last_soma_spikes_local[i][k], spike_times_dend_global[fired_ID]);
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
                                for  (size_t k = 0; k < last_soma_spikes_local[i].size(); k++)
                                {
                                    double dt = internal_time - last_soma_spikes_local[i][k];

                                    if (dt < STDP_WINDOW)
                                    {
           //                             printf("LTP from neuron %d onto %d; somatic spike at %f; dendritic spike: %f\n", Id_RA_local[i], fired_ID,
                            //                                last_soma_spikes_local[i][k], spike_times_dend_global[fired_ID]);
                                        LTP(weights_local[i][fired_ID], dt);
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
                    if ( (active_supersynapses_local[i].size() == Nss) && (!remodeled_local[i]) )
                    {
					    this->axon_remodeling(i);

				    }
                }
	        	
				// update fired arrays and indicators
				some_RA_neuron_fired_dend_local = 0;
            	some_RA_neuron_fired_dend_global = 0;
            	
				RA_neurons_fired_dend_global.clear();
            	RA_neurons_fired_dend_realID.clear();
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

		//if (std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].end(), 0.0) >= 1)
		//{
		//	std::cout << "Recent bursts of neuron " << Id_RA_local[i] << " :" << std::endl;
			
		//	for (size_t j = 0; j < num_bursts_in_recent_trials[i].size(); j++)
		//		std::cout << num_bursts_in_recent_trials[i][j] << "\t";
			
		//	std::cout << std::endl;
		//}
	}

	

	this->update_Ei();

//printf("internal time = %f\n", internal_time);
	//printf("t*timeStep = %f\n", t*timeStep);
}

void PoolParallel::update_Ei()
{
	for (int i = 0; i < N_RA_local; i++)
	{
	// if not a training neuron and not mature but maturation process was triggered
		if ( (Id_RA_local[i] >= N_TR) && (maturation_triggered_local[i] == 1) && (mature_local[i] == 0) )
		{
			gaba_potential_local[i] = gaba_potential_local[i] - GABA_DOWN; // decrease gaba potential
					
			if (gaba_potential_local[i] <= E_GABA_MATURE) // if gaba potential hit the bottom, make neuron mature
			{
				gaba_potential_local[i] = E_GABA_MATURE;
				mature_local[i] = 1;
			}

		}
		// if not a training neuron and if maturation is not triggered
		else if ( (Id_RA_local[i] >= N_TR) && (maturation_triggered_local[i] == 0) )
		{
			double average_rate = std::accumulate(num_bursts_in_recent_trials[i].begin(), 
												 num_bursts_in_recent_trials[i].end(), 0.0) / static_cast<double>(RATE_WINDOW);
			
			// check trigger maturation if rate exceeds threshold
			if (average_rate >= BURST_RATE_THRESHOLD)
				maturation_triggered_local[i] = 1;

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
        std::vector<int>::iterator pos = std::find(active_supersynapses_local[i].begin(),
                                active_supersynapses_local[i].end(), syn_ID);
        // if not erase it
        if (pos == active_supersynapses_local[i].end())
        {
            active_to_erase.push_back(active_synapses_local[i][j]);
            active_local[i][syn_ID] = false;
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

    remodeled_local[i] = true;
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

    if ( (w >= ACTIVATION) && (!active_local[i][j]) && (!remodeled_local[i]) )
    {
       // printf("Activated synapse from %d onto %d; weight = %f\n", Id_RA_local[i], j, w);
        active_local[i][j] = true;
        active_synapses_local[i].push_back(j);
    }

    if ( (w >= SUPERSYNAPSE_THRESHOLD) && (!supersynapses_local[i][j]) && (static_cast<int>(active_supersynapses_local[i].size()) < Nss) )
    {
       // printf("Activated supersynapse from %d onto %d; weight = %f\n", Id_RA_local[i], j, w);
        supersynapses_local[i][j] = true;
        active_supersynapses_local[i].push_back(j);
    }

    if ( (w < SUPERSYNAPSE_THRESHOLD) && (supersynapses_local[i][j]) )
    {
       // printf("Deactivated supersynapse from %d onto %d; weight = %f\n", Id_RA_local[i], j, w);
        supersynapses_local[i][j] = false;
        remodeled_local[i] = false;
        std::vector<int>::iterator pos = std::find(active_supersynapses_local[i].begin(),
                                                    active_supersynapses_local[i].end(), j);

        if (pos!= active_supersynapses_local[i].end())
            active_supersynapses_local[i].erase(pos);
        else
            std::cout << "Supersynapse " << Id_RA_local[i] << " -> " << j << " to be erased is not found!" << std::endl;
    }

    if ( (w < ACTIVATION) && (active_local[i][j]) )
    {
      //  printf("Deactivated synapse from %d onto %d; w = %f\n", Id_RA_local[i], j, w);
        active_local[i][j] = false;

        std::vector<int>::iterator pos = std::find(active_synapses_local[i].begin(),
                                                    active_synapses_local[i].end(), j);

        if (pos!= active_synapses_local[i].end())
            active_synapses_local[i].erase(pos);
        else
            std::cout << "Active synapse " << Id_RA_local[i] << " -> " << j << " to be erased is not found!" << std::endl;


    }

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

void PoolParallel::LTP(double &w, double t)
{
    //printf("LTP\n");
	if (t <= T_P)
    {
        //if ((w + R * A_P * G_P * ((1 + F_0) * t / T_P - F_0))<0)
          //  printf("LTP. Linear interval. Weight became negative. w = %f\n", w);
        //double temp = R * A_P * G_P * ((1 + F_0) * t / T_P - F_0);
        //if (temp<0)
         //   printf("Negative LTP!!! t = %f; temp = %f\n", t, temp);
        //else
         //   printf("Positive LTP!!! t = %f; temp = %f\n", t, temp);
		//w = w + temp;
		if (w >= SUPERSYNAPSE_THRESHOLD)
		{
			w = w + R * A_P_SUPER * G_P * ((1 + F_0) * t / T_P - F_0);
		}
		else
		{
			w = w + R * A_P * G_P * ((1 + F_0) * t / T_P - F_0);
		}
	

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
   		if (w >= SUPERSYNAPSE_THRESHOLD)
		{

   			w = w + R * A_P_SUPER * G_P * exp(-(t - T_P) / TAU_P);
		}
		else
		{
   			w = w + R * A_P * G_P * exp(-(t - T_P) / TAU_P);
		}

        //std::cout << "w = " << w << std::endl;
    }

    if (w < 0)
        w = 0;

    if (w > G_MAX)
        w = G_MAX;
}

void PoolParallel::LTD(double &w, double t)
{
	if (t <= T_P)
	{

        //if ((w - R * A_D * w * ((1 - F_0) * t / T_D + F_0))<0)
            //printf("LTD. Linear interval. Weight became negative. w = %f\n", w);
		//w = w - R * A_D * w * ((1 - F_0) * t / T_D + F_0);
		//w = w - R * A_D * ((1 - F_0) * t / T_D + F_0);
		if (w >= SUPERSYNAPSE_THRESHOLD)
		{

			w = w - R * (A_P_SUPER * G_P * F_0 + (A_D_SUPER - A_P_SUPER * G_P * F_0) * t / T_D);
		}
		else
		{
			w = w - R * (A_P * G_P * F_0 + (A_D - A_P * G_P * F_0) * t / T_D);
        }
		
		if (w < 0)
            w = 0;

       // std::cout << "w = " << w << std::endl;
	}
	else
	{
        //if ((w - R * A_D * w * exp(-(t - T_D) / TAU_D))<0)
          //  printf("LTD. Exponential interval. Weight became negative. w = %f\n", w);

		//w = w - R * A_D * w * exp(-(t - T_D) / TAU_D);
		if (w >= SUPERSYNAPSE_THRESHOLD)
		{
			w = w - R * A_D_SUPER * exp(-(t - T_D) / TAU_D);
		}
		else
		{
			w = w - R * A_D * exp(-(t - T_D) / TAU_D);
		}

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
            spikes_in_trial_soma_global[i].resize(spike_num_soma_global[i]);
            spikes_in_trial_dend_global[i].resize(spike_num_dend_global[i]);

            //printf("Master, spike_num_global[%d] = %d\n", i, spike_num_global[i]);
        }


        // Copy from master's local arrays
        for (int i = 0; i < N_RA_local; i++)
        {
            spikes_in_trial_soma_global[i] = spikes_in_trial_soma_local[i];
            spikes_in_trial_dend_global[i] = spikes_in_trial_dend_local[i];
        }

			
    // Gather from others
		int N = N_RA_sizes[0]; // number of RA neurons in the processes with lower rank

        for (int i = 1; i < MPI_size; i++)
        {

            for (int j = 0; j < N_RA_sizes[i]; j++)
            {
                int count;
				int receive_index = N + j;

                if (spike_num_soma_global[receive_index] != 0)
                {
                    MPI_Recv(&spikes_in_trial_soma_global[receive_index][0],
                        spike_num_soma_global[receive_index], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv spikes in trial; from i = %d  count = %d\n", i, count);
                }

                if (spike_num_dend_global[receive_index] != 0)
                {
                    MPI_Recv(&spikes_in_trial_dend_global[receive_index][0],
                        spike_num_dend_global[receive_index], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

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
        supersyn_sizes_local[i] = active_supersynapses_local[i].size();
        syn_sizes_local[i] = active_synapses_local[i].size();
        spike_num_soma_local[i] = spikes_in_trial_soma_local[i].size();
        spike_num_dend_local[i] = spikes_in_trial_dend_local[i].size();

        //printf("Rank = %d, supersyn_sizes_local[%d] = %d\n", MPI_rank, Id_RA_local[i], supersyn_sizes_local[i]);
        //printf("Rank = %d, syn_sizes_local[%d] = %d\n", MPI_rank, Id_RA_local[i], syn_sizes_local[i]);
        //printf("Rank = %d, spike_num_local[%d] = %d\n", MPI_rank, Id_RA_local[i], spike_num_local[i]);
    }

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
    
    MPI_Gatherv(&maturation_triggered_local[0], N_RA_local, MPI_INT,
        &maturation_triggered_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);
    
	MPI_Gatherv(&gaba_potential_local[0], N_RA_local, MPI_DOUBLE,
        &gaba_potential_global[0], recvcounts_RA, displs_RA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if (MPI_rank == 0)
    {
        for (int i = 0; i < N_RA; i++)
        {
            //printf("Master; supersyn_sizes_global[%d] = %d\n", i, supersyn_sizes_global[i]);
            active_supersynapses_global[i].resize(supersyn_sizes_global[i]);
            active_synapses_global[i].resize(syn_sizes_global[i]);
            spikes_in_trial_soma_global[i].resize(spike_num_soma_global[i]);
            spikes_in_trial_dend_global[i].resize(spike_num_dend_global[i]);

            //printf("Master, supersyn_sizes_global[%d] = %d\n", i, supersyn_sizes_global[i]);
            //printf("Master, syn_sizes_global[%d] = %d\n", i, syn_sizes_global[i]);
            //printf("Master, spike_num_global[%d] = %d\n", i, spike_num_global[i]);
        }


        // Copy from master's local arrays
        for (int i = 0; i < N_RA_local; i++)
        {
            active_supersynapses_global[i] = active_supersynapses_local[i];
            active_synapses_global[i] = active_synapses_local[i];
            spikes_in_trial_soma_global[i] = spikes_in_trial_soma_local[i];
            spikes_in_trial_dend_global[i] = spikes_in_trial_dend_local[i];

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
               weights_global[i][j] = weights_local[i][j];

               //printf("Master; weights_global[%d][%d] = %f\n", i, j,
                 //   weights_global[i][j]);
            }

        }
    // Gather from others
		int N = N_RA_sizes[0]; // number of RA neurons in the processes with lower rank

        for (int i = 1; i < MPI_size; i++)
        {

            for (int j = 0; j < N_RA_sizes[i]; j++)
            {
                int count;
				int receive_index = N + j;
                MPI_Recv(&weights_global[receive_index][0],
                                        N_RA, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_DOUBLE, &count);
                //printf("Recv weights; from i = %d  count = %d\n", i, count);
                if (supersyn_sizes_global[receive_index] != 0)
                {
                    MPI_Recv(&active_supersynapses_global[receive_index][0],
                        supersyn_sizes_global[receive_index], MPI_INT, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv supersynapses; from i = %d  count = %d\n", i, count);
                }

                if (syn_sizes_global[receive_index] != 0)
                {
                    MPI_Recv(&active_synapses_global[receive_index][0],
                        syn_sizes_global[receive_index], MPI_INT, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv synapses; from i = %d  count = %d\n", i, count);
                }

                if (spike_num_soma_global[receive_index] != 0)
                {
                    MPI_Recv(&spikes_in_trial_soma_global[receive_index][0],
                        spike_num_soma_global[receive_index], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv spikes in trial; from i = %d  count = %d\n", i, count);
                }

                if (spike_num_dend_global[receive_index] != 0)
                {
                    MPI_Recv(&spikes_in_trial_dend_global[receive_index][0],
                        spike_num_dend_global[receive_index], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

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
                MPI_Send(&active_supersynapses_local[i][0],
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

void PoolParallel::write_chain_test(int num_trials, std::vector<int>& num_dend_spikes, std::vector<double>& mean_burst_time, std::vector<double>& std_burst_time, const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;

        int size;

        out.open(filename, std::ios::out | std::ios::binary);

        // write number of neurons to each file

        out.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));
        out.write(reinterpret_cast<char *>(&num_trials), sizeof(num_trials));

        for (int i = 0; i < N_RA; i++)
        {
            out.write(reinterpret_cast<char *>(&num_dend_spikes[i]), sizeof(num_dend_spikes[i]));
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

        for (int i = 1; i <= N_RA; i++)
        {
            // write neuron's ID number
            out_RA_RA.write(reinterpret_cast<char *>(&i), sizeof(i));

            // write number of connections to RA neurons

            size = active_synapses_global[i-1].size();
            out_RA_RA.write(reinterpret_cast<char *>(&size), sizeof(size));

            for (int j = 0; j < size; j++)
            {
                int k = active_synapses_global[i-1][j];
                out_RA_RA.write(reinterpret_cast<char *>(&k), sizeof(k));
                out_RA_RA.write(reinterpret_cast<char *>(&weights_global[i-1][k]), sizeof(weights_global[i-1][k]));

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

        for (int i = 1; i <= N_RA; i++)
        {
            // write neuron's ID number
            out_RA_RA.write(reinterpret_cast<char *>(&i), sizeof(i));

            // write number of connections to RA neurons

            size = active_supersynapses_global[i-1].size();
            out_RA_RA.write(reinterpret_cast<char *>(&size), sizeof(size));

            for (int j = 0; j < size; j++)
            {
                int k = active_supersynapses_global[i-1][j];
                out_RA_RA.write(reinterpret_cast<char *>(&k), sizeof(k));
                out_RA_RA.write(reinterpret_cast<char *>(&weights_global[i-1][k]), sizeof(weights_global[i-1][k]));

            }

        }
        out_RA_RA.close();
    }
}

void PoolParallel::write_weights(const char * filename)
{
	if (MPI_rank == 0)
	{
        std::ofstream output;

        output.open(filename, std::ios::out | std::ios::binary);

        output.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));

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


void PoolParallel::write_time_info(const char* timeInfoFile)
{

	if (MPI_rank == 0)
	{
		std::ofstream out_timeInfo;

		out_timeInfo.open(timeInfoFile, std::ios::binary | std::ios::out);

		// write simulation information
		out_timeInfo.write(reinterpret_cast<char*>(&trial_number), sizeof(trial_number)); // trial number
		out_timeInfo.write(reinterpret_cast<char*>(&internal_time), sizeof(internal_time)); // internal time of each neuron

		out_timeInfo.write(reinterpret_cast<char*>(&network_time), sizeof(network_time)); // time of the network

		out_timeInfo.close();
	}



}

void PoolParallel::write_sim_info(const char* simInfoFile, int synapses_trials_update, int weights_trials_update)
{

	if (MPI_rank == 0)
	{
		std::ofstream out_simInfo;

		out_simInfo.open(simInfoFile, std::ios::binary | std::ios::out);

		// write simulation information
		out_simInfo.write(reinterpret_cast<char*>(&trial_duration), sizeof(trial_duration)); // trial duration in ms
		out_simInfo.write(reinterpret_cast<char*>(&synapses_trials_update), sizeof(synapses_trials_update)); // how often synaptic information is updated

		out_simInfo.write(reinterpret_cast<char*>(&weights_trials_update), sizeof(weights_trials_update)); // how often weights are updated

		out_simInfo.close();
	}



}

void PoolParallel::write_num_synapses(const char* fileSynapses)
{
	if (MPI_rank == 0)
	{
		std::ofstream out_synapses;

		out_synapses.open(fileSynapses, std::ios::binary | std::ios::out | std::ios::app);


		// count active synapses and supersynapses
		int active = 0;
		int super = 0;

		for (int i = 0; i < N_RA; i++)
		{
			active += (int) active_synapses_global[i].size(); 
			super += (int) active_supersynapses_global[i].size(); 

		}
		
        out_synapses.write(reinterpret_cast<char*>(&trial_number), sizeof(trial_number)); // write current trial numbe
		out_synapses.write(reinterpret_cast<char*>(&active), sizeof(active)); // write number of active synapses
		out_synapses.write(reinterpret_cast<char*>(&super), sizeof(super)); // write number of supersynapses
	
        out_synapses.close();
	}

}

void PoolParallel::write_coordinates(const char* xy_RA, const char* xy_I)
{

    if (MPI_rank == 0)
    {
        std::ofstream out_RA;
        std::ofstream out_I;

        // open files
        out_RA.open(xy_RA, std::ios::binary | std::ios::out);
        out_I.open(xy_I, std::ios::binary | std::ios::out);

        // write number of neurons
        out_RA.write(reinterpret_cast<char*>(&N_RA), sizeof(N_RA));
        out_I.write(reinterpret_cast<char*>(&N_I), sizeof(N_I));

        // write coordinates
        for (int i = 0; i < N_RA; i++)
        {
            out_RA.write(reinterpret_cast<char*>(&xx_RA[i]), sizeof(xx_RA[i]));
            out_RA.write(reinterpret_cast<char*>(&yy_RA[i]), sizeof(yy_RA[i]));
        }

        for (int i = 0; i < N_I; i++)
        {
            out_I.write(reinterpret_cast<char*>(&xx_I[i]), sizeof(xx_I[i]));
            out_I.write(reinterpret_cast<char*>(&yy_I[i]), sizeof(yy_I[i]));
        }
        // close files
        out_RA.close();
        out_I.close();
	}
}

void PoolParallel::write_soma_time_info(const char* filename)
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
                //out.write(reinterpret_cast<char *>(&spikes_in_trial_soma_global[i][j]), sizeof(double));
	        double relative_spike_time = spikes_in_trial_soma_global[i][j] - (trial_number - 1) * trial_duration;
       		//printf("relative_spike_time_soma = %f\n", relative_spike_time);
            out.write(reinterpret_cast<char *>(&relative_spike_time), sizeof(double));
	        }
        }
        //out.write(reinterpret_cast<char *>(spike_times), N_RA*sizeof(double));

        out.close();
    }
}

void PoolParallel::write_dend_time_info(const char* filename)
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
	        double relative_spike_time = spikes_in_trial_dend_global[i][j] - (trial_number - 1) * trial_duration;
        	out.write(reinterpret_cast<char *>(&relative_spike_time), sizeof(double));
			//printf("Neuron %d; relative spike time = %f\n", i, relative_spike_time);
	    }
	}
        //out.write(reinterpret_cast<char *>(spike_times), N_RA*sizeof(double));

        out.close();
    }
}

void PoolParallel::write_interneuron_time_info(const char* filename)
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
	        double relative_spike_time = spikes_in_trial_interneuron_global[i][j] - (trial_number - 1) * trial_duration;
        	out.write(reinterpret_cast<char *>(&relative_spike_time), sizeof(double));
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
        for (int i = 1; i <= N_RA; i++)
        {
            size = syn_ID_RA_I_global[i-1].size();

            out_RA_I.write(reinterpret_cast<char *>(&i), sizeof(i));
            out_RA_I.write(reinterpret_cast<char *>(&size), sizeof(size)); // write neuron's ID

            for (int j = 0; j < size; j++)
            {
                int k = syn_ID_RA_I_global[i-1][j];
                double G = weights_RA_I_global[i-1][j];

                out_RA_I.write(reinterpret_cast<char *>(&k), sizeof(k));
                out_RA_I.write(reinterpret_cast<char *>(&G), sizeof(G));

            }

        }

        // write connections from I to RA
        for (int i = 1; i <= N_I; i++)
        {
            out_I_RA.write(reinterpret_cast<char *>(&i), sizeof(i)); // write neuron's ID number

            size = syn_ID_I_RA_global[i-1].size();
            out_I_RA.write(reinterpret_cast<char *>(&size), sizeof(size)); // write number of targets a neuron has
            for (int j = 0; j < size; j++)
            {
                    int k = syn_ID_I_RA_global[i-1][j];
                    double G = weights_I_RA_global[i-1][j];
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

void PoolParallel::write_pajek_super(const char * filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;
        out.open(filename, std::ios::out);

        out << "*Vertices " << N_RA << "\n";

        for (int i = 0; i < N_TR; i++)
        {

            out << i + 1 << " \"" << i << "\" " << xx_RA[i]/SIDE << " " << yy_RA[i]/SIDE << " ic Green\n";
        }

        for (int i = N_TR; i < N_RA; i++)
        {

            out << i + 1 << " \"" << i << "\" " << xx_RA[i]/SIDE << " " << yy_RA[i]/SIDE << " ic Yellow\n";
        }
        out << "*Arcs\n";
        for (int i = 0; i < N_RA; i++)
        {
            for (int j = 0; j < active_supersynapses_global[i].size(); j++)
                {
                    int syn_ID = active_supersynapses_global[i][j];
                    out << i + 1 << " " << syn_ID + 1 << " " << weights_global[i][syn_ID] << "\n";
                }
        }
    }
}

void PoolParallel::write_pajek_active(const char * filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;
        out.open(filename, std::ios::out);

        out << "*Vertices " << N_RA << "\n";

        for (int i = 0; i < N_TR; i++)
        {

            out << i + 1 << " \"" << i << "\" " << xx_RA[i]/SIDE << " " << yy_RA[i]/SIDE << " ic Green\n";
        }

        for (int i = N_TR; i < N_RA; i++)
        {

            out << i + 1 << " \"" << i << "\" " << xx_RA[i]/SIDE << " " << yy_RA[i]/SIDE << " ic Yellow\n";
        }
        out << "*Arcs\n";
        for (int i = 0; i < N_RA; i++)
        {
            for (int j = 0; j < active_synapses_global[i].size(); j++)
                {
                    int syn_ID = active_synapses_global[i][j];
                    out << i + 1 << " " << syn_ID + 1 << " " << weights_global[i][syn_ID] << "\n";
                }
        }
    }
}

void PoolParallel::write_pajek_all(const char * filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;
        out.open(filename, std::ios::out);

        out << "*Vertices " << N_RA << "\n";

        for (int i = 0; i < N_TR; i++)
        {

            out << i + 1 << " \"" << i << "\" " << xx_RA[i]/SIDE << " " << yy_RA[i]/SIDE << " ic Green\n";
        }

        for (int i = N_TR; i < N_RA; i++)
        {

            out << i + 1 << " \"" << i << "\" " << xx_RA[i]/SIDE << " " << yy_RA[i]/SIDE << " ic Yellow\n";
        }

        out << "*Arcs\n";
        for (int i = 0; i < N_RA; i++)
        {
            for (int j = 0; j < N_RA; j++)
                {
                    if (j != i)
                    {
                        out << i + 1 << " " << j + 1 << " " << weights_global[i][j] << "\n";
                    }
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

        for (int i = 0; i < N_RA; i++)
        {
            out.write(reinterpret_cast<char *>(&gaba_potential_global[i]), sizeof(gaba_potential_global[i]));
            out.write(reinterpret_cast<char *>(&maturation_triggered_global[i]), sizeof(maturation_triggered_global[i]));
            out.write(reinterpret_cast<char *>(&mature_global[i]), sizeof(mature_global[i]));
	    
	    }
		out.close();
	}

}
