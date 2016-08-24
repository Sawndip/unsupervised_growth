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

PoolParallel::PoolParallel(double beta, double beta_s, double Tp, double Td, double tauP, double tauD, double Ap, double Ad, double Ap_super, double Ad_super, double activation, double super_threshold, 
                        double Gmax, int N_ra, int Nic, int NiInC, int N_ss, int N_tr) : BETA(beta), BETA_SUPERSYNAPSE(beta_s), 
                        A_P(Ap), A_D(Ad), T_P(Tp), T_D(Td), TAU_P(tauP), TAU_D(tauD), A_P_SUPER(Ap_super), A_D_SUPER(Ad_super), ACTIVATION(activation), SUPERSYNAPSE_THRESHOLD(super_threshold), G_MAX(Gmax),
				        N_RA(N_ra), num_inh_clusters(Nic), num_inh_in_cluster(NiInC), Nss(N_ss), N_TR(N_tr)
{

    MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);



    // white noise
    double mu_soma = 30; // dc component of white noise to somatic compartment
    double sigma_soma = 100; // variance of white noise to somatic compartment
    double mu_dend = 50; // dc component of white noise to dendritic compartment
    double sigma_dend = 185; // variance of white noise to dendritic compartment

	double delay = 0;
	int num_related_spikes = 0;

	N_I = num_inh_clusters*num_inh_clusters*num_inh_in_cluster;

    N_RA_local = N_RA/MPI_size;
    N_I_local = N_I/MPI_size;


	// if root assign all excessive neurons
	if (MPI_rank == 0)
	{
        if (N_TR > N_RA)
            printf("Number of training neurons exceeds total number of RA neurons!\n");
        else
        {
            if (N_TR > N_RA_local)
                printf("Number of training neurons exceeds number of neurons assigned to master process");
        }

		N_RA_local += N_RA % MPI_size;
		N_I_local += N_I % MPI_size;

	}

	printf("My rank is %d; N_RA_local = %d; N_I_local = %d\n", MPI_rank, N_RA_local, N_I_local);


	HVCRA_local = new HH2_final_pool[N_RA_local];
	HVCI_local = new HHI_final_pool[N_I_local];

    Id_RA_local = new unsigned[N_RA_local];
    Id_I_local = new unsigned[N_I_local];

    update_Ge_RA_local = new double[N_RA];
    update_Gi_RA_local = new double[N_RA];
	update_g_local = new double[N_RA];
    update_Ge_I_local = new double[N_I];

    update_Ge_RA_global = new double[N_RA];
    update_Gi_RA_global = new double[N_RA];
	update_g_global = new double[N_RA];
    update_Ge_I_global = new double[N_I];

	//active_local = new bool*[N_RA];
	//active_supersynapses_local = new bool*[N_RA_local];
	weights_RA_I_local = new std::vector<double>[N_RA_local];
	weights_I_RA_local = new std::vector<double>[N_I_local];
	syn_ID_RA_I_local = new std::vector<unsigned>[N_RA_local];
	syn_ID_I_RA_local = new std::vector<unsigned>[N_I_local];

    weights_global = new double*[N_RA];
    weights_RA_I_global = new std::vector<double>[N_RA];
    weights_I_RA_global = new std::vector<double>[N_I];
    syn_ID_RA_I_global = new std::vector<unsigned>[N_RA];
    syn_ID_I_RA_global = new std::vector<unsigned>[N_I];

    active_global = new bool*[N_RA];
    supersynapses_global = new bool*[N_RA];
	active_supersynapses_global = new std::vector<unsigned>[N_RA];
	active_synapses_global = new std::vector<unsigned>[N_RA];

    weights_local = new double*[N_RA_local];
    active_local = new bool*[N_RA_local];
    supersynapses_local = new bool*[N_RA_local];
    active_supersynapses_local = new std::vector<unsigned>[N_RA_local];
	active_synapses_local = new std::vector<unsigned>[N_RA_local];

	last_soma_spikes_local = new std::deque<double>[N_RA];
	last_soma_spikes_global = new std::deque<double>[N_RA];

    spike_times_soma_global = new double[N_RA];
    spikes_in_trial_soma_global = new std::vector<double>[N_RA];
    spike_times_dend_global = new double[N_RA];
    spikes_in_trial_dend_global = new std::vector<double>[N_RA];

    spike_times_soma_local = new double[N_RA_local];
    spikes_in_trial_soma_local = new std::vector<double>[N_RA_local];
    spike_times_dend_local = new double[N_RA_local];
    spikes_in_trial_dend_local = new std::vector<double>[N_RA_local];

	remodeled_local = new bool[N_RA_local];
	mature_local = new int[N_RA_local];

	mature_global = new int[N_RA];

	for (int i = 0; i < N_RA_local; i++)
	{
        weights_local[i] = new double[N_RA];
        active_local[i] = new bool[N_RA];
        supersynapses_local[i] = new bool[N_RA];
        spike_times_soma_local[i] = -100.0;
        spike_times_dend_local[i] = -100.0;
        remodeled_local[i] = false;
		mature_local[i] = 0;

        for (int j = 0; j < N_RA; j++)
        {
            weights_local[i][j] = 0.0;
            active_local[i][j] = false;
            supersynapses_local[i][j] = false;
        }

		if (MPI_rank == 0)
			Id_RA_local[i] = i;
		else
			Id_RA_local[i] = MPI_rank*N_RA_local + N_RA % MPI_size + i;

	}

    for (int i = 0; i < N_I_local; i++)
	{
		if (MPI_rank == 0)
			Id_I_local[i] = i;
		else
			Id_I_local[i] = MPI_rank*N_I_local + N_I % MPI_size + i;

	}


	for (int i = 0; i < N_RA; i++)
	{
		active_global[i] = new bool[N_RA];
		supersynapses_global[i] = new bool[N_RA];
		weights_global[i] = new double[N_RA];

	}

	for (int i = 0; i < N_RA; i++)
	{
		for (int j = 0; j < N_RA; j++)
		{
			active_global[i][j] = false;
			supersynapses_global[i][j] = false;
			weights_global[i][j] = 0.0;

		}
		spike_times_soma_global[i] = -100.0;
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

	delete[] spike_times_soma_global;
	delete[] spike_times_dend_global;
	delete[] spike_times_soma_local;
	delete[] spike_times_dend_local;

	delete[] active_global;
	delete[] supersynapses_global;
	delete[] weights_global;
	delete[] remodeled_local;
	delete[] mature_local;
	delete[] mature_global;

	delete[] weights_RA_I_local;
	delete[] weights_I_RA_local;
	delete[] syn_ID_RA_I_local;
	delete[] syn_ID_I_RA_local;

	delete[] weights_RA_I_global;
	delete[] weights_I_RA_global;
	delete[] syn_ID_RA_I_global;
	delete[] syn_ID_I_RA_global;

	delete[] update_Ge_RA_local;
	delete[] update_Gi_RA_local;
	delete[] update_g_local;
	delete[] update_Ge_I_local;

	delete[] update_Ge_RA_global;
	delete[] update_Gi_RA_global;
	delete[] update_g_global;
	delete[] update_Ge_I_global;
//	MPI_Finalize();
}


const double PoolParallel::DEMOTIVATION_WINDOW = 200; // demotivation window in ms
// coordinates
const double PoolParallel::CLUSTER_SIZE = 5; // cluster size
const double PoolParallel::MIN_INTERNEURON_DISTANCE = 1; // minimum distance between neurons
const double PoolParallel::LAMBDA_RA2I = 5; // spatial scale of probability of connections decay
const double PoolParallel::A_RA2I = 4.0;
const double PoolParallel::B_RA2I = 0.030;
const double PoolParallel::MEAN_RA2I = 20.0;
const double PoolParallel::SIGMA_RA2I = 30.0;


const double PoolParallel::LAMBDA_I2RA = 4; // spatial scale of probability of connections decay
const double PoolParallel::CONNECT_CONST_I2RA = 4.0;

const double PoolParallel::SIDE = 100; // length of HVC side



// developmental GABA switch
const double PoolParallel::T_GABA = 10000;
const double PoolParallel::E_GABA_IMMATURE = -50;
const double PoolParallel::E_GABA_MATURE = -80;
const int PoolParallel::N_MATURATION = 100;

// constants for STDP-rules
const int PoolParallel::NUM_SOMA_SPIKES = 5; // number of last somatic spikes to store (necessary for LTP rule)
const double PoolParallel::G_P = 0.1;

const double PoolParallel::R = 1;
const double PoolParallel::F_0 = 1.0;

const double PoolParallel::LTP_WINDOW = 50;
const double PoolParallel::LTD_WINDOW = 50;

const double PoolParallel::DELAY_WINDOW = 35;

// NMDA
const double PoolParallel::g_KICK = 0.22; // glutamate kick for NMDA receptors

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

                    for (unsigned j = 0; j < xx_I.size(); j++)
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

void PoolParallel::initialize_inhibitory_clusters()
{
	if (MPI_rank == 0)
	{	
        double xx; // temporary x-coordinate
        double yy; // temporary y-coordinate

        bool close; // are neurons too close or not

        // set coordinates for inhibitory cluster centers HVC(I) neurons on a square lattice
        for (int i = 0; i < num_inh_clusters; i++)
        {
            for (int k = 0; k < num_inh_clusters; k++)
            {
				xx = (double) (i+1) * SIDE / (num_inh_clusters+1);
				yy = (double) (k+1) * SIDE / (num_inh_clusters+1);

				xx_I_cluster.push_back(xx);
				yy_I_cluster.push_back(yy);

            }
        }

		// set coordinates for interneurons inside each cluster

		for (int i = 0; i < num_inh_clusters*num_inh_clusters; i++)
		{
			for (int j = 0; j < num_inh_in_cluster; j++)
			{
                do
                {
					close = false;
					xx = xx_I_cluster[i] + 2 * (generator.random(2*MIN_INTERNEURON_DISTANCE) - MIN_INTERNEURON_DISTANCE);
					yy = yy_I_cluster[i] + 2 * (generator.random(2*MIN_INTERNEURON_DISTANCE) - MIN_INTERNEURON_DISTANCE);

					// check distances to all previous I neurons

                    for (unsigned k = 0; k < xx_I.size(); k++)
                    {
                        if (distance(xx,yy,xx_I[k],yy_I[k]) < MIN_INTERNEURON_DISTANCE)
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
	}
}


void PoolParallel::initialize_RA_for_inh_clusters()
{
	if (MPI_rank == 0)
	{
     	// set coordinates for HVC(RA) neurons
		// for each inhibitory neuron

		bool close; // indicator whether neurons are too close to each other ro ot
		double xx, yy; // temporary coordinates of neurons

		int total_num_inh_clusters = num_inh_clusters * num_inh_clusters;

		int average = N_RA / (total_num_inh_clusters); // average number of RA neurons per single inhibitory cluster
		int N_remaining = N_RA - average * total_num_inh_clusters; // number of RA neurons that need to be distributed in addition to average amount
		
		N_cluster = new int[total_num_inh_clusters]; // array with number of RA neurons per cluster

		for (int i = 0; i < total_num_inh_clusters; i++)
		{
			if (N_remaining > 0)
			{
				N_cluster[i] = average + 1; // number of RA neurons in cluster
				N_remaining -= 1;
			}
			else
			{
				N_cluster[i] = average;
			}

		}

		//for (int i = 0; i < total_num_inh_clusters; i++)
		//	printf("N_cluster = %d\n", N_cluster[i]);
		

		//	if (N_remaining < N_cluster)
		//		N_cluster = N_remaining;
		int ind = 0;
		
		for (int i = 0; i < total_num_inh_clusters; i++)
		{
			// for all RA neurons in the cluster
			for (int j = 0; j < N_cluster[i]; j++)
			{
				do
				{
					close = false;

					xx = xx_I_cluster[i] + generator.random(2*CLUSTER_SIZE) - CLUSTER_SIZE;
					yy = yy_I_cluster[i] + generator.random(2*CLUSTER_SIZE) - CLUSTER_SIZE;

					
					for (int m = 0; m < num_inh_in_cluster; m++)
					{

						double dtoI = distance(xx, yy, xx_I[i*num_inh_in_cluster+m], yy_I[i*num_inh_in_cluster+m]);
					
						// if too close to inhibitory neuron or if beyond cluster radius
						if ((dtoI < MIN_INTERNEURON_DISTANCE)||(dtoI > CLUSTER_SIZE))
						{
							close = true;
						
							//printf("Close; xx = %f; yy = %f; xx_I_cluster[i] = %f; yy_I_cluster[i] = %f; xx_I = %f; yy_I = %f; distance = %f\n", xx, yy, xx_I_cluster[i], yy_I_cluster[i], xx_I[i+m], yy_I[i+m], dtoI);
						}
					}

					// check distance to all other RA neurons in the cluster
					if (!close)
					{
						for (int k = ind; k < ind + j; k++)
						{
							double dtoRA = distance(xx, yy, xx_RA[k], yy_RA[k]);

							if (dtoRA < MIN_INTERNEURON_DISTANCE)
							{
								close = true;
								break;
							}
						}

					}

				} while (close);

				xx_RA.push_back(xx);
				yy_RA.push_back(yy);

				//printf("Pushed xx = %f yy = %f\n", xx, yy);

			}
			ind += N_cluster[i];
		}
	}
}

void PoolParallel::initialize_connections_for_inhibitory_clusters(double Gei_mean, double Gei_var, double Gie_mean, double Gie_var)
{
    if (MPI_rank == 0)
    {

        // connections for HVC(RA) neurons
		for (int i = 0; i < N_RA; i++)
		{
	 		for (int j = 0; j < N_I; j++)
	        {
				// nearby connections
				if (distance(xx_I[j], yy_I[j], xx_RA[i], yy_RA[i]) < CLUSTER_SIZE)
				{

		        	 double G = Gei_mean + generator.normal_distribution() * Gei_var;

		             weights_RA_I_global[i].push_back(G);
		             syn_ID_RA_I_global[i].push_back(j);

				}
			}
		}
		
		
		// distant connections
		int ind = 0;

		for (int i = 0; i < num_inh_clusters*num_inh_clusters; i++)
		{
			for (int j = 0; j < N_cluster[i]; j++)
			{
				int target_cluster = i;

				while (target_cluster == i)
				{
					target_cluster = (int) generator.random(num_inh_clusters*num_inh_clusters);

				}
				
				int target = target_cluster*num_inh_in_cluster + (int) generator.random(2);

		    	double G = Gei_mean + generator.normal_distribution() * Gei_var;

		    	weights_RA_I_global[ind + j].push_back(G);
		    	syn_ID_RA_I_global[ind + j].push_back(target);

			}

			ind += N_cluster[i];

		}
		
		/*
		int ind = 0;

		for (int i = 0; i < N_I; i++)
		{
			int neuron_with_distant_connection = (int) generator.random(N_cluster[i]);
				
			int target = i;

			while (target == i)
			{
				target = (int) generator.random(N_I);
	
			}

		    double G = Gei_mean + generator.normal_distribution() * Gei_var;

		    weights_RA_I_global[ind + neuron_with_distant_connection].push_back(G);
		    syn_ID_RA_I_global[ind + neuron_with_distant_connection].push_back(target);


			ind += N_cluster[i];	

		}
		*/

			
		
		// connections for HVC(I) neurons


		for (int i = 0; i < N_I; i++)
     	{
			for (int j = 0; j < N_RA; j++)
			{
				if (distance(xx_I[i], yy_I[i], xx_RA[j], yy_RA[j]) < CLUSTER_SIZE)
				{

		        	 double G = Gie_mean + generator.normal_distribution() * Gie_var;

		             weights_I_RA_global[i].push_back(G);
		             syn_ID_I_RA_global[i].push_back(j);

				}

			}
		}
	}

}

void PoolParallel::initialize_equal_clusters()
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

                    for (unsigned j = 0; j < xx_I.size(); j++)
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
		// for each inhibitory neuron
		
		int average = N_RA / N_I; // average number of RA neurons per single inhibitory neuron
		int N_remaining = N_RA - average*N_I; // number of RA neurons that need to be distributed in addition to average amount
		
		N_cluster = new int[N_I]; // array with number of RA neurons per cluster

		for (int i = 0; i < N_I; i++)
		{
			if (N_remaining > 0)
			{
				N_cluster[i] = average + 1; // number of RA neurons in cluster
				N_remaining -= 1;
			}
			else
			{
				N_cluster[i] = average;
			}

		}

		//for (int i = 0; i < N_I; i++)
		//	printf("N_cluster = %d\n", N_cluster[i]);

		//	if (N_remaining < N_cluster)
		//		N_cluster = N_remaining;
		for (int i = 0; i < N_I; i++)
		{
			// for all RA neurons in the cluster
			for (int j = 0; j < N_cluster[i]; j++)
			{
				do
				{
					close = false;

					xx = xx_I[i] + generator.random(2*CLUSTER_SIZE) - CLUSTER_SIZE;
					yy = yy_I[i] + generator.random(2*CLUSTER_SIZE) - CLUSTER_SIZE;

					double dtoI = distance(xx, yy, xx_I[i], yy_I[i]);
					
					// if too close to inhibitory neuron or if beyond cluster radius
					if ((dtoI < MIN_INTERNEURON_DISTANCE)||(dtoI > CLUSTER_SIZE))
						close = true;

					// check distance to all other RA neurons in the cluster
					if (!close)
					{
						for (int k = 0; k < j; k++)
						{
							double dtoRA = distance(xx, yy, xx_RA[k], yy_RA[k]);

							if (dtoRA < MIN_INTERNEURON_DISTANCE)
							{
								close = true;
								break;
							}
						}

					}

				} while (close);

				xx_RA.push_back(xx);
				yy_RA.push_back(yy);

			}
			

		}

	
    }
}
void PoolParallel::initialize_clusters()
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

                    for (unsigned j = 0; j < xx_I.size(); j++)
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
		// for each inhibitory neuron
		
		int average = N_RA / N_I; // average number of RA neurons per single inhibitory neuron
		int N_remaining = N_RA; // number of RA neurons that need to be distributed in addition to average amount
		
		int* N_cluster = new int[N_I]; // array with number of RA neurons per cluster

		for (int i = 0; i < N_I; i++)
		{
			N_cluster[i] = 2*average - (int) generator.random(2*average); // number of RA neurons in cluster
			

			if (i == N_I - 1)
			{
				if (N_remaining >= 3*average)
				{
					N_cluster[i] = 2*average;
					N_remaining -= N_cluster[i];
					
					for (int j = 0; j < N_remaining; j++)
					{
						N_cluster[j] += 1;
					}
				}
				else
					N_cluster[i] = N_remaining;
			}

			N_remaining -= N_cluster[i];
			
		}
		for (int i = 0; i < N_I; i++)
			printf("N_cluster = %d\n", N_cluster[i]);

		//	if (N_remaining < N_cluster)
		//		N_cluster = N_remaining;
		for (int i = 0; i < N_I; i++)
		{
			// for all RA neurons in the cluster
			for (int j = 0; j < N_cluster[i]; j++)
			{
				do
				{
					close = false;

					xx = xx_I[i] + generator.random(2*CLUSTER_SIZE) - CLUSTER_SIZE;
					yy = yy_I[i] + generator.random(2*CLUSTER_SIZE) - CLUSTER_SIZE;

					double dtoI = distance(xx, yy, xx_I[i], yy_I[i]);
					
					// if too close to inhibitory neuron or if beyond cluster radius
					if ((dtoI < MIN_INTERNEURON_DISTANCE)||(dtoI > CLUSTER_SIZE))
						close = true;

					// check distance to all other RA neurons in the cluster
					if (!close)
					{
						for (int k = 0; k < j; k++)
						{
							double dtoRA = distance(xx, yy, xx_RA[k], yy_RA[k]);

							if (dtoRA < MIN_INTERNEURON_DISTANCE)
							{
								close = true;
								break;
							}
						}

					}

				} while (close);

				xx_RA.push_back(xx);
				yy_RA.push_back(yy);

			}
			

		}

		delete[] N_cluster;
    }
}

void PoolParallel::initialize_test_allI2RA_connections(double Gie)
{
    for (int i = 0; i < N_I_local; i++)
    {
        for (int j = 0; j < N_RA; j++)
        {
            weights_I_RA_local[i].push_back(Gie);
            syn_ID_I_RA_local[i].push_back(j);
        }

    }

    if (MPI_rank == 0)
    {
        for (int i = 0; i < N_I; i++)
        {
            for (int j = 0; j < N_RA; j++)
            {
                weights_I_RA_global[i].push_back(Gie);
                syn_ID_I_RA_global[i].push_back(j);
            }

        }
    }

}

void PoolParallel::initialize_test_allRA2I_connections(double Gei)
{
    for (int i = 0; i < N_RA_local; i++)
    {
        for (int j = 0; j < N_I; j++)
        {
            weights_RA_I_local[i].push_back(Gei);
            syn_ID_RA_I_local[i].push_back(j);
        }

    }

    if (MPI_rank == 0)
    {
        for (int i = 0; i < N_RA; i++)
        {
            for (int j = 0; j < N_I; j++)
            {
                weights_RA_I_global[i].push_back(Gei);
                syn_ID_RA_I_global[i].push_back(j);
            }

        }
    }

}


void PoolParallel::read_from_file(const char* RA_xy, const char* I_xy, const char* RA_RA_all, const char* RA_RA_active, const char* RA_RA_super, const char* RA_I, const char * I_RA, const char* mature)
{
	if (MPI_rank == 0)
	{
	
		std::ifstream inp_RA_xy, inp_I_xy, inp_RA_RA, inp_RA_RA_super, inp_RA_I, inp_I_RA, inp_mature;
	
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
	
		this->read_all_connections_from_net(RA_RA_all, &weights_global);	
	
		// read RA to RA supersynapses
		this->read_connections_from_net(RA_RA_super, &active_supersynapses_global, &supersynapses_G);	
		this->read_connections_from_net(RA_RA_active, &active_synapses_global, &active_G);
	/*
		for (int i = 0; i < N_RA; i++)
		{
			printf("Neuron %d has %d supersynapses\n", i, (int) active_supersynapses_global[i].size());

			for (int j = 0; j < (int) active_supersynapses_global[i].size(); j++)
			{
				printf("target id: %d\tsynaptic weight: %f\n", active_supersynapses_global[i][j], supersynapses_G[i][j]);
			}

		}
	*/
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
		
		// read maturation indicators

		inp_mature.open(mature, std::ios::binary | std::ios::in);
		inp_mature.read(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));

		for (int i = 0; i < N_RA; i++)
		{	
			inp_mature.read(reinterpret_cast<char *>(&mature_global[i]), sizeof(mature_global[i]));
			
			//if (mature_global[i] > 0)
				//printf("Neuron %d is mature\n", i);
			//else
				//printf("Neuron %d is immature\n", i);
		}

		
		inp_RA_I.close();
		inp_I_RA.close();
		inp_mature.close();
	}
}

void PoolParallel::print_simulation_parameters()
{
	if (MPI_rank == 0)
	{
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

	}


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

	//for (unsigned i = 0; i < items.size(); i++)
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

void PoolParallel::read_connections_from_net(const char* filename, std::vector<unsigned int>** target_id, std::vector<double>** target_G)
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

	//for (unsigned i = 0; i < items.size(); i++)
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

void PoolParallel::send_RARA_connections()
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

		sendcounts_syn_num_RA[0] = N_RA_local;
		displs_syn_num_RA[0] = 0;

		for (int i = 1; i < MPI_size; i++)
		{
			sendcounts_syn_num_RA[i] = N_RA / MPI_size;
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

	// send maturation indicator
	MPI_Scatterv(&mature_global[0], sendcounts_syn_num_RA, displs_syn_num_RA, MPI_INT,
		&mature_local[0], N_RA_local, MPI_INT, 0, MPI_COMM_WORLD);
	

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

		sendcounts_syn_num_RA[0] = N_RA_local;
		sendcounts_syn_num_I[0] = N_I_local;
		displs_syn_num_RA[0] = 0;
		displs_syn_num_I[0] = 0;

		for (int i = 1; i < MPI_size; i++)
		{
			sendcounts_syn_num_RA[i] = N_RA / MPI_size;
			sendcounts_syn_num_I[i] = N_I / MPI_size;
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

void PoolParallel::initialize_connections_for_clusters(double Gei_mean, double Gei_var, double Gie_mean, double Gie_var)
{

    if (MPI_rank == 0)
    {

        // connections for HVC(RA) neurons
		for (int i = 0; i < N_RA; i++)
		{
	 		for (int j = 0; j < N_I; j++)
	        {
				// nearby connections
				if (distance(xx_I[j], yy_I[j], xx_RA[i], yy_RA[i]) < CLUSTER_SIZE)
				{

		        	 double G = Gei_mean + generator.normal_distribution() * Gei_var;

		             weights_RA_I_global[i].push_back(G);
		             syn_ID_RA_I_global[i].push_back(j);

				}
			}
		}
		
		
		// distant connections
		int ind = 0;

		for (int i = 0; i < N_I; i++)
		{
			for (int j = 0; j < N_cluster[i]; j++)
			{
				int target = i;

				while (target == i)
				{
					target = (int) generator.random(N_I);

				}

		    	double G = Gei_mean + generator.normal_distribution() * Gei_var;

		    	weights_RA_I_global[ind + j].push_back(G);
		    	syn_ID_RA_I_global[ind + j].push_back(target);

			}

			ind += N_cluster[i];

		}
		
		/*
		int ind = 0;

		for (int i = 0; i < N_I; i++)
		{
			int neuron_with_distant_connection = (int) generator.random(N_cluster[i]);
				
			int target = i;

			while (target == i)
			{
				target = (int) generator.random(N_I);
	
			}

		    double G = Gei_mean + generator.normal_distribution() * Gei_var;

		    weights_RA_I_global[ind + neuron_with_distant_connection].push_back(G);
		    syn_ID_RA_I_global[ind + neuron_with_distant_connection].push_back(target);


			ind += N_cluster[i];	

		}
		*/

			
		
		// connections for HVC(I) neurons


		for (int i = 0; i < N_I; i++)
     	{
			for (int j = 0; j < N_RA; j++)
			{
				if (distance(xx_I[i], yy_I[i], xx_RA[j], yy_RA[j]) < CLUSTER_SIZE)
				{

		        	 double G = Gie_mean + generator.normal_distribution() * Gie_var;

		             weights_I_RA_global[i].push_back(G);
		             syn_ID_I_RA_global[i].push_back(j);

				}

			}
		}
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

    for (int i = 0; i < N_RA_local; i++)
    {
        HVCRA_local[i].reset();
        spikes_in_trial_soma_local[i].clear();
        spikes_in_trial_dend_local[i].clear();

    }
    for (int i = 0; i < N_I_local; i++)
        HVCI_local[i].reset();
}

void PoolParallel::randomize_after_trial()
{
    for (int i = 0; i < N_RA; i++)
    {
        spikes_in_trial_soma_global[i].clear();
        spikes_in_trial_dend_global[i].clear();
    }

    for (int i = 0; i < N_RA_local; i++)
    {
	HVCRA_local[i].set_to_rest();

        spikes_in_trial_soma_local[i].clear();
        spikes_in_trial_dend_local[i].clear();

    }

    for (int i = 0; i < N_I_local; i++)
        HVCI_local[i].set_to_rest();
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

void PoolParallel::set_generator4neurons()
{
    for (int i = 0; i < N_RA_local; i++)
        HVCRA_local[i].set_noise_generator(&generator);

    for (int i = 0; i < N_I_local; i++)
        HVCI_local[i].set_noise_generator(&generator);
}

void PoolParallel::set_dynamics(double interval, double tS)
{
	size = (int) round(interval/tS) + 1;
	timeStep = tS;
	trial_duration = interval;

    //printf("size = %d\n", size);

	for (int i = 0; i < N_RA_local; i++)
		HVCRA_local[i].set_dynamics(interval, tS);

	for (int i = 0; i < N_I_local; i++)
		HVCI_local[i].set_dynamics(interval, tS);
}


void PoolParallel::enable_motivation_noise()
{
    for (unsigned i = 0; i < N_RA_local; i++)
    	HVCRA_local[i].set_white_noise_distribution_soma(mu_soma, sigma_soma);

    for (unsigned i = 0; i < N_RA_local; i++)
    	HVCRA_local[i].set_white_noise_distribution_dend(mu_dend, sigma_dend);
}

void PoolParallel::disable_motivation_noise()
{
    for (unsigned i = 0; i < N_RA_local; i++)
    	HVCRA_local[i].set_white_noise_distribution_soma(0, sigma_soma);

    for (unsigned i = 0; i < N_RA_local; i++)
    	HVCRA_local[i].set_white_noise_distribution_dend(0, sigma_dend);
}

void PoolParallel::set_white_noise_distribution_soma(double mu, double sigma)
{
    mu_soma = mu;
    sigma_soma = sigma;

    for (unsigned i = 0; i < N_RA_local; i++)
        HVCRA_local[i].set_white_noise_distribution_soma(mu, sigma);
}

void PoolParallel::set_white_noise_distribution_dend(double mu, double sigma)
{
    mu_dend = mu;
    sigma_dend = sigma;

    for (unsigned i = 0; i < N_RA_local; i++)
        HVCRA_local[i].set_white_noise_distribution_dend(mu, sigma);
}

void PoolParallel::set_white_noise_soma()
{
    for (unsigned i = 0; i < N_RA_local; i++)
        HVCRA_local[i].set_white_noise_soma();
}

void PoolParallel::set_white_noise_dend()
{
    for (unsigned i = 0; i < N_RA_local; i++)
        HVCRA_local[i].set_white_noise_dend();
}

void PoolParallel::set_no_noise()
{
    for (unsigned i = 0; i < N_RA_local; i++)
        HVCRA_local[i].set_no_noise();

    for (unsigned i = 0; i < N_I_local; i++)
        HVCI_local[i].set_no_noise();
}

void PoolParallel::set_no_noise_RA()
{
    for (unsigned i = 0; i < N_RA_local; i++)
        HVCRA_local[i].set_no_noise();

}

void PoolParallel::set_no_noise_I()
{
    for (unsigned i = 0; i < N_I_local; i++)
        HVCI_local[i].set_no_noise();

}

void PoolParallel::set_testing_current()
{
    if (MPI_rank == 0)
    {
        double start = 2*T_P + current_injection_time;
        std::function<double (double)> f = std::bind(&training_current, start, _1);

        for (int i = 0; i < 5; i++)
            HVCRA_local[N_TR+i].set_dend_current(f);
    }
}

void PoolParallel::trial(bool training)
{
	int some_RA_neuron_fired_soma_local;
	int some_RA_neuron_fired_dend_local;

	int some_I_neuron_fired_local;
	int some_RA_neuron_fired_soma_global;
	int some_RA_neuron_fired_dend_global;

	int some_I_neuron_fired_global;

	std::vector<unsigned> RA_neurons_fired_soma_local;
	std::vector<unsigned> RA_neurons_fired_soma_global;
	std::vector<unsigned> RA_neurons_fired_soma_realID;

	std::vector<unsigned> RA_neurons_fired_dend_local;
    std::vector<unsigned> RA_neurons_fired_dend_global;
    std::vector<unsigned> RA_neurons_fired_dend_realID;

	std::vector<unsigned> I_neurons_fired;

    trial_number++;
    internal_time = (trial_number - 1) * trial_duration;
    // if training trial
    if (training)
        this->set_training_current();

    //this->set_testing_current();


	bool depolarization = false;
	bool demotivation = false;
	// evolve dynamics
	for (unsigned t = 1; t < size; t++)
	{
		internal_time += timeStep;
		
		if ((t*timeStep <= trial_duration - DEMOTIVATION_WINDOW)&&(!depolarization))
		{
		    this->enable_motivation_noise();
		    depolarization = true;
		}
		
		if ((t*timeStep >= trial_duration - DEMOTIVATION_WINDOW)&&(!demotivation))
		{
		    this->disable_motivation_noise();
		    demotivation = true;
		}

		some_RA_neuron_fired_soma_local = 0;
		some_RA_neuron_fired_soma_global = 0;

		some_RA_neuron_fired_dend_local = 0;
        some_RA_neuron_fired_dend_global = 0;

        some_I_neuron_fired_local = 0;
        some_I_neuron_fired_global = 0;

        RA_neurons_fired_soma_local.clear();
        RA_neurons_fired_soma_global.clear();
        RA_neurons_fired_soma_realID.clear();

        RA_neurons_fired_dend_local.clear();
        RA_neurons_fired_dend_global.clear();
        RA_neurons_fired_dend_realID.clear();

        I_neurons_fired.clear();


		// initialize update arrays
		for (int i = 0; i < N_RA; i++)
		{
            update_g_local[i] = 0.0;
		    update_Ge_RA_local[i] = 0.0;
            update_Gi_RA_local[i] = 0.0;
            update_Ge_RA_global[i] = 0.0;
            update_Gi_RA_global[i] = 0.0;
		}

		for (int i = 0; i < N_I; i++)
		{
            update_Ge_I_local[i] = 0.0;
            update_Ge_I_global[i] = 0.0;
		}


		for (int i = 0; i < N_RA_local; i++)
		{
            int syn_num = active_synapses_local[i].size();

            // set GABA potential
            if ((mature_local[i] > 0)||(syn_num >= N_MATURATION)||(Id_RA_local[i] < N_TR))
                HVCRA_local[i].set_Ei(E_GABA_MATURE);
            else
                //HVCRA_local[i].set_Ei(E_GABA(syn_num));
		HVCRA_local[i].set_Ei(E_GABA_IMMATURE);
            // RK4 step
            HVCRA_local[i].R4_step_no_target_update();

            if (HVCRA_local[i].get_fired_soma())
            {
                some_RA_neuron_fired_soma_local = 1;
                spike_times_soma_local[i] = internal_time;
                spikes_in_trial_soma_local[i].push_back(internal_time);
                RA_neurons_fired_soma_local.push_back(i);
                RA_neurons_fired_soma_realID.push_back(Id_RA_local[i]);


                last_soma_spikes_local[i].push_back(internal_time);

                while ((!last_soma_spikes_local[i].empty())&&
                    (last_soma_spikes_local[i].back() - last_soma_spikes_local[i].front() > LTP_WINDOW))
                    last_soma_spikes_local[i].pop_front();

                //for (int j = 0; j < last_soma_spikes_local[i].size(); j++)
                //{
                  //  printf("My rank = %d; Soma spikes of neuron %d:  spike_time = %f\n", MPI_rank, Id_RA_local[i],
                    //    last_soma_spikes_local[i][j]);

                //}

                //printf("My rank = %d; Soma neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], spike_times_soma_local[i]);

                //RA_neurons_fired_soma_realID.push_back(Id_RA_local[i]);
            }

            if (HVCRA_local[i].get_fired_dend())
            {
                some_RA_neuron_fired_dend_local = 1;
                spike_times_dend_local[i] = internal_time;
                spikes_in_trial_dend_local[i].push_back(internal_time);
                //printf("My rank = %d; RA neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], spike_times_local[i]);
                RA_neurons_fired_dend_local.push_back(i);
                RA_neurons_fired_dend_realID.push_back(Id_RA_local[i]);

                //printf("My rank = %d; Dend neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], spike_times_dend_local[i]);

            }

		}

		for (int i = 0; i < N_I_local; i++)
		{
            HVCI_local[i].R4_step_no_target_update();

            if (HVCI_local[i].get_fired())
            {
                //printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
                some_I_neuron_fired_local = 1;
                I_neurons_fired.push_back(i);
            }
		}

		// fill update array for targets
		// for fired RA neurons
		if (some_RA_neuron_fired_soma_local == 1)
		{
            // loop over all fired neurons
            for (int i = 0; i < RA_neurons_fired_soma_local.size(); i++)
            {
                int fired_ID = RA_neurons_fired_soma_local[i]; // local ID of fired neuron

                // loop over all inhibitory targets of fired neurons
                int num_I_targets = syn_ID_RA_I_local[fired_ID].size();
                for (int j = 0; j < num_I_targets; j++)
                {
                    int syn_ID = syn_ID_RA_I_local[fired_ID][j];
                    update_Ge_I_local[syn_ID] += weights_RA_I_local[fired_ID][j];

                }

                // loop over all excitatory targets
                int num_RA_targets = active_synapses_local[fired_ID].size();
                for (int j = 0; j < num_RA_targets; j++)
                {
                    int syn_ID = active_synapses_local[fired_ID][j];
                    update_Ge_RA_local[syn_ID] += weights_local[fired_ID][syn_ID];
                }

                // if neuron is saturated apply LTD only to supersynapses and supply glutamate only to supersynapses
                if (active_supersynapses_local[fired_ID].size() == Nss)
                {
                    for (int k = 0; k < active_supersynapses_local[fired_ID].size(); k++)
                    {
                        int target_ID = active_supersynapses_local[fired_ID][k];

						//update_g_local[target_ID] += g_KICK; // update glutamate of supersynaptic target

                        double dt = internal_time - spike_times_dend_global[target_ID];

                        if ((Id_RA_local[fired_ID] != target_ID)&&(dt < LTD_WINDOW))
                        {
                        
         //                   printf("LTD from neuron %d onto %d; somatic spike: %f; dendritic spike: %f\n", Id_RA_local[fired_ID], target_ID,
			    	//spike_times_soma_local[fired_ID], spike_times_dend_global[target_ID]);
                            LTD(weights_local[fired_ID][target_ID], dt);
                            update_synapse(fired_ID, target_ID, weights_local[fired_ID][target_ID]);
                        }
                    }
                }
                // if not saturated apply LTD rule with the last dendritic spike of all neurons and add glutamate to all neuron except the fired one
                else
                {
                    for (int j = 0; j < N_RA; j++)
                    {
                        double dt = internal_time - spike_times_dend_global[j];

                        if (Id_RA_local[fired_ID] != j)
                        {
							//update_g_local[j] += g_KICK; // update glutamate

							if (dt < LTD_WINDOW)
							{
                           // printf("LTD from neuron %d onto %d; somatic spike: %f; dendritic spike: %f\n", Id_RA_local[fired_ID], j,
		//	    	spike_times_soma_local[fired_ID], spike_times_dend_global[j]);
                           		LTD(weights_local[fired_ID][j], dt);
                           		update_synapse(fired_ID, j, weights_local[fired_ID][j]);
							}
                        }
                    }
                }

            }
		}

		// for fired I neurons
		if (some_I_neuron_fired_local == 1)
		{
            // loop over all fired neurons
            for (int i = 0; i < I_neurons_fired.size(); i++)
            {
                int fired_ID = I_neurons_fired[i];

                int num_RA_targets = syn_ID_I_RA_local[fired_ID].size();
                // loop over all targets of fired neuron
                for (int j = 0; j < num_RA_targets; j++)
                {
                    int syn_ID = syn_ID_I_RA_local[fired_ID][j];
                    update_Gi_RA_local[syn_ID] += weights_I_RA_local[fired_ID][j];
                    //printf("Rank = %d; i = %d; update_Gi_RA_local[%d] = %f; weights_I_RA_local[%d][%d] = %f\n", MPI_rank, i, syn_ID,
                     //   update_Gi_RA_local[syn_ID], weights_I_RA_local[fired_ID][j], fired_ID, j);
                }

            }
		}

        // get if any neurons fired in some process
        MPI_Allreduce(&some_RA_neuron_fired_soma_local, &some_RA_neuron_fired_soma_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&some_RA_neuron_fired_dend_local, &some_RA_neuron_fired_dend_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&some_I_neuron_fired_local, &some_I_neuron_fired_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        if (some_I_neuron_fired_global > 0)
        {
            // sum update array and send to all processes

            MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            for (int i = 0; i < N_RA_local; i++)
            {
                HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
            }
        }

        //if (some_RA_neuron_fired_global == 1)
        //    printf("Rank %d; some_RA_neuron_fired_global: %d\n", MPI_rank, some_RA_neuron_fired_global);

        // if somatic compartment of any neuron in the pool fired, update synaptic conductances
        if (some_RA_neuron_fired_soma_global > 0)
        {
            // sum all update arrays and send to all processes

            MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            //MPI_Allreduce(&update_g_local[0], &update_g_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            // now update excitatory conductances and glutamate input of all neurons
            for (int i = 0; i < N_RA_local; i++)
            {
                HVCRA_local[i].raiseE(update_Ge_RA_global[Id_RA_local[i]]); // update conductance
				//HVCRA_local[i].raiseGlutamate(update_g_global[Id_RA_local[i]]); // update glutamate
				//if (update_g_global[Id_RA_local[i]]!=0)
				//	HVCRA_local[i].set_glutamate(g_KICK);
            }

            for (int i = 0; i < N_I_local; i++)
            {
                HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);
            }

            // send fired neurons
            /*
            int num_RA_fired_soma_local = RA_neurons_fired_soma_local.size();
            int num_RA_fired_soma_global;

            int* recvcounts = new int[MPI_size];
            int* displs = new int[MPI_size];
            MPI_Status status;

            // get total number of fired neurons
            MPI_Allreduce(&num_RA_fired_soma_local, &num_RA_fired_soma_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            RA_neurons_fired_soma_global.resize(num_RA_fired_soma_global);
            */
            //printf("Rank %d; fired RA num: %d\n", MPI_rank, num_RA_fired_local);

            //if (MPI_rank == 0)
            //{
             //   printf("Master; fired RA num global: %d\n", num_RA_fired_global);
             //   printf("Master; RA_neurons_fired_global.size(): %d\n", RA_neurons_fired_global.size());
            //}

            // get array with number of fired neurons in each process
            /*
            MPI_Allgather(&num_RA_fired_soma_local, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, MPI_COMM_WORLD);

            displs[0] = 0;
            for (int i = 1; i < MPI_size; i++)
            {
                displs[i] = displs[i-1] + recvcounts[i-1];
            }


            //for (int i = 0; i < RA_neurons_fired_realID.size(); i++)
            //        printf("Rank %d; fired RA neuron: %d\n", MPI_rank, RA_neurons_fired_realID[i]);

            // get fired neurons
            MPI_Allgatherv(&RA_neurons_fired_soma_realID[0], num_RA_fired_soma_local, MPI_INT,
                &RA_neurons_fired_soma_global[0], recvcounts, displs, MPI_INT, MPI_COMM_WORLD);

            delete [] recvcounts;
            delete [] displs;
            */
            //if (MPI_rank == 0)
            //{
                //printf("Master; RA_neurons_fired_global.size(): %d\n", RA_neurons_fired_global.size());
                //for (int i = 0; i < RA_neurons_fired_global.size(); i++)
                  //  printf("Rank %d; fired RA neuron: %d\n", MPI_rank, RA_neurons_fired_global[i]);
            //}

            // change spike times global array
            /*
            for (int i = 0; i < RA_neurons_fired_soma_global.size(); i++)
            {
                spike_times_soma_global[RA_neurons_fired_soma_global[i]] = internal_time;
                last_soma_spikes_global[RA_neurons_fired_soma_global[i]].push_back(internal_time);

                // if array already has > NUM_SOMA_SPIKES delete the earliest one
                if (last_soma_spikes_global[i].size() > NUM_SOMA_SPIKES)
                    last_soma_spikes_global[i].pop_front();

            }
            */

        }

        // if dendritic compartment of any neuron in the pool fired, update weights
        if (some_RA_neuron_fired_dend_global > 0)
        {
            // send all fired neurons
            int num_RA_fired_dend_local = RA_neurons_fired_dend_local.size();
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

            // send spike times
            /*
            recvcounts[0] = N_RA/MPI_size + N_RA%MPI_size;
            displs[0] = 0;
            for (int i = 1; i < MPI_size; i++)
            {
                recvcounts[i] = N_RA/MPI_size;
                displs[i] = displs[i-1] + recvcounts[i-1];
            }

            MPI_Allgatherv(&spike_times_local[0], N_RA_local, MPI_DOUBLE, &spike_times_global[0],
                recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
            */
            delete [] recvcounts;
            delete [] displs;

            /*
            if (MPI_rank == 0)
            {   for (int i = 0; i < N_RA; i++)
                    printf("neuron %d; spike_time_dend = %f\n", i ,spike_times_dend_global[i]);
            }
            */
            // apply LTP rule for last NUM_SOMA_SPIKES and dendritic spike of dired neurons
            for (int i = 0; i < N_RA_local; i++)
            {
                // if neuron is saturated apply LTP only if dendritic spike occured in supersynapse
                if (active_supersynapses_local[i].size() == Nss)
                {
                    for (int j = 0; j < RA_neurons_fired_dend_global.size(); j++)
                    {
                        int fired_ID = RA_neurons_fired_dend_global[j];

                        std::vector<unsigned>::iterator pos = std::find(active_supersynapses_local[i].begin(),
                                        active_supersynapses_local[i].end(), fired_ID);

                        if (pos!=active_supersynapses_local[i].end())
                        {
                            for  (int k = 0; k < last_soma_spikes_local[i].size(); k++)
                            {
                                double dt = internal_time - last_soma_spikes_local[i][k];

                                if (dt < LTP_WINDOW)
                                {
                                    LTP(weights_local[i][fired_ID], dt);
         //                           printf("LTP from neuron %d onto %d; somatic spike at %f; dendritic spike at %f\n", Id_RA_local[i], fired_ID,
             //                                               last_soma_spikes_local[i][k], spike_times_dend_global[fired_ID]);
                                    update_synapse(i, fired_ID, weights_local[i][fired_ID]);
                                }
                            }
                        }
                    }

                }
                // if not saturated apply LTP for all dendritic spikes
                else
                {
                    for (int j = 0; j < RA_neurons_fired_dend_global.size(); j++)
                    {
                        int fired_ID = RA_neurons_fired_dend_global[j];
                        // don't allow self-to-self connections
                        if (fired_ID != Id_RA_local[i])
                        {
                            // loop over last somatic spikes
                            for  (int k = 0; k < last_soma_spikes_local[i].size(); k++)
                            {
                                double dt = internal_time - last_soma_spikes_local[i][k];

                                if (dt < LTP_WINDOW)
                                {
           //                         printf("LTP from neuron %d onto %d; somatic spike at %f; dendritic spike: %f\n", Id_RA_local[i], fired_ID,
                            //                                last_soma_spikes_local[i][k], spike_times_dend_global[fired_ID]);
                                    LTP(weights_local[i][fired_ID], dt);
                                    update_synapse(i, fired_ID, weights_local[i][fired_ID]);
                                }
                            }
                        }
                    }
                }

            }

            /*
            // loop over all fired neurons
            for (int i = 0; i < RA_neurons_fired_dend_global.size(); i++)
            {
                int fired_ID = RA_neurons_fired_dend_global[i];
                int rank, ind;

                // calculate process and index in array where fired neuron is located
                if (fired_ID < N_RA/MPI_size + N_RA%MPI_size)
                {
                    rank = 0;
                    ind = fired_ID;
                }
                else
                {
                    rank = (fired_ID - N_RA%MPI_size) / (N_RA/MPI_size);
                    ind = (fired_ID - N_RA/MPI_size - N_RA%MPI_size) % (N_RA/MPI_size);
                }
                //printf("Process with fired neuron: %d; ind = %d\n", rank, ind);

                // if process has fired neuron
                if (MPI_rank == rank)
                {
                    for (int k = 0; k < N_RA_local; k++)
                    {
                        int real_ID = Id_RA_local[k];

                        if (k==ind)
                        {
                          // update weights from fired neuron onto all other neurons
                          // if saturated apply STDP rule only to active targets which are supersynapses

                            if (active_supersynapses_local[ind].size() == Nss)
                            {
                                for (int j = 0; j < active_synapses_local[ind].size(); j++)
                                {
                                    int syn_ID = active_synapses_local[ind][j];

                                    STDP(ind, syn_ID, spike_times_global[fired_ID], spike_times_global[syn_ID]);
                                }
                            }
                            // if not saturated apply STDP to all neurons
                            else
                            {
                                for (int j = 0; j < N_RA; j++)
                                {
                                    STDP(ind, j, spike_times_global[fired_ID], spike_times_global[j]);
                                }
                            }

                        }

                            // update weights from all other neurons onto fired
                        else
                        {
                            // if neuron is saturated, then check if fired neuron is among supersynapses
                            if (active_supersynapses_local[k].size() == Nss)
                            {
                                std::vector<unsigned>::iterator pos = std::find(active_supersynapses_local[k].begin(),
                                                    active_supersynapses_local[k].end(), fired_ID);

                                // if fired neuron is supersynapse, apply STDP rules
                                if (pos != active_supersynapses_local[k].end())
                                    STDP(k, fired_ID, spike_times_global[real_ID], spike_times_global[fired_ID]);
                            }
                            // if not saturated just apply STDP rule
                            else
                                STDP(k, fired_ID, spike_times_global[k], spike_times_global[fired_ID]);

                        }

                    }

                }

                // if process doesn't have fired neuron just update weights from all other neurons onto fired
                else
                {
                    for (int k = 0; k < N_RA_local; k++)
                    {
                        int real_ID = Id_RA_local[k];
                        // if neuron is saturated, then check if fired neuron is among supersynapses
                        if (active_supersynapses_local[k].size() == Nss)
                        {
                            std::vector<unsigned>::iterator pos = std::find(active_supersynapses_local[k].begin(),
                                            active_supersynapses_local[k].end(), fired_ID);

                        // if fired neuron is supersynapse, apply STDP rules
                            if (pos != active_supersynapses_local[k].end())
                                    STDP(k, fired_ID, spike_times_global[real_ID], spike_times_global[fired_ID]);
                        }
                        // if not saturated just apply STDP rule
                        else
                        {
                            STDP(k, fired_ID, spike_times_global[real_ID], spike_times_global[fired_ID]);

                        }
                    }
                }


            }

            */
            // check if we need axon remodeling

            for (int i = 0; i < N_RA_local; i++)
            {
                if ((active_supersynapses_local[i].size()==Nss)&&(!remodeled_local[i]))
                {
					this->axon_remodeling(i);
					mature_local[i] = 1;
				}
            }
        }


        MPI_Barrier(MPI_COMM_WORLD);
    }
    this->potentiation_decay();
    //printf("After potentiation decay")
    this->update_all_synapses();
	//printf("internal time = %f\n", internal_time);
	//printf("t*timeStep = %f\n", t*timeStep);
}




void PoolParallel::STDP(int i, int j, double ti, double tj)
{
    // if i spiked before j, apply LTP i to j connection, LTD else


    if (Id_RA_local[i]!=j)
    {
        double dt = fabs(ti - tj);
        //printf("Rank %d; from %d onto %d; t1 = %f, t2 = %f\n", MPI_rank, Id_RA_local[i], j, ti, tj);
        if (ti<=tj)
            LTP(weights_local[i][j], dt);
        else
            LTD(weights_local[i][j], dt);

        update_synapse(i, j, weights_local[i][j]);
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
    for (int j = 0; j < active_synapses_local[i].size(); j++)
    {
        int syn_ID = active_synapses_local[i][j];

        // check if active synapse is among supersynapses
        std::vector<unsigned>::iterator pos = std::find(active_supersynapses_local[i].begin(),
                                active_supersynapses_local[i].end(), syn_ID);
        // if not erase it
        if (pos == active_supersynapses_local[i].end())
        {
            active_synapses_local[i].erase(active_synapses_local[i].begin() + j);
            active_local[i][syn_ID] = false;
        }
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
            this->update_synapse(i, j, weights_local[i][j]);
        }
    }

}

void PoolParallel::update_synapse(int i, int j, double w)
{
    if ((w > ACTIVATION)&&(!active_local[i][j])&&(!remodeled_local[i]))
    {
       // printf("Activated synapse from %d onto %d; weight = %f\n", Id_RA_local[i], j, w);
        active_local[i][j] = true;
        active_synapses_local[i].push_back(j);
    }

    if ((w > SUPERSYNAPSE_THRESHOLD)&&(!supersynapses_local[i][j])&&(active_supersynapses_local[i].size()<Nss))
    {
       // printf("Activated supersynapse from %d onto %d; weight = %f\n", Id_RA_local[i], j, w);
        supersynapses_local[i][j] = true;
        active_supersynapses_local[i].push_back(j);
    }

    if ((w < SUPERSYNAPSE_THRESHOLD)&&(supersynapses_local[i][j]))
    {
       // printf("Deactivated supersynapse from %d onto %d; weight = %f\n", Id_RA_local[i], j, w);
        supersynapses_local[i][j] = false;
        remodeled_local[i] = false;
        std::vector<unsigned>::iterator pos = std::find(active_supersynapses_local[i].begin(),
                                                    active_supersynapses_local[i].end(), j);

        if (pos!= active_supersynapses_local[i].end())
            active_supersynapses_local[i].erase(pos);
        else
            printf("Supersynapse to be erased is not found!\n");
    }

    if ((w < ACTIVATION)&&(active_local[i][j]))
    {
      //  printf("Deactivated synapse from %d onto %d; w = %f\n", Id_RA_local[i], j, w);
        active_local[i][j] = false;

        std::vector<unsigned>::iterator pos = std::find(active_synapses_local[i].begin(),
                                                    active_synapses_local[i].end(), j);

        if (pos!= active_synapses_local[i].end())
            active_synapses_local[i].erase(pos);
        else
            printf("Active synapse to be erased is not found!\n");


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
		prob = A_RA2I * exp(-d / LAMBDA_RA2I) + B_RA2I * exp(-(d - MEAN_RA2I)*(d - MEAN_RA2I)/(2*SIGMA_RA2I*SIGMA_RA2I));
		//printf("p = %f\n", prob);
 	}
	return prob;
}

double PoolParallel::p_RA2I_distant(int i_RA, int j_I)
{
	double prob;
    double d;
    d = distance(xx_RA[i_RA], yy_RA[i_RA], xx_I[j_I], yy_I[j_I]);

	if (d < MIN_INTERNEURON_DISTANCE)
		return 0;
	else
	{
		prob = B_RA2I * exp(-(d - MEAN_RA2I)*(d - MEAN_RA2I)/(2*SIGMA_RA2I*SIGMA_RA2I));
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
		prob = CONNECT_CONST_I2RA * exp(-d / LAMBDA_I2RA);
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

double PoolParallel::E_GABA(int n)
{
    return E_GABA_MATURE + (E_GABA_IMMATURE - E_GABA_MATURE) * exp(-n/(N_MATURATION));
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

    int *recvcounts = new int[MPI_size];
    int *displs = new int[MPI_size];

    if (MPI_rank == 0)
    {
        recvcounts[0] = N_RA_local;
        displs[0] = 0;

        for (int i = 1; i < MPI_size; i++)
        {
            recvcounts[i] = N_RA / MPI_size;
            displs[i] = displs[i-1] + recvcounts[i-1];
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

    MPI_Gatherv(&supersyn_sizes_local[0], N_RA_local, MPI_INT,
        &supersyn_sizes_global[0], recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gatherv(&syn_sizes_local[0], N_RA_local, MPI_INT,
        &syn_sizes_global[0], recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);


    MPI_Gatherv(&spike_num_soma_local[0], N_RA_local, MPI_INT,
        &spike_num_soma_global[0], recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gatherv(&spike_num_dend_local[0], N_RA_local, MPI_INT,
        &spike_num_dend_global[0], recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

	// gather maturation indicator
    MPI_Gatherv(&mature_local[0], N_RA_local, MPI_INT,
        &mature_global[0], recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
    
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
        int offset = N_RA / MPI_size;


        for (int i = 1; i < MPI_size; i++)
        {
            for (int j = 0; j < offset; j++)
            {
                int count;

                MPI_Recv(&weights_global[N_RA_local + (i-1)*offset + j][0],
                                        N_RA, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_DOUBLE, &count);
                //printf("Recv weights; from i = %d  count = %d\n", i, count);
                if (supersyn_sizes_global[N_RA_local + (i-1)*offset + j] != 0)
                {
                    MPI_Recv(&active_supersynapses_global[N_RA_local + (i-1)*offset + j][0],
                        supersyn_sizes_global[N_RA_local + (i-1)*offset + j], MPI_UNSIGNED, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_UNSIGNED, &count);
                    //printf("Recv supersynapses; from i = %d  count = %d\n", i, count);
                }

                if (syn_sizes_global[N_RA_local + (i-1)*offset + j] != 0)
                {
                    MPI_Recv(&active_synapses_global[N_RA_local + (i-1)*offset + j][0],
                        syn_sizes_global[N_RA_local + (i-1)*offset + j], MPI_UNSIGNED, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_UNSIGNED, &count);
                    //printf("Recv synapses; from i = %d  count = %d\n", i, count);
                }

                if (spike_num_soma_global[N_RA_local + (i-1)*offset + j] != 0)
                {
                    MPI_Recv(&spikes_in_trial_soma_global[N_RA_local + (i-1)*offset + j][0],
                        spike_num_soma_global[N_RA_local + (i-1)*offset + j], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv spikes in trial; from i = %d  count = %d\n", i, count);
                }

                if (spike_num_dend_global[N_RA_local + (i-1)*offset + j] != 0)
                {
                    MPI_Recv(&spikes_in_trial_dend_global[N_RA_local + (i-1)*offset + j][0],
                        spike_num_dend_global[N_RA_local + (i-1)*offset + j], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

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
                        supersyn_sizes_local[i], MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);

            if (syn_sizes_local[i] != 0)
                MPI_Send(&active_synapses_local[i][0],
                        syn_sizes_local[i], MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);

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


    delete [] recvcounts;
    delete [] displs;
    delete [] supersyn_sizes_global;
    delete [] supersyn_sizes_local;
    delete [] syn_sizes_global;
    delete [] syn_sizes_local;
    delete [] spike_num_soma_local;
    delete [] spike_num_soma_global;
    delete [] spike_num_dend_local;
    delete [] spike_num_dend_global;

}

void PoolParallel::get_neuronRA_location(unsigned n, int* rank, int* shift)
{
    int master_N_RA = N_RA / MPI_size + N_RA % MPI_size;

    if (n < master_N_RA)
    {
        *rank = 0;
	*shift = n;
    }
    else
    {
        *rank = (n - master_N_RA) / (N_RA / MPI_size) + 1;
	*shift = (n - master_N_RA) % (N_RA / MPI_size);
    }


}

void PoolParallel::get_neuronI_location(unsigned n, int* rank, int* shift)
{
    int master_N_I = N_I / MPI_size + N_I % MPI_size;

    if (n < master_N_I)
    {
        *rank = 0;
	*shift = n;
    }
    else
    {
        *rank = (n - master_N_I) / (N_I / MPI_size) + 1;
        *shift = (n - master_N_I) % (N_I / MPI_size);
    }


}
void PoolParallel::statistics()
{

    int local_pool_spikes_soma = 0;
    int total_pool_spikes_soma = 0;
    int local_pool_spikes_dend = 0;
    int total_pool_spikes_dend = 0;

    std::vector<unsigned> training_I_targets;
    std::vector<unsigned> training_I_targets_global;
    
    for (int i = 0; i < N_RA_local; i++)
    {
    // get statistics for all pool neurons first
        if (Id_RA_local[i] >= N_TR)
	{
	    //printf("Rank = %d, ID_RA_local[i] = %d\n", MPI_rank, Id_RA_local[i]);
	    local_pool_spikes_soma += HVCRA_local[i].get_spike_number_soma();
	    local_pool_spikes_dend += HVCRA_local[i].get_spike_number_dend();
	}
	else
    // get spikes of the neurons that receive inhibitory kicks due to training neurons
	    for (int j = 0; j < syn_ID_RA_I_local[i].size(); j++)
	    {
                training_I_targets.push_back(syn_ID_RA_I_local[i][j]);

	    }
    }
    
    MPI_Reduce(&local_pool_spikes_soma, &total_pool_spikes_soma, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_pool_spikes_dend, &total_pool_spikes_dend, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    int size_I_targets = training_I_targets.size();

    if (MPI_rank == 0)
    {
        printf("\nTraining_I_targets: \n");
       
        for (int i = 0; i < training_I_targets.size(); i++)
	    printf("%u\t", training_I_targets[i]);

        printf("\n");

	size_I_targets = training_I_targets.size();
    }

    // Send training I targets to all processes
    
    MPI_Bcast(&size_I_targets, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    training_I_targets.resize(size_I_targets);

    MPI_Bcast(&training_I_targets[0], size_I_targets, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    //printf("Rank = %d; training_I_targets.last() = %u\n", MPI_rank, training_I_targets.back());

    std::vector<unsigned> training_RA_targets, training_RA_targets_global;

    int size_RA_targets, size_RA_targets_global;

    int local_I_targets_spikes = 0;
	int total_I_targets_spikes = 0;
	int local_I_spikes = 0;
	int total_I_spikes = 0;

	// get total number of spikes of all inhibitory neurons

	for (int i = 0; i < N_I_local; i++)
		local_I_spikes += HVCI_local[i].get_spike_number();

	MPI_Reduce(&local_I_spikes, &total_I_spikes, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	// get number of spikes for inhibitory targets of training neurons and get IDs of RA neurons that are targeted by these interneurons

    for (int i = 0; i < training_I_targets.size(); i++)
    {
    	int rank; // rank of the process with I neuron
	int shift; // location of the neuron in the local array

        this->get_neuronI_location(training_I_targets[i], &rank, &shift);
  
        if (MPI_rank == rank)
		{
            for (int j = 0; j < syn_ID_I_RA_local[shift].size(); j++)
	    	{
                training_RA_targets.push_back(syn_ID_I_RA_local[shift][j]);

	    	}
			local_I_targets_spikes += HVCI_local[shift].get_spike_number();
		}
    }

    size_RA_targets = training_RA_targets.size();

    //printf("Rank = %d; size_RA_targets = %d\n", MPI_rank, size_RA_targets);

    MPI_Allreduce(&size_RA_targets, &size_RA_targets_global, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    
	//printf("Rank = %d; local_I_targets_spikes = %d\n", MPI_rank, local_I_targets_spikes);

	MPI_Reduce(&local_I_targets_spikes, &total_I_targets_spikes, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    //printf("Rank = %d; total size = %d\n", MPI_rank, size_RA_targets_global);


    training_RA_targets_global.resize(size_RA_targets_global);

    int* displs = new int[MPI_size];
    int* recvcounts = new int[MPI_size];

    MPI_Allgather(&size_RA_targets, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, MPI_COMM_WORLD);

    displs[0] = 0;
    for (int i = 1; i < MPI_size; i++)
    {
        displs[i] = displs[i-1] + recvcounts[i-1];
    }

    //if (MPI_rank == 0)
    //{
     //   printf("Master; recvcounts[MPI_size-1] = %d\n", recvcounts[MPI_size-1]);
    //}

    MPI_Allgatherv(&training_RA_targets[0], size_RA_targets, MPI_UNSIGNED, &training_RA_targets_global[0], recvcounts, displs, MPI_UNSIGNED, MPI_COMM_WORLD);

    delete [] recvcounts;
    delete [] displs;

    /*if (training_RA_targets.size() > 0)
    {
        printf("Training_RA_targets: \n");
       
        for (int i = 0; i < training_RA_targets.size(); i++)
	    printf("%u\t", training_RA_targets[i]);

    	 printf("\n");
    }*/

    // get rid of similar neurons
    std::vector<unsigned>::iterator iter;

    std::sort(training_RA_targets_global.begin(), training_RA_targets_global.end());

    iter = std::unique(training_RA_targets_global.begin(), training_RA_targets_global.end());
    training_RA_targets_global.resize(std::distance(training_RA_targets_global.begin(), iter));

    if (MPI_rank == 0)
    {

        printf("\nMaster process; Training_RA_targets: \n");
       
        for (int i = 0; i < training_RA_targets_global.size(); i++)
	    printf("%u\t", training_RA_targets_global[i]);

    	 printf("\n");

    }
    

    int local_kicked_spikes_soma = 0;
    int total_kicked_spikes_soma = 0;
    int local_kicked_spikes_dend = 0;
    int total_kicked_spikes_dend = 0;


    for (int i = 0; i < training_RA_targets_global.size(); i++)
    {
    	int rank; // rank of the process with RA neuron
	int shift; // location of the neuron in the local array

        this->get_neuronRA_location(training_RA_targets_global[i], &rank, &shift);
 	//printf("My rank = %d ; rank = %d; shift = %d\n", MPI_rank, rank, shift); 
        if (MPI_rank == rank)
	{
            local_kicked_spikes_soma += HVCRA_local[shift].get_spike_number_soma();	
            local_kicked_spikes_dend += HVCRA_local[shift].get_spike_number_dend();
	
	}
    }

    MPI_Reduce(&local_kicked_spikes_soma, &total_kicked_spikes_soma, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
    MPI_Reduce(&local_kicked_spikes_dend, &total_kicked_spikes_dend, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    //printf("Rank = %d; local_kicked_spikes_soma = %d\n", MPI_rank, local_kicked_spikes_soma);
    
    bool related_spike = false;
	int related_spikes = 0;

    if (MPI_rank == 0)
    {
	printf("\nInhibitory neurons:\n");
	printf("Total number of inhibitory targets of training neurons: %d\n", training_I_targets.size());
	printf("Total number of interneurons: %d\n", N_I);
	printf("Total number of spikes of all interneurons excluding targets of training neurons: %d\n", total_I_spikes - total_I_targets_spikes);
	printf("Total number of spikes of inhibitory targets: %d\n", total_I_targets_spikes);
	printf("Average number of interneuron (excluding targets) spikes per neuron per trial: %f\n", (double)(total_I_spikes - total_I_targets_spikes) / (trial_number*(N_I - training_I_targets.size())));
	printf("Average number of inhibitory target spikes per neuron per trial: %f\n", (double) total_I_targets_spikes / (trial_number * training_I_targets.size()));
       printf("\nExcitatory neurons:\n");
	   printf("Total number of not kicked neurons: %d\n", N_RA - N_TR - training_RA_targets_global.size());
	printf("Total number of kicked neurons: %d\n", training_RA_targets_global.size());
	//printf("Total somatic spikes of pool neurons: %d\n", total_pool_spikes_soma);
        //printf("Total dendritic spikes of pool neurons: %d\n", total_pool_spikes_dend);
	printf("Total somatic spikes of not kicked neurons: %d\n", total_pool_spikes_soma - total_kicked_spikes_soma);
	printf("Total dendritic spikes of not kicked neurons: %d\n", total_pool_spikes_dend - total_kicked_spikes_dend);
	printf("Total somatic spikes of kicked neurons: %d\n", total_kicked_spikes_soma);
        printf("Total dendritic spikes of kicked neurons: %d\n", total_kicked_spikes_dend);
		
		for (int i = 0; i < training_RA_targets_global.size(); i++)
		{

			int RA_ID = training_RA_targets_global[i];
			for (int j = 0; j < spikes_in_trial_dend_global[RA_ID].size(); j++)
			{
				double time_difference = spikes_in_trial_dend_global[RA_ID][j] - current_injection_time - (trial_number-1)*trial_duration;
				
				if ((time_difference < DELAY_WINDOW)&&(time_difference > 0))
				{
					printf("Spike related to the kick. Neuron: %d. Time difference: %f\n", RA_ID, time_difference);
					
					delay += time_difference;
					num_related_spikes += 1;
					related_spike = true;
				}
			}
			//printf("\n");

		}
		
		//printf("\nCurrent injection time: %f\n", current_injection_time + (trial_number - 1)*trial_duration);
		printf("Number of related dendritic spikes: %d\n", num_related_spikes);
		printf("Average time delay between kick and spike: %f\n", delay / num_related_spikes);
		printf("\n");
	

	printf("Average number of dendritic spikes per trial per not kicked neuron: %f\n", (double)
			(total_pool_spikes_dend - total_kicked_spikes_dend) / (trial_number*(N_RA - N_TR - training_RA_targets_global.size())));

 
	printf("Average number of dendritic spikes per kicked neuron: %f\n", (double)
			(total_kicked_spikes_dend) / (trial_number*training_RA_targets_global.size()));
	
	printf("\n");
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

        //for (unsigned i = 0; i < N_RA; i++)
        //    printf("spike_times[%d] = %f\n", i, spike_times[i]);
        // write spike times
        for (unsigned i = 0; i < N_RA; i++)
        {
            int spike_array_size = spikes_in_trial_soma_global[i].size();
            //printf("Neuron %d; number of somatic spikes in trial: %d\n", i, spike_array_size);
            out.write(reinterpret_cast<char *>(&spike_array_size), sizeof(int));

            for (int j = 0; j < spike_array_size; j++)
	    {
                //out.write(reinterpret_cast<char *>(&spikes_in_trial_soma_global[i][j]), sizeof(double));
	        double relative_spike_time = spikes_in_trial_soma_global[i][j] - (trial_number - 1) * trial_duration;
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

        //for (unsigned i = 0; i < N_RA; i++)
        //    printf("spike_times[%d] = %f\n", i, spike_times[i]);
        // write spike times
        for (unsigned i = 0; i < N_RA; i++)
        {
            int spike_array_size = spikes_in_trial_dend_global[i].size();
            //printf("Neuron %d; number of dendritic spikes in trial: %d\n", i, spike_array_size);
            out.write(reinterpret_cast<char *>(&spike_array_size), sizeof(int));
	    
            for (int j = 0; j < spike_array_size; j++)
	    {
                //out.write(reinterpret_cast<char *>(&spikes_in_trial_dend_global[i][j]), sizeof(double));
	        double relative_spike_time = spikes_in_trial_dend_global[i][j] - (trial_number - 1) * trial_duration;
        	out.write(reinterpret_cast<char *>(&relative_spike_time), sizeof(double));
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
    int master_N_RA = N_RA/MPI_size + N_RA%MPI_size;
    //printf("master_N_RA = %d; n = %d\n", master_N_RA, n);

    if (n >= N_RA)
    {
        printf("Selected neuron ID doesn't exist in the pool! Instead wrote neuron 0 to the file.\n");

        if (MPI_rank == 0)
            HVCRA_local[0].writeToFile(filename);
    }
    else
    {
        if (n < master_N_RA)
        {
            //printf("master_N_RA = %d; n = %d\n", master_N_RA, n);

            if (MPI_rank == 0)
            {
                HVCRA_local[n].writeToFile(filename);
                //printf("Master. master_N_RA = %d; n = %d\n", master_N_RA, n);
            }
        }
        else
        {
            int rank = (n - master_N_RA) / (N_RA / MPI_size);
            int ind = (n - master_N_RA) % (N_RA / MPI_size);
            if (MPI_rank == rank + 1)
            {
                HVCRA_local[ind].writeToFile(filename);
                //printf("Rank = %d; ind = %d\n", MPI_rank, ind);
            }

        }
    }
}

void PoolParallel::write_I(const char* filename, int n)
{
    int master_N_I = N_I/MPI_size + N_I%MPI_size;
    //printf("master_N_RA = %d; n = %d\n", master_N_RA, n);

    if (n >= N_I)
    {
        printf("Selected neuron ID doesn't exist in the pool! Instead wrote neuron 0 to the file.\n");

        if (MPI_rank == 0)
            HVCI_local[0].writeToFile(filename);
    }
    else
    {
        if (n < master_N_I)
        {
            //printf("master_N_RA = %d; n = %d\n", master_N_RA, n);

            if (MPI_rank == 0)
            {
                HVCI_local[n].writeToFile(filename);
                //printf("Master. master_N_I = %d; n = %d\n", master_N_I, n);
            }
        }
        else
        {
            int rank = (n - master_N_I) / (N_I / MPI_size);
            int ind = (n - master_N_I) % (N_I / MPI_size);
            if (MPI_rank == rank + 1)
            {
                HVCI_local[ind].writeToFile(filename);
                //printf("ind = %d\n", ind);
            }

        }
    }
}

void PoolParallel::write_mature(const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;

        out.open(filename, std::ios::out | std::ios::binary );
        out.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));

        for (int i = 0; i < N_RA; i++)
        {
            out.write(reinterpret_cast<char *>(&mature_global[i]), sizeof(mature_global[i]));
	    
	    }
		out.close();
	}

}
