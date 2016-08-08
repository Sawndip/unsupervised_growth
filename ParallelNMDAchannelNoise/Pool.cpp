#include "Pool.h"
#include "HHI_final_pool.h"
#include "HH2_final_pool.h"
#include "poisson_noise.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "training_current.h"
#include <omp.h>

Pool::Pool(int N_tr, int N_ra, int N_i, int Ns) : N_TR(N_tr),
				N_RA(N_ra), N_I(N_i), Nss(Ns)
{
	HVCRA = new HH2_final_pool[N_RA];
	HVCI = new HHI_final_pool[N_I];

	active = new bool*[N_RA];
	active_supersynapses = new bool*[N_RA];
	weights = new double*[N_RA];
	weights_RA_I = new std::vector<double>[N_RA];
	weights_I_RA = new std::vector<double>[N_I];
	syn_ID_RA_I = new std::vector<unsigned>[N_RA];
	syn_ID_I_RA = new std::vector<unsigned>[N_I];
	spike_times = new double[N_RA];
	spikes_in_trial = new std::vector<double>[N_RA];
	remodeled = new bool[N_RA];

	supersynapses = new std::vector<unsigned>[N_RA];
	active_synapses = new std::vector<unsigned>[N_RA];

	for (int i = 0; i < N_RA; i++)
    {
		active[i] = new bool[N_RA];
		active_supersynapses[i] = new bool[N_RA];
		weights[i] = new double[N_RA];

	}

	for (int i = 0; i < N_RA; i++)
	{
		for (int j = 0; j < N_RA; j++)
		{
			active[i][j] = false;
			active_supersynapses[i][j] = false;
			weights[i][j] = 0.0;

		}
		spike_times[i] = 0.0;
		remodeled[i] = false;

	}
}

Pool::~Pool()
{
	delete[] HVCRA;
	delete[] HVCI;

	for (int i = 0; i < N_RA; i++)
	{
		delete[] active[i];
		delete[] active_supersynapses[i];
		delete[] weights[i];
	}
	delete[] active;
	delete[] active_supersynapses;
	delete[] weights;
	delete[] spike_times;
	delete[] remodeled;
}


const double Pool::MIN_INTERNEURON_DISTANCE = 1; // minimum distance between neurons
const double Pool::LAMBDA = 50; // spatial scale of probability of connections decay
const double Pool::CONNECT_CONST = 0.5;
const double Pool::ACTIVATION = 0.01; // activation threshold for synapses
const double Pool::SIDE = 200; // length of HVC side
const double Pool::SUPERSYNAPSE_THRESHOLD = 0.02; // threshold for supersynaptic connection

		// developmental GABA switch
const double Pool::T_GABA = 10000;
const double Pool::E_GABA_IMMATURE = -50;
const double Pool::E_GABA_MATURE = -80;
const int Pool::N_MATURATION = 5;
		// constants for STDP-rules
const double Pool::G_MAX = 0.275; // constant for maximum weight value
const double Pool::BETA = 0.99; // constant for potentiation decay

const double Pool::A_P = 0.01;
const double Pool::G_P = 0.1;
const double Pool::T_P = 5;
const double Pool::TAU_P = 20;

const double Pool::A_D = 0.0105;
const double Pool::T_D = 5.25;
const double Pool::TAU_D = 20;

const double Pool::R = 1;
const double Pool::F_0 = 0.25;



void Pool::set_generator(Poisson_noise* g)
{
	generator = g;
}

void Pool::initialize_pool(double Gie_max, double Gei_max)
{
	double xx; // temporary x-coordinate
	double yy; // temporary y-coordinate

	bool close; // are neurons too close or not

	// set coordinates for HVC(RA) neurons
	for (int i = 0; i < N_RA; i++)
		{
			do
			{
				close = false;
				xx = generator->random(SIDE);
				yy = generator->random(SIDE);
				// check distances to all previous RA neurons
				for (int j = 0; j < i; j++)
				{
					if ((fabs(xx - xx_RA[j]) < MIN_INTERNEURON_DISTANCE)
						&&(fabs(yy - yy_RA[j]) < MIN_INTERNEURON_DISTANCE))
					{
						close = true;
						break;
					}

				}
			} while(close);

			xx_RA.push_back(xx);
			yy_RA.push_back(yy);
			HVCRA[i].set_coordinates(xx, yy);

		}

	// set coordinates for HVC(I) neurons
	for (int i = 0; i < (int) sqrt(N_I); i++)
	{
		for (int k = 0; k < (int) sqrt(N_I); k++)
		{
			do
				{

					close = false;
					xx = (double) i * SIDE / sqrt(N_I) + generator->random(SIDE / (5*sqrt(N_I)));
					yy = (double) k * SIDE / sqrt(N_I) + generator->random(SIDE / (5*sqrt(N_I)));

					// check distances to all RA neurons
					for (int j = 0; j < N_RA; j++)
					{
						if ((fabs(xx - xx_RA[j]) < MIN_INTERNEURON_DISTANCE)
							&&(fabs(yy - yy_RA[j]) < MIN_INTERNEURON_DISTANCE))
						{
							close = true;
							break;
						}

					}
					// check distances to all previous I neurons
					if (!close)
					{
						for (int j = 0; j < i; j++)
						{
							if ((fabs(xx - xx_I[j]) < MIN_INTERNEURON_DISTANCE)
								&&(fabs(yy - yy_I[j]) < MIN_INTERNEURON_DISTANCE))
							{
								close = true;
								break;
							}
						}
					}

				} while(close);
			xx_I.push_back(xx);
			yy_I.push_back(yy);
			HVCI[i+k].set_coordinates(xx, yy);
		}
	}
	// make synaptic connections for HVC(RA) neurons
	for (int i = 0; i < N_RA; i++)
	{
		// all-to-all connections to other HVC(RA) neurons

		// connections to HVC(I) neurons
		for (int j = 0; j < N_I; j++)
		{
			if (generator->random(1) < p(i,j))
			{
				double G = generator->random(Gei_max);
				weights_RA_I[i].push_back(G);
				syn_ID_RA_I[i].push_back(j);
			}
		}
	}

	// make connections for HVC(I) neurons
	for (int i = 0; i < N_I; i++)
	{
		for (int j = 0; j < N_RA; j++)
		{
			if (generator->random(1) < p(j,i))
			{
				double G = generator->random(Gie_max);
				weights_I_RA[i].push_back(G);
				syn_ID_I_RA[i].push_back(j);
			}
		}
	}
}


void Pool::test_initialization_no_connections2I(double Gie_max)
{
    this->initialize_coordinates();

	// make all-to-all connections for HVC(I) neurons
	for (int i = 0; i < N_I; i++)
	{
		for (int j = 0; j < N_RA; j++)
		{
            weights_I_RA[i].push_back(Gie_max);
            syn_ID_I_RA[i].push_back(j);

		}
	}

}

void Pool::test_initialization_no_connections()
{
    this->initialize_coordinates();

}

void Pool::test_initialization(double Gie_max, double Gei_max)
{
    // initialize coordinates

    this->initialize_coordinates();
    
    // make synaptic connections for HVC(RA) neurons
	for (int i = 0; i < N_RA; i++)
	{

		// all-to-all connections to HVC(I) neurons
		for (int j = 0; j < N_I; j++)
		{
            weights_RA_I[i].push_back(Gei_max);
            syn_ID_RA_I[i].push_back(j);

		}
	}

	// make all-to-all connections for HVC(I) neurons
	for (int i = 0; i < N_I; i++)
	{
		for (int j = 0; j < N_RA; j++)
		{
            weights_I_RA[i].push_back(Gie_max);
            syn_ID_I_RA[i].push_back(j);

		}
	}

	// initialize spike times

    /*
	for (int i = 0; i < N_RA; i++)
	{
	    if (i < N_TR)
            spike_times[i] = -10;
        else
            spike_times[i] = 0;
	}
	*/


}

void Pool::initialize_coordinates()
{

	double xx; // temporary x-coordinate
	double yy; // temporary y-coordinate

	bool close; // are neurons too close or not

	// set coordinates for HVC(RA) neurons
	for (int i = 0; i < N_RA; i++)
		{
			do
			{
				close = false;
				xx = generator->random(SIDE);
				yy = generator->random(SIDE);
				// check distances to all previous RA neurons
				for (int j = 0; j < i; j++)
				{
					if ((fabs(xx - xx_RA[j]) < MIN_INTERNEURON_DISTANCE)
						&&(fabs(yy - yy_RA[j]) < MIN_INTERNEURON_DISTANCE))
					{
						close = true;
						break;
					}

				}
			} while(close);

			xx_RA.push_back(xx);
			yy_RA.push_back(yy);
			HVCRA[i].set_coordinates(xx, yy);

		}

	// set coordinates for HVC(I) neurons
	for (int i = 0; i < (int) sqrt(N_I); i++)
	{
		for (int k = 0; k < (int) sqrt(N_I); k++)
		{
			do
				{

					close = false;
					xx = (double) i * SIDE / sqrt(N_I) + generator->random(SIDE / (5*sqrt(N_I)));
					yy = (double) k * SIDE / sqrt(N_I) + generator->random(SIDE / (5*sqrt(N_I)));

					// check distances to all RA neurons
					for (int j = 0; j < N_RA; j++)
					{
						if ((fabs(xx - xx_RA[j]) < MIN_INTERNEURON_DISTANCE)
							&&(fabs(yy - yy_RA[j]) < MIN_INTERNEURON_DISTANCE))
						{
							close = true;
							break;
						}

					}
					// check distances to all previous I neurons
					if (!close)
					{
						for (int j = 0; j < i; j++)
						{
							if ((fabs(xx - xx_I[j]) < MIN_INTERNEURON_DISTANCE)
								&&(fabs(yy - yy_I[j]) < MIN_INTERNEURON_DISTANCE))
							{
								close = true;
								break;
							}
						}
					}

				} while(close);
			xx_I.push_back(xx);
			yy_I.push_back(yy);
			HVCI[i+k].set_coordinates(xx, yy);
		}
	}
}

void Pool::visualize_RA_connections()
{
	for (int i = 0; i < N_RA; i++)
	{
		printf("HVC(RA) neuron %d\n", i);
		printf("Number of connections to I neurons: %d\n", weights_RA_I[i].size());
		for (int j = 0; j < weights_RA_I[i].size(); j++)
		{
			printf("Neuron: %d\tSynaptic conductance: %f\n", syn_ID_RA_I[i][j], weights_RA_I[i][j]);

		}
	}
}

void Pool::visualize_I_connections()
{
	for (int i = 0; i < N_I; i++)
	{
		printf("HVC(I) neuron %d\n", i);
		printf("Connections to HVC(RA) neurons: %d\n", weights_I_RA[i].size());
		for (int j = 0; j < weights_I_RA[i].size(); j++)
		{
			printf("Target neuron %d \t Synaptic conductance %f\n",
					syn_ID_I_RA[i][j], weights_I_RA[i][j]);
		}

	}
}

void Pool::reset_after_trial()
{


	for (int i = 0; i < N_RA; i++)
    {
        spikes_in_trial[i].clear();
		HVCRA[i].reset();
    }
	for (int i = 0; i < N_I; i++)
		HVCI[i].reset();
}


void Pool::set_training_current()
{
    std::function<double (double)> f = &training_current;
	for (int i = 0; i < N_TR; i++)
		HVCRA[i].set_dend_current(f);
}

void Pool::set_generator4neurons()
{
	for (int i = 0; i < N_RA; i++)
		HVCRA[i].set_noise_generator(generator);

	for (int i = 0; i < N_I; i++)
		HVCI[i].set_noise_generator(generator);
}
void Pool::set_dynamics(double interval, double tS)
{
	size = (int) round(interval/tS);
	timeStep = tS;

	for (int i = 0; i < N_RA; i++)
		HVCRA[i].set_dynamics(interval, tS);

	for (int i = 0; i < N_I; i++)
		HVCI[i].set_dynamics(interval, tS);
}

void Pool::ground_state(unsigned N_trial)
{
	for (unsigned i = 0; i < N_trial; i++)
	{
		this->trial();
		this->reset_after_trial();
	}
}

void Pool::set_no_noise()
{
    for (unsigned i = 0; i < N_RA; i++)
        HVCRA[i].set_no_noise();

    for (unsigned i = 0; i < N_I; i++)
        HVCI[i].set_no_noise();
}

void Pool::set_no_noise_RA()
{
    for (unsigned i = 0; i < N_RA; i++)
        HVCRA[i].set_no_noise();

}

void Pool::set_no_noise_I()
{
    for (unsigned i = 0; i < N_I; i++)
        HVCI[i].set_no_noise();

}

void Pool::trial()
{
	bool some_RA_neuron_fired;
	bool some_I_neuron_fired;

	std::vector<unsigned> RA_neurons_fired;
	std::vector<unsigned> I_neurons_fired;

    trial_number++;

	// evolve dynamics
	for (unsigned t = 1; t < size; t++)
	{
		internal_time += timeStep;

		some_RA_neuron_fired = false;
		some_I_neuron_fired = false;

		RA_neurons_fired.clear();
		I_neurons_fired.clear();

		#pragma omp parallel for num_threads(8) shared(some_RA_neuron_fired)
		for (unsigned i = 0; i < N_RA; i++)
		{
		    //printf("My thread number = %d", omp_get_thread_num());
			//printf("Before set Ei\n");
			//HVCRA[i].set_Ei(E_GABA(internal_time));
			int syn_num = active_synapses[i].size();
	

      
            if (supersynapses[i].size() < Nss)
            {
                //HVCRA[i].set_Ei(E_GABA_IMMATURE);
                //if (i < N_TR)
                  //  std::cout << "Training neuron " << i << "with Ei = " << E_GABA(syn_num) << std::endl;
                HVCRA[i].set_Ei(E_GABA(syn_num));
               //printf("syn_num = %d\n", syn_num);
               //printf("E_GABA(syn_num) = %f\n", E_GABA(syn_num));
            }
            else
                HVCRA[i].set_Ei(E_GABA_MATURE);
			//printf("Before RK4 HVCRA step\n");
			HVCRA[i].R4_step_no_target_update();

			// if fired
			//printf("Number of spikes for neuron %d: %d\n", i, HVCRA[i].get_spikes());

			if (HVCRA[i].get_fired())
			{
			    #pragma omp critical
			    {
			        RA_neurons_fired.push_back(i);
			    }

				//printf("Neuron %d fired\n", i);

				//spike_times[i] = internal_time - (trial_number-1)*size *timeStep;
				spike_times[i] = internal_time;
				spikes_in_trial[i].push_back(internal_time);
				//printf("spike_times[%d] = %f\n", i, spike_times[i]);
				if (!some_RA_neuron_fired)
                   			 some_RA_neuron_fired = true;

			}
		}
		#pragma omp parallel for num_threads(8) shared(some_I_neuron_fired)
		for (unsigned i = 0; i < N_I; i++)
		{
			HVCI[i].R4_step_no_target_update();
			//if (some_RA_neuron_fired)
			//	printf("Before RK4 HVCI step\n");

			if (HVCI[i].get_fired())
			{
			    if (!some_I_neuron_fired)
        		            some_I_neuron_fired = true;
				#pragma omp critical
				{
                		    I_neurons_fired.push_back(i);
				}

			}

		}

		// if some HVC(I) neuron fired
		if (some_I_neuron_fired)
		{
			// loop over all fired inhibitory neurons
			for (unsigned i = 0; i < I_neurons_fired.size(); i++)
			{
				unsigned n = I_neurons_fired[i];
				// loop over all targets of fired neuron
				for (unsigned j = 0; j < weights_I_RA[n].size(); j++)
				{
					unsigned k = syn_ID_I_RA[n][j];
					HVCRA[k].raiseI(weights_I_RA[n][j]); // increase inhibitory conductance
					//printf("Here_1\n");
				}
			}
		}

		// if some HVC(RA) neuron fired
		if (some_RA_neuron_fired)
		{
		    std::sort(RA_neurons_fired.begin(), RA_neurons_fired.end());

			// loop over all fired HVC(RA) neurons
			for (unsigned i = 0; i < RA_neurons_fired.size(); i++)
			{
				// update conductances of inhibitory targets
				unsigned n = RA_neurons_fired[i];
				// loop over all HVC(I) targets of fired HVC(RA) neuron
				for (unsigned j = 0; j < weights_RA_I[n].size(); j++)
				{
					unsigned k = syn_ID_RA_I[n][j];
					HVCI[k].raiseE(weights_RA_I[n][j]);
				}

				// update conductances of excitatory targets
				// loop over all targets of fired neurons
				for (unsigned j = 0; j < active_synapses[n].size(); j++)
				{
					unsigned k = active_synapses[n][j];
                    printf("active_synapses[%d][%d] = %d\n", n,j,k);
					HVCRA[k].raiseE(weights[n][k]);
				}

				// apply STDP rules
				// if not saturated
				if (supersynapses[n].size() < Nss)
				{
					STDP(n, i, RA_neurons_fired.begin()); // STDP-rules for all fired neurons
				    //std::cout << "STDP not saturated" << std::endl;
                }
				// if saturated
				else
					STDP_saturated(n);
			}
			//printf("Before update\n");
			//printf("After remodeling\n");
		}
		this->update_synapses(); // update states of synapses
        //printf("After update\n");
        this->axon_remodeling(); // apply axon remodeling

	}
    this->potentiation_decay();
    this->update_synapses();
    this->axon_remodeling();

    for (unsigned jj = 0; jj < N_TR; jj++)
    {
        int ii = active_synapses[jj].size();
        std::cout << "Training neuron " << jj << "with Ei = " << E_GABA(ii) << std::endl;
    }
}

void Pool::potentiation_decay()
{
    for (unsigned i = 0; i < N_RA; i++ )
    {
        for (unsigned j = 0; j < N_RA; j++)
        {
            weights[i][j] *= BETA;
        }
    }
}

double Pool::p(int i_RA, int j_I)
{
	double prob;
	double d;
	d = distance(xx_RA[i_RA], yy_RA[i_RA], xx_I[j_I], yy_I[j_I]);
	if (d < MIN_INTERNEURON_DISTANCE)
		return 0;
	else
	{
		prob = CONNECT_CONST * exp(-d / LAMBDA);
		//printf("p = %f\n", prob);
	}
	return prob;
}

void Pool::LTP(double &w, double t)
{
	if (t <= T_P)
    {
        //if ((w + R * A_P * G_P * ((1 + F_0) * t / T_P - F_0))<0)
          //  printf("LTP. Linear interval. Weight became negative. w = %f\n", w);

		w = w + R * A_P * G_P * ((1 + F_0) * t / T_P - F_0);
		if (w < 0)
            w = 0;
        

        //std::cout << "w = " << w << std::endl;
    }
	else
    {
        //if ((w + R * A_P * G_P * exp(-(t - T_P) / TAU_P))<0)
          //  printf("LTP. Exponential interval. Weight became negative. w = %f\n", w);

    	w = w + R * A_P * G_P * exp(-(t - T_P) / TAU_P);
    	if (w < 0)
            w = 0;


        //std::cout << "w = " << w << std::endl;
    }
}

void Pool::LTD(double &w, double t)
{
	if (t <= T_P)
	{

        //if ((w - R * A_D * w * ((1 - F_0) * t / T_D + F_0))<0)
            //printf("LTD. Linear interval. Weight became negative. w = %f\n", w);
		w = w - R * A_D * w * ((1 - F_0) * t / T_D + F_0);
        if (w < 0)
            w = 0;

       // std::cout << "w = " << w << std::endl;
	}
	else
	{
        //if ((w - R * A_D * w * exp(-(t - T_D) / TAU_D))<0)
          //  printf("LTD. Exponential interval. Weight became negative. w = %f\n", w);

		w = w - R * A_D * w * exp(-(t - T_D) / TAU_D);

        if (w < 0)
            w = 0;

        //std::cout << "w = " << w << std::endl;
	}
}

double Pool::E_GABA(double t)
{
	return E_GABA_MATURE + (E_GABA_IMMATURE - E_GABA_MATURE) * exp(-t/T_GABA);
}

double Pool::E_GABA(int n)
{
    return E_GABA_MATURE + (E_GABA_IMMATURE - E_GABA_MATURE) * exp(-n/(N_MATURATION));
}

void Pool::update_synapses()
{
	for (unsigned i = 0; i < N_RA; i++)
    {
		for (unsigned j = 0; j < N_RA; j++)
		{
			if (i!=j)
			{
				// potentiation decay
				//weights[i][j] *= BETA;
				// check if weights crossed thresholds
				if ((weights[i][j] >= ACTIVATION)&&(!active[i][j])&&(!remodeled[i]))
						{
							active_synapses[i].push_back(j);
							active[i][j] = true;
						}
				// if synapse exceeded sypersynaptic threshold and neuron is not saturated
				if ((weights[i][j] > SUPERSYNAPSE_THRESHOLD)&&(!active_supersynapses[i][j])
					&&(supersynapses[i].size() < Nss))
					{
						supersynapses[i].push_back(j);
						active_supersynapses[i][j] = true;
					}

				// if synapse crossed deactivation threshold
				if ((weights[i][j] < ACTIVATION)&&(active[i][j])&&(!remodeled[i]))
					{
						active[i][j] = false;

						std::vector<unsigned>::iterator pos =
							std::find(active_synapses[i].begin(), active_synapses[i].end(), j);

						if (pos!=active_synapses[i].end())
							active_synapses[i].erase(pos);
						else
							printf("Synapse not found among active!\n");
					}
				
				
					// is synapse decreased below sypersynaptic threshold
				if ((weights[i][j] < SUPERSYNAPSE_THRESHOLD)&&(active_supersynapses[i][j]))
				{
						// if it was saturated and already remodeled, set remodeled to false so
						// that it can be remodeled in the future
					if ((supersynapses[i].size() == Nss)&&(remodeled[i]))
						remodeled[i] = false;

					std::vector<unsigned>::iterator pos =
						std::find(supersynapses[i].begin(), supersynapses[i].end(), j);

					if (pos!=supersynapses[i].end())
					{
						supersynapses[i].erase(pos); // delete synapse from array of supersynapses
						active_supersynapses[i][j] = false;

					}
					else
						printf("Synapse not found among active supersynapses!\n");
				}
				
			}
		}
    }
}

void Pool::STDP(unsigned n, unsigned index_in_fired_array, std::vector<unsigned>::iterator neurons_fired_iter)
{
	bool fired2fired; // indicator that we consider neurons that fired together and we already updated weights for them
    
    //printf("index_in_fired_array = %d\n", index_in_fired_array);
    //for (int i = 0; i < neurons_fired.size(); i++)
    //    printf("neurons_fired[%d] = %d", i, neurons_fired[i]);

	for (int j = 0; j < N_RA; j++)
	{
		fired2fired = false;

		if (j < n)
		{
			for (int k = 0; k < index_in_fired_array; k++)
            {   
                //std::cout << "neurons_fired = " << *(neurons_fired_iter+k) << std::endl;
				if (j == *(neurons_fired_iter+k))
					fired2fired = true;
            }
		}


		// if not the same neuron and if not another fired neuron already taken into account
		if ((j!=n)&&(!fired2fired))
		{
			STDP_one2one(n, j); // modify all connectionos from neuron n
			// if neuron j is not saturated or forms supersynapse on neuron n
			// modify all connections to neuron n
			if ((supersynapses[j].size() < Nss)||(active_supersynapses[j][n]))
				STDP_one2one(j, n);
		}
	}
}

void Pool::STDP_one2one(unsigned i, unsigned j)
{
	double dt;

	// apply LTD or LTP rules
	// if neuron i spikes before neuron j
	if (spike_times[i] <= spike_times[j])
	{
		dt = fabs(spike_times[i] - spike_times[j]);
        //printf("dt = %f\n", dt);
		LTP(weights[i][j], dt);

		if (weights[i][j] > G_MAX)
			weights[i][j] = G_MAX;

        if (weights[i][j] < 0)
            printf("Weight became negative! Weight[%d][%d] = %f\n", i, j, weights[i][j]);

        if (isinf(weights[i][j]))
            std::cout << "weight became infinite :" << i << " " << j << std::endl;
       
        if (isnan(weights[i][j]))
            std::cout << "weight became NaN :" << i << " " << j << std::endl;
	}
	// if neuron i spikes after neuron j
	else
	{
		dt = fabs(spike_times[i] - spike_times[j]);
		//printf("dt = %f\n", dt);
        LTD(weights[i][j], dt);

		if (weights[i][j] < 0)
            printf("Weight became negative! Weight[%d][%d] = %f\n", i, j, weights[i][j]);
        
        if (isinf(weights[i][j]))
            std::cout << "weight became infinite :" << i << " " << j << std::endl;
       
        if (isnan(weights[i][j]))
            std::cout << "weight became NaN :" << i << " " << j << std::endl;
	}
    //std::cout << "weight[" << i << "][" << j < "] = " << weights[i][j] << std::endl;
    //printf("weights[%d][%d] = %f\n", i,j,weights[i][j]);
}

void Pool::STDP_saturated(unsigned n)
{
	double dt;
	// loop among all supersynapses (targets of our saturated neuron n)
	for (int i = 0; i < supersynapses[n].size(); i++)
	{
			int supersynapse_ID = supersynapses[n][i];
            
            printf("Neuron %d has supersynapse onto %d\n", n, supersynapses[n][i]);
			STDP_one2one(n, supersynapse_ID);
	}

	// loop among all other neurons in the pool
	for (int i = 0; i < N_RA; i++)
	{
		// if not the same neuron and if neuron i is not saturated or forms supersynapse
		// with neuron n
		if ((i!=n)&&((supersynapses[i].size() < Nss)||(active_supersynapses[i][n])))
			STDP_one2one(i, n);
	}
}

void Pool::axon_remodeling()
{
	// loop over all neurons
	for (unsigned i = 0; i < N_RA; i++)
	{
		// if neuron became saturated and not already remodeled
		if ((supersynapses[i].size() == Nss)&&(!remodeled[i]))
		{
			// loop over all targets of the saturated neuron
			for (unsigned j = 0; j < active_synapses[i].size(); j++)
			{
				unsigned n = active_synapses[i][j];
				if (n!=i)
				{
					// find active synapse among supersynapses
					std::vector<unsigned>::iterator pos =
						std::find(supersynapses[i].begin(), supersynapses[i].end(), n);

					// if not found, deactivate synapse
					if (pos == supersynapses[i].end())
					{
						active[i][n] = false;
						active_synapses[i].erase(active_synapses[i].begin() + j);
					}


				}

			}
			remodeled[i] = true;
		}
	}
}

void Pool::write_coordinates(const char* xy_RA, const char* xy_I)
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

void Pool::write_invariable_synapses(const char* RA_I, const char* I_RA)
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
		size = syn_ID_RA_I[i-1].size();

		out_RA_I.write(reinterpret_cast<char *>(&i), sizeof(i));
		out_RA_I.write(reinterpret_cast<char *>(&size), sizeof(size)); // write neuron's ID

		for (int j = 0; j < size; j++)
		{
			int k = syn_ID_RA_I[i-1][j];
			double G = weights_RA_I[i-1][j];

			out_RA_I.write(reinterpret_cast<char *>(&k), sizeof(k));
			out_RA_I.write(reinterpret_cast<char *>(&G), sizeof(G));

		}

	}

	// write connections from I to RA
	for (int i = 1; i <= N_I; i++)
	{
		out_I_RA.write(reinterpret_cast<char *>(&i), sizeof(i)); // write neuron's ID number

		size = syn_ID_I_RA[i-1].size();
		out_I_RA.write(reinterpret_cast<char *>(&size), sizeof(size)); // write number of targets a neuron has
		for (int j = 0; j < size; j++)
		{
				int k = syn_ID_I_RA[i-1][j];
				double G = weights_I_RA[i-1][j];
				out_I_RA.write(reinterpret_cast<char *>(&k), sizeof(k)); // write targets ID
				out_I_RA.write(reinterpret_cast<char *>(&G), sizeof(G)); // write targets conductance

		}
	}
	// close files
	out_I_RA.close();
	out_RA_I.close();
}


void Pool::test_connections()
{
    for (int i = 0; i < 5; i++)
    {
        active_synapses[i].push_back((i+1)*2);
        active_synapses[2*i].push_back(i+1);
        weights[i][(i+1)*2] = i;
        weights[2*i][i+1] = i;
    }

}

void Pool::write_variable_synapses(const char * RA_RA)
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

		size = active_synapses[i-1].size();
		out_RA_RA.write(reinterpret_cast<char *>(&size), sizeof(size));

		for (int j = 0; j < size; j++)
		{
			int k = active_synapses[i-1][j];
			out_RA_RA.write(reinterpret_cast<char *>(&k), sizeof(k));
			out_RA_RA.write(reinterpret_cast<char *>(&weights[i-1][k]), sizeof(weights[i-1][k]));

		}

	}
	out_RA_RA.close();
}

void Pool::write_weights(const char * filename)
{
	std::ofstream output;

	output.open(filename, std::ios::out | std::ios::binary);

	output.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));

	for (int i = 0; i < N_RA; i++)
	{
		for (int j = 0; j < N_RA; j++)
		{
			//printf("weigths[%d][%d] = %1.10f\n", i, j, weights[i][j]);
			output.write(reinterpret_cast<char *>(&weights[i][j]), sizeof(weights[i][j]));

		}
	}
	output.close();
}

void Pool::write_time_info(const char* filename)
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
        int spike_array_size = spikes_in_trial[i].size();
        out.write(reinterpret_cast<char *>(&spike_array_size), sizeof(int));
        for (int j = 0; j < spike_array_size; j++)
            out.write(reinterpret_cast<char *>(&spikes_in_trial[i][j]), sizeof(double));
    }
    //out.write(reinterpret_cast<char *>(spike_times), N_RA*sizeof(double));

    out.close();

}

void Pool::write_RA(const char* filename, int n)
{
    if (n < N_RA)
        HVCRA[n].writeToFile(filename);
    else
    {
        printf("Selected neuron ID doesn't exist in the pool! Instead wrote neuron 0 to the file.\n");
        HVCRA[0].writeToFile(filename);
    }

}

void Pool::write_I(const char* filename, int n)
{
    if (n < N_I)
        HVCI[n].writeToFile(filename);
    else
    {

        printf("Selected neuron ID doesn't exist in the pool! Instead wrote neuron 0 to the file.\n");
        HVCI[0].writeToFile(filename);
    }
}

void Pool::write_state(const char* filename)
{
    std::ofstream out;
    
    std::vector<double>::iterator it, it_begin, it_end;
    out.open(filename, std::ios::binary | std::ios::out);

    out.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));
    
    for (unsigned i = 0; i< N_RA ; i++)
    {
        HVCRA[i].get_state(it_begin, it_end);
        std::cout << "state[1] = " <<  *(it_end-1) << std::endl;
        if (std::isinf(*(it_end-1))||(std::isnan(*(it_end-1))))
            std::cout << "Dynamics blew up! Neuron " << i << std::endl;
        else
        {
            for (it = it_begin; it!= it_end; it++)
                out.write(reinterpret_cast<char *>(&(*it)), sizeof(*it));
        }

    }
    out.close();


}
