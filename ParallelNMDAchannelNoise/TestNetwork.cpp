#include "TestNetwork.h"
#include <sys/stat.h>
#include <omp.h>
#include <fstream>

const double TestNetwork::NETWORK_UPDATE_FREQUENCY = 0.1; // synchronization frequency of neural dynamics in network
const double TestNetwork::STDP_WINDOW = 100.0; // burst events outside of this window do not lead to STDP
const double TestNetwork::WAIT_TIME = 150.0; // time to wait before current can be injected to training neurons
const double TestNetwork::TRAINING_KICK = 3.0; // strength of excitatory conductance pulse delivered to training neurons 

TestNetwork::TestNetwork()
{
	N_RA = 0;
	N_TR = 0;
	N_I = 0;
	
	synaptic_params.Nss = 0; // maximum number of allowed outgoing supersynapses per neuron

    synaptic_params.A_P = 0; // constant for LTP
    synaptic_params.A_D = 0; // constant for LTD
    synaptic_params.T_0 = 0; // shift of STDP rule to the right
    synaptic_params.T_P = 0; // best time for LTP response
    synaptic_params.T_D = 0; // best time for LTD response
    synaptic_params.TAU_P = 0; // time scale of LTP efficiency
    synaptic_params.TAU_D = 0; // time scale of LTD efficiency
    
    synaptic_params.BETA = 0; // potentiation decay constant for active and silent synapses
    synaptic_params.BETA_SUPERSYNAPSE = 0; // potentiation decay constant for super synapses

    synaptic_params.ACTIVATION_THRESHOLD = 0; // threshold for synapse activation
    synaptic_params.SUPERSYNAPSE_THRESHOLD = 0; // threshold for supersynapse activation
    synaptic_params.WEIGHT_MAX = 0; // maximum synaptic weight   
    
    gaba_params.E_GABA_MATURE = 0; // GABA reverse potential for mature neurons
    gaba_params.E_GABA_IMMATURE = 0; // GABA reverse potential for newborn neurons

    gaba_params.GABA_RATE_THRESHOLD = 0; // burst rate threshold for maturation. if exceeded, GABA reverse potential decreases by GABA_DOWN
    gaba_params.DEATH_RATE_THRESHOLD = 0; // burst rate threshold for neuronal replacement

    gaba_params.RATE_WINDOW_SHORT = 0; // time window in which "instantaneous" firing rate for a short period is measured
    gaba_params.RATE_WINDOW_LONG = 0; // time window in which "average" firing rate for a long period is measured

    gaba_params.GABA_DOWN = 0; // decrease in GABA potential
    
    // initialize noise generators
    for (int i = 0; i < MAX_NUM_THREADS; i++)
		noise_generators[i].set_seed(1991 + 1000 * i);	
}

void TestNetwork::set_synaptic_parameters(const struct SynapticParameters& syn_params){synaptic_params = syn_params;};
void TestNetwork::set_gaba_parameters(const struct GabaParameters& g_params){gaba_params = g_params;};

void TestNetwork::set_noise_parameters(double mean_s, double sigma_s, double mean_d, double sigma_d)
{
	white_noise_mean_soma = mean_s;
	white_noise_std_soma = sigma_s;
	white_noise_mean_dend = mean_d;
	white_noise_std_dend = sigma_d;
}

void TestNetwork::resize_arrays()
{
	HVCRA.resize(N_RA);
	HVCI.resize(N_I);
	
	// HVC-RA <-> HVC-I connections
	i_syn_ID_RA2I.resize(N_RA);
	i_syn_ID_I2RA.resize(N_I);
	
	d_weights_RA2I.resize(N_RA);
	d_weights_I2RA.resize(N_I);
	
	// HVC-RA -> HVC-RA connections
	i_activeSynapses_id.resize(N_RA);
    i_superSynapses_id.resize(N_RA);

	d_weights.resize(N_RA);
	
	b_active.resize(N_RA);
    b_super.resize(N_RA);

	// arrays with spike times in trial
    d_spikes_in_trial_soma.resize(N_RA);
    d_spikes_in_trial_dend.resize(N_RA);
	d_spikes_in_trial_interneurons.resize(N_I);
	
    // arrays for STDP and maturation processes
    d_last_dend_spike.resize(N_RA);
    d_E_gaba.resize(N_RA);
    b_remodeled.resize(N_RA);
    
    for (int i = 0; i < N_RA; i++)
    {
        b_active[i].resize(N_RA);
        b_super[i].resize(N_RA);
        d_weights[i].resize(N_RA);
		
		std::fill(b_active[i].begin(),  b_active[i].end(), false);
		std::fill(b_super[i].begin(),   b_super[i].end(), false);
        std::fill(d_weights[i].begin(), d_weights[i].end(), 0.0);
    }
    
    std::fill(b_remodeled.begin(), b_remodeled.end(), false);
    
    // set all training neurons mature and remaining pool neurons immature
    std::fill(d_E_gaba.begin() + N_TR, d_E_gaba.end(), gaba_params.E_GABA_IMMATURE);
	std::fill(d_E_gaba.begin(), d_E_gaba.begin() + N_TR, gaba_params.E_GABA_MATURE);

	// number of bursts in recent trials used in firing rate calculations
    i_num_bursts_in_recent_trials.resize(N_RA);

    std::fill(i_num_bursts_in_recent_trials.begin(), i_num_bursts_in_recent_trials.end(), intBuffer(gaba_params.RATE_WINDOW_LONG));

    for (int i = 0; i < N_RA; i++)
        for (int j = 0; j < gaba_params.RATE_WINDOW_LONG; j++)
            i_num_bursts_in_recent_trials[i].push_back(0);
}

void TestNetwork::create_second_layer_growth_sim(int N_tr, int num_groups, int num_neurons_in_group, double Gei, double dGie)
{
	N_RA = num_groups * num_neurons_in_group + N_tr;
	N_TR = N_tr;
	N_I = 1;
	
	// resize arrays
	this->resize_arrays();
	
	// 1st layer has N_TR neurons; 2nd layer has N_RA - N_TR neurons, which are
	// divided into groups that receive different inhibition 
	
	// connect first training neuron to first interneuron
	i_syn_ID_RA2I[0].push_back(0);
	d_weights_RA2I[0].push_back(Gei);
	
	// connect interneuron to HVC-RA neurons in the second layer
	for (int i = 0; i < num_groups; i++)
	{
		for (int j = 0; j < num_neurons_in_group; j++)
		{
			i_syn_ID_I2RA[0].push_back(N_TR + i*num_neurons_in_group + j);
			d_weights_I2RA[0].push_back(dGie * static_cast<double>(i));
		}
	}
	
	// initialize noise generators
	for (int i = 0; i < N_RA; i++)
	{
		HVCRA[i].set_noise_generator(&noise_generators[i / MAX_NUM_THREADS]);
		
		HVCRA[i].set_white_noise(white_noise_mean_soma, white_noise_std_soma,
								 white_noise_mean_dend, white_noise_std_dend);
	
	}	
	
	HVCI[0].set_noise_generator(&noise_generators[0]);
	//HVCI[0].set_poisson_noise();
}

void TestNetwork::simulate(double trial_duration, double timestep, int save_freq, std::string outputDir)
{
	// set dynamics for neurons
	for (int i = 0; i < N_RA; i++)
		HVCRA[i].set_dynamics(timestep);
		
	for (int i = 0; i < N_I; i++)
		HVCI[i].set_dynamics(timestep);
		
    std::string fileTimeWeights = outputDir + "weights_time.bin";
    std::string fileNeuronStates = outputDir + "neuron_states.bin";
   

	// delete files if they exist
	struct stat buf;
	
	if ( stat(fileNeuronStates.c_str(), &buf) == 0 )
		std::remove(fileNeuronStates.c_str());

    if ( stat(fileTimeWeights.c_str(), &buf) == 0 )
		std::remove(fileTimeWeights.c_str());

    
    int trial_number = 1;
    
    while (true)
    {
        std::cout << "Trial " << trial_number << std::endl;

        this->trial(trial_duration, timestep);
		
		//std::cout << "trial is done" << std::endl;
		
        if (trial_number % save_freq == 0)
        {
            this->write_time_weights(trial_number, fileTimeWeights.c_str());
            this->write_neuron_states(trial_number, fileNeuronStates.c_str());
            
        }
        
        this->reset_after_trial(trial_duration);
        trial_number++;
        
    }
}

void TestNetwork::trial(double trial_duration, double timestep)
{
	double d_internal_time = 0;
    double d_network_time = NETWORK_UPDATE_FREQUENCY;
    
    std::vector<int> i_RA_fired_soma; // RA neurons that fired somatic spikes
    std::vector<int> i_RA_fired_dend; // RA neurons that fired dendritic spikes
    std::vector<int> i_I_fired; // fired interneurons
    
    bool b_training_excited = false; // indicator that training neurons were excited in this trial
 
	double d_training_kick_time = noise_generators[0].random(trial_duration - 2*WAIT_TIME) + WAIT_TIME;
	int i_training_kick_time = static_cast<int>(d_training_kick_time / timestep);
 
    // evolve dynamics
    for (int t = 0; t < static_cast<int>(trial_duration / timestep); t++)
	{
		d_internal_time += timestep;
		
		if ( ( !b_training_excited ) && ( t == i_training_kick_time ) )
		{
			for (int i = 0; i < N_TR; i++)
				HVCRA[i].raiseE(TRAINING_KICK);
				
			b_training_excited = true;
		}
		
		#pragma omp parallel for
		for (int i = 0; i < N_RA; i++)
		{
			//std::cout << "Number of threads = " << omp_get_num_threads() << std::endl;
			
		    HVCRA[i].set_Ei(d_E_gaba[i]);
		    
		    HVCRA[i].set_noise_generator(&noise_generators[omp_get_thread_num()]);
            
            // Debraband step
            HVCRA[i].Debraband_step_no_target_update();
            
            // if some neuron produced somatic spike
            if ( HVCRA[i].get_fired_soma() )
            {
                d_spikes_in_trial_soma[i].push_back(d_internal_time);
                
                # pragma omp critical
                {
					i_RA_fired_soma.push_back(i);
				}


                //for (int j = 0; j < last_soma_spikes_local[i].size(); j++)
                //{
                  //  printf("My rank = %d; Soma spikes of neuron %d:  spike_time = %f\n", MPI_rank, Id_RA_local[i],
                    //    last_soma_spikes_local[i][j]);

                //}

                //printf("My rank = %d; Soma neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], spike_times_soma_local[i]);

                //RA_neurons_fired_soma_realID.push_back(Id_RA_local[i]);
           
            } // end of if_get_fired_soma looop
            
            // if some neuron produced dendritic spike, store this neuron in array
            if (HVCRA[i].get_fired_dend())
            {
                // apply LTD
                // if neuron is saturated apply LTD only to supersynapses
                if ( static_cast<int>(i_superSynapses_id[i].size()) == synaptic_params.Nss )
                {
                    for (size_t k = 0; k < i_superSynapses_id[i].size(); k++)
                    {
                        int target_ID = i_superSynapses_id[i][k];
	
						//std::cout << "Rank = " << MPI_rank << " target_ID = " << target_ID << std::endl;
						double dt = d_internal_time - d_last_dend_spike[target_ID];

						//std::cout << "dt = " << dt << std::endl;

						if ( ( i != target_ID ) && ( dt < STDP_WINDOW ) )
						{
							
							  //      printf("Rank = %d; LTD from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f\n", MPI_rank, i, target_ID, internal_time,
								//			last_dend_spike[target_ID]);
							LTD_burst(d_weights[i][target_ID], -dt);
							update_synapse(i, target_ID);
						}
					}
				
                }
                // if not saturated apply LTD rule with the last dendritic spike of all neurons
                else
                {
                    for (int j = 0; j < N_RA; j++)
                    {
						if ( i != j )
                        {
							double dt = d_internal_time - d_last_dend_spike[j];

							//std::cout << "dt = " << dt << std::endl;
                        
							if ( dt < STDP_WINDOW )
							{
                            //    printf("Rank = %d; LTD from neuron %d onto %d; somatic spike at %f; dendritic spike at %f\n", MPI_rank, i, j, internal_time,
							//			last_dend_spike[j]);
								double old_w = d_weights[i][j];
							
                           		LTD_burst(d_weights[i][j], -dt);
                           		update_synapse(i, j);
                           		
                           		std::cout << "LTD " << i << " -> " << j << " (" << d_internal_time << ", " 
												<< d_last_dend_spike[j] << ") dw = " << d_weights[i][j] - old_w 
												<< " w = " << d_weights[i][j] << "\n";
										
							}
						}
                    }
                }
                
                # pragma omp critical
                {
					i_RA_fired_dend.push_back(i);
                }
                
                d_spikes_in_trial_dend[i].push_back(d_internal_time);
                
                //printf("My rank = %d; Dend neuron %d fired; spike_time = %f\n", MPI_rank, i, spikes_in_trial_dend[i].back());
            }
        }
		
		#pragma omp parallel for
        for (int i = 0; i < N_I; i++)
		{
			HVCI[i].set_noise_generator(&noise_generators[omp_get_thread_num()]);
            HVCI[i].DP8_step_no_target_update();
            
            if ( HVCI[i].get_fired() )
            {
				# pragma omp critical
                {
					i_I_fired.push_back(i);
				}
				
				d_spikes_in_trial_interneurons[i].push_back(d_internal_time);
            }
	    }
        
        if ( d_internal_time >= d_network_time )
        {
            // update all targets of HVC(RA) neurons
            for (size_t i = 0; i < i_RA_fired_soma.size(); i++)
            {
                // update interneurons
                int neuron_id = i_RA_fired_soma[i];
                
                size_t num_targets = i_syn_ID_RA2I[neuron_id].size();
                
                for (size_t j = 0; j < num_targets; j++)
					HVCI[i_syn_ID_RA2I[neuron_id][j]].raiseE(d_weights_RA2I[neuron_id][j]);
                
                // update all excitatory targets
                size_t num_active_synapses = i_activeSynapses_id[neuron_id].size();
                
                for (size_t j = 0; j < num_active_synapses; j++)
                {
                    int target_ID = i_activeSynapses_id[neuron_id][j];
                    HVCRA[target_ID].raiseE(d_weights[neuron_id][target_ID]);
                }
            }

            // update all targets of HVC(I) neurons
            for (size_t i = 0; i < i_I_fired.size(); i++)
            { 
				int neuron_id = i_I_fired[i];
                
                size_t num_targets = i_syn_ID_I2RA[neuron_id].size();
                
                for (size_t j = 0; j < num_targets; j++)
					HVCRA[i_syn_ID_I2RA[neuron_id][j]].raiseI(d_weights_I2RA[neuron_id][j]);
			}
			
			// for all neurons that had a dendritic spike
			// update last dendritic spike time
            for (size_t i = 0; i < i_RA_fired_dend.size(); i++)
				d_last_dend_spike[i_RA_fired_dend[i]] = d_internal_time;
			
            for (size_t i = 0; i < i_RA_fired_dend.size(); i++)
            {
                // apply LTP rule for recent dendritic bursts of other neurons and dendritic spike of fired neurons
                for (int j = 0; j < N_RA; j++)
                {
                    // if neuron is saturatedapply LTP only if dendritic spike occured in supersynapse
                    if ( static_cast<int>(i_superSynapses_id[j].size()) == synaptic_params.Nss  )
                    {
                        std::vector<int>::iterator pos = std::find(i_superSynapses_id[j].begin(),
                                        i_superSynapses_id[j].end(), i_RA_fired_dend[i]);

                        if ( pos != i_superSynapses_id[j].end() )
                        {
							//std::cout << "Rank = " << MPI_rank << " RA_fired_dend = " << RA_fired_dend[i] << std::endl;
							double dt = d_last_dend_spike[i_RA_fired_dend[i]] - d_last_dend_spike[j];

							if ( dt < STDP_WINDOW )
							{
								// if time is smaller than T_0 apply LTD, otherwise apply LTP
								if ( dt <= synaptic_params.T_0 )
								{
									LTD_burst(d_weights[j][i_RA_fired_dend[i]], dt);
									update_synapse(j, i_RA_fired_dend[i]);
								}
								else
								{
									//printf("Rank = %d; LTP from neuron %d onto %d; somatic spike at %f; dendritic spike at %f\n", MPI_rank, j, RA_fired_dend[i],
									//                      spikes_in_trial_soma[j][k], last_dend_spike[RA_fired_dend[i]]);
									LTP_burst(d_weights[j][i_RA_fired_dend[i]], dt);
									update_synapse(j, i_RA_fired_dend[i]);
								}

							}
						}
                        
                    }

                    // if not saturated apply LTP for all dendritic spikes
                    else
                    {
                        // don't allow self-to-self connections
                        if ( j != i_RA_fired_dend[i] )
                        {
							//std::cout << "Rank = " << MPI_rank << " RA_fired_dend = " << RA_fired_dend[i] << std::endl;
							double dt = d_last_dend_spike[i_RA_fired_dend[i]] - d_last_dend_spike[j];

							if ( dt < STDP_WINDOW )
							{
								// if time is smaller than T_0 apply LTD, otherwise apply LTP
								if ( dt <= synaptic_params.T_0 )
								{
									LTD_burst(d_weights[j][i_RA_fired_dend[i]], dt);
									update_synapse(j, i_RA_fired_dend[i]);
								}
								else
								{
									
									
									
									double old_w = d_weights[j][i_RA_fired_dend[i]];
									
									LTP_burst(d_weights[j][i_RA_fired_dend[i]], dt);
									
									std::cout << "LTP " << j << " -> " << i_RA_fired_dend[i] << " (" << d_last_dend_spike[j] << ", " 
												<< d_last_dend_spike[i_RA_fired_dend[i]] << ") dw = " << d_weights[j][i_RA_fired_dend[i]] - old_w 
												<< " w = " << d_weights[j][i_RA_fired_dend[i]] << "\n";
										
									update_synapse(j, i_RA_fired_dend[i]);
								}
							}						
                        }
                    }   
                }
            } // end of RA_fired_dend loop

            i_RA_fired_soma.clear();
            i_RA_fired_dend.clear();
            i_I_fired.clear();
        
            d_network_time += NETWORK_UPDATE_FREQUENCY;
            // check if we need axon remodeling

            for (int i = 0; i < N_RA; i++)
            {
                if ( ( static_cast<int>(i_superSynapses_id[i].size()) == synaptic_params.Nss ) && ( !b_remodeled[i] ) )
                {
                    std::cout << "Axon remodeling for neuron " << i << std::endl;
                    this->axon_remodeling(i);
                }
            }
		} // end of if internal_time > network_time

    } // end of time loop
    this->potentiation_decay();
    
    //printf("After potentiation decay")
    
    this->update_all_synapses();

    for (int i = 0; i < N_RA; i++)
        i_num_bursts_in_recent_trials[i].push_front(static_cast<int>(d_spikes_in_trial_dend[i].size()));

    this->update_Ei();
}

void TestNetwork::reset_after_trial(double trial_duration)
{
    for (int i = 0; i < N_RA; i++)
    {
	    HVCRA[i].reset_time();

        d_spikes_in_trial_soma[i].clear();
        d_spikes_in_trial_dend[i].clear();
        
        
    }

	std::for_each(d_last_dend_spike.begin(), d_last_dend_spike.end(), [=](double &t){ t -= trial_duration; });

    for (int i = 0; i < N_I; i++)
    {
        HVCI[i].reset_time();
		d_spikes_in_trial_interneurons[i].clear();
	}
}

void TestNetwork::LTP_burst(double &w, double dt)
{
	// dt in LTP_burst is postsynaptic burst time - presynaptic burst time  
	if (dt <= synaptic_params.T_0)
	{
		std::cerr << "Time t = " << dt << " is smaller than T_0 in LTP!" << std::endl;
		return;
	}
	
	if (dt <= synaptic_params.T_0 + synaptic_params.T_P)
		w = w + synaptic_params.A_P * (dt - synaptic_params.T_0) / synaptic_params.T_P;
	else
		w = w + synaptic_params.A_P * exp(-(dt - synaptic_params.T_0 - synaptic_params.T_P) / synaptic_params.TAU_P);
		
	if (w > synaptic_params.WEIGHT_MAX)
        w = synaptic_params.WEIGHT_MAX;

    if (w < 0)
        w = 0;
}

void TestNetwork::LTD_burst(double &w, double dt)
{
	// dt in LTD_burst is postsynaptic burst time - presynaptic burst time  
    if (dt > synaptic_params.T_0)
    {
		std::cerr << "Time in LTD dt = " << dt << " is bigger than T0 = " << synaptic_params.T_0 << std::endl;
		return;
	}
	
	if (dt >= synaptic_params.T_0 - synaptic_params.T_D)
		w = w + w * synaptic_params.A_D * (dt - synaptic_params.T_0) / synaptic_params.T_D;
	else
		w = w - w * synaptic_params.A_D * exp((dt - synaptic_params.T_0 + synaptic_params.T_D) / synaptic_params.TAU_D);
		
	if (w < 0)
		w = 0;
}


void TestNetwork::update_Ei()
{
	for (int i = N_TR; i < N_RA; i++)
	{
	// if not a training neuron and not mature but rate exceeds threshold
        double average_burst_rate_short_term = std::accumulate(i_num_bursts_in_recent_trials[i].begin(), i_num_bursts_in_recent_trials[i].begin() + gaba_params.RATE_WINDOW_SHORT, 0.0) 
                                                                                                                / static_cast<double>(gaba_params.RATE_WINDOW_SHORT);
        
        double average_burst_rate_long_term = std::accumulate(i_num_bursts_in_recent_trials[i].begin(), i_num_bursts_in_recent_trials[i].end(), 0.0) 
                                                                                                                / static_cast<double>(gaba_params.RATE_WINDOW_LONG);
        
        //else if ( (trial_number >= RATE_WINDOW_LONG) && (average_burst_rate_long_term <= SILENT_RATE_THRESHOLD) )
        //{
        //    E_gaba[i] += GABA_UP;

         //   if (E_gaba[i] >= E_GABA_IMMATURE)
         //       E_gaba[i] = E_GABA_IMMATURE;
        //}

		if ( average_burst_rate_short_term >= gaba_params.GABA_RATE_THRESHOLD )
		{
			d_E_gaba[i] = d_E_gaba[i] - gaba_params.GABA_DOWN * static_cast<int>(d_spikes_in_trial_dend[i].size()); // decrease gaba potential
					
			if ( d_E_gaba[i] <= gaba_params.E_GABA_MATURE ) // if gaba potential hit the bottom, make neuron mature
				d_E_gaba[i] = gaba_params.E_GABA_MATURE;
			
		}
	}
}

void TestNetwork::potentiation_decay()
{
    for (int i = 0; i < N_RA; i++)
    {
        for (int j = 0; j < N_RA; j++)
        {
            if (d_weights[i][j] >= synaptic_params.SUPERSYNAPSE_THRESHOLD)
                d_weights[i][j] -= synaptic_params.BETA_SUPERSYNAPSE;
                
            else
                d_weights[i][j] -= synaptic_params.BETA;
                
			if (d_weights[i][j] < 0)
				d_weights[i][j] = 0;
        }
    }
}

void TestNetwork::axon_remodeling(int i)
{
	// erase all active synapses if they are not among supersynapses

	std::vector<int> active_shifts_to_erase; // active synapses to erase

	// find all active synapses needed to be erased
    for (size_t j = 0; j < i_activeSynapses_id[i].size(); j++)
    {
        int syn_ID = i_activeSynapses_id[i][j];

        // check if active synapse is among supersynapses
        std::vector<int>::iterator pos = std::find(i_superSynapses_id[i].begin(),
                                i_superSynapses_id[i].end(), syn_ID);
        // if not erase it
        if ( pos == i_superSynapses_id[i].end() )
        {
            active_shifts_to_erase.push_back(j); // save shift of synapse to erase
            b_active[i][syn_ID] = false;
        }
    }

	// sort array with shifts
	std::sort(active_shifts_to_erase.begin(), active_shifts_to_erase.end());

	// erase synapses
	for (size_t j = 0; j < active_shifts_to_erase.size(); j++)
	{
		i_activeSynapses_id[i][active_shifts_to_erase[j]] = i_activeSynapses_id[i].back();
		i_activeSynapses_id[i].pop_back();
	}

    b_remodeled[i] = true;
}

void TestNetwork::update_all_synapses()
{
    //printf("Process %d; Updating all synapses after potentiation decay\n", MPI_rank);
    for (int i = 0; i < N_RA; i++)
    {
        for (int j = 0; j < N_RA; j++)
        {
            if (i != j)
                this->update_synapse(i, j);	
        }
    }
}

void TestNetwork::update_synapse(int i, int j)
{
	double w = d_weights[i][j];

    if ( (w >= synaptic_params.ACTIVATION_THRESHOLD) && ( !b_active[i][j] ) && ( !b_remodeled[i] ) )
    {
        //std::cout << "Rank = " << MPI_rank << " Activated synapse " << i << " -> " << j << " weight = " << w << std::endl;;
        b_active[i][j] = true;
        i_activeSynapses_id[i].push_back(j);
    }

    if ( ( w >= synaptic_params.SUPERSYNAPSE_THRESHOLD ) && ( !b_super[i][j] ) && ( static_cast<int>(i_superSynapses_id[i].size()) < synaptic_params.Nss ) )
    {
        //std::cout << "Rank = " << MPI_rank << " Activated supersynapse " << i << " -> " << j << " weight = " << w << std::endl;;
        b_super[i][j] = true;
        i_superSynapses_id[i].push_back(j);
    }

    if ( ( w < synaptic_params.SUPERSYNAPSE_THRESHOLD ) && ( b_super[i][j] ) )
    {
        //std::cout << "Rank = " << MPI_rank << " Deactivated supersynapse " << i << " -> " << j << " weight = " << w << std::endl;;
        b_super[i][j] = false;
        b_remodeled[i] = false;
        std::vector<int>::iterator pos = std::find(i_superSynapses_id[i].begin(),
                                                    i_superSynapses_id[i].end(), j);

        if ( pos != i_superSynapses_id[i].end() )
            i_superSynapses_id[i].erase(pos);
        else
            std::cout << "Supersynapse " << i << " -> " << j << " to be erased is not found!" << std::endl;
    }

    if ( ( w < synaptic_params.ACTIVATION_THRESHOLD ) && ( b_active[i][j] ) )
    {
        //std::cout << "Rank = " << MPI_rank << " Deactivated synapse " << i << " -> " << j << " weight = " << w << std::endl;;
        b_active[i][j] = false;

        std::vector<int>::iterator pos = std::find(i_activeSynapses_id[i].begin(),
                                                    i_activeSynapses_id[i].end(), j);

        if ( pos != i_activeSynapses_id[i].end() )
            i_activeSynapses_id[i].erase(pos);
        else
            std::cout << "Active synapse " << i << " -> " << j << " to be erased is not found!" << std::endl;
    }

}




void TestNetwork::write_neuron_states(int trial_number, const char* filename)
{
    std::ofstream out;

    // if file doesn't exist, create new file and write some parameters
    // otherwise append to the end
    
    struct stat buf;
	
	if ( stat(filename, &buf) == 0 )
	   out.open(filename, std::ios::binary | std::ios::out | std::ios::app);
	else
    {
        out.open(filename, std::ios::binary | std::ios::out);
        out.write(reinterpret_cast<char*>(&N_RA), sizeof(N_RA)); // write current trial number
    }
     

    out.write(reinterpret_cast<char*>(&trial_number), sizeof(trial_number)); // write current trial number
    
    // calculate average firing rate
    std::vector<double> mean_firing_rate(N_RA);

    for (int i = 0; i < N_RA; i++)
        mean_firing_rate[i] = std::accumulate(i_num_bursts_in_recent_trials[i].begin(), i_num_bursts_in_recent_trials[i].begin() + gaba_params.RATE_WINDOW_SHORT, 0.0) 
																/ static_cast<double>(gaba_params.RATE_WINDOW_SHORT);

    for (int i = 0; i < N_RA; i++)
    {
		int i_remodeled = static_cast<int>(b_remodeled[i]);
		
        out.write(reinterpret_cast<char*>(&i_remodeled), sizeof(int)); 
        out.write(reinterpret_cast<char*>(&d_E_gaba[i]), sizeof(double)); 
        out.write(reinterpret_cast<char*>(&mean_firing_rate[i]), sizeof(double)); 
	}
	
    out.close();
}


void TestNetwork::write_time_weights(int trial_number, const char* filename)
{
	std::ofstream out;

	// if file doesn't exist, create new file and write some parameters
    // otherwise append to the end
    
    struct stat buf;
	
	if ( stat(filename, &buf) == 0 )
	   out.open(filename, std::ios::binary | std::ios::out | std::ios::app);
	else
    {
        out.open(filename, std::ios::binary | std::ios::out);
        out.write(reinterpret_cast<char*>(&N_RA), sizeof(N_RA)); // write current trial number
    }

    out.write(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));
    
    for (int i = 0; i < N_RA; i++)
        for (int j = 0; j < N_RA; j++)
            out.write(reinterpret_cast<char *>(&d_weights[i][j]), sizeof(d_weights[i][j]));

    
    out.close();
}
