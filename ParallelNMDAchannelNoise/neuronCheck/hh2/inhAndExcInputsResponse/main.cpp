#include "HHI_buffer.h" 
#include "HH2_test.h"
#include <iostream>
#include "poisson_noise.h"
#include <vector>
#include <algorithm>
#include <sys/stat.h>
#include <string>

int main()
{
	double Gie = 0.0; // strength of inhibitory kick
	
	double Gei = 0.25;
	
	double G_TRAINING_KICK = 3.0; // strength of excitatory kick delivered to training neurons
	double T_KICK = 100; // time when training neurons receive excitatory kick
	
	double std_delays = 0; // standard deviation of axonal delays
	double mean_delays = 0; // mean axonal delay

	int num_source_neurons = 10; // number of neurons in first layer (sources)
	int num_target_neurons = 200; // number of neurons in second layer (targets)

	// white noise parameters
	double white_noise_mean_soma = 0.000000;
	double white_noise_std_soma = 0.10000; // was 0.10
	double white_noise_mean_dend = 0.000000;
	double white_noise_std_dend = 0.20000; // was 0.10
	
	// trial parameters
	double TIMESTEP = 0.02;
	double TRIAL_DURATION = 200;
	
	std::vector<std::vector<std::pair<double,int>>> delivery_queue_RA_RA; // queue with spike delivery times for HVC-RA -> HVC-RA interactions
		
	double E_GABA_MATURE = -80.000000;
	double E_GABA_IMMATURE = -80.000000;

	double E_REST_MATURE = -80.000000;
	double E_REST_IMMATURE = -55.000000;

	double AD_MATURE = 10000.000000;
	double AD_IMMATURE = 10000.000000;

	double GK_MATURE = 8.000000;
	double GK_IMMATURE = 8.000000;

	double GNA_MATURE = 60.000000;
	double GNA_IMMATURE = 60.000000;

	double RC_MATURE = 55.000000;
	double RC_IMMATURE = 55.0000;
	
	double GCA_MATURE = 55.0;
	double GCA_IMMATURE = 0.0;
	
	double GCAK_MATURE = 150.0;
	double GCAK_IMMATURE = 150.0;
	
	double GSL_MATURE = 0.1;
	double GSL_IMMATURE = 0.1;
	
	double GDL_MATURE = 0.1;
	double GDL_IMMATURE = 0.1;

	double age = 100;
	
	double gaba_potential = E_GABA_MATURE + (E_GABA_IMMATURE - E_GABA_MATURE) * exp(-age);
	double rest_potential = E_REST_MATURE + (E_REST_IMMATURE - E_REST_MATURE) * exp(-age);
	
	double Gk = GK_MATURE + (GK_IMMATURE - GK_MATURE) * exp(-age);
	double GNa = GNA_MATURE + (GNA_IMMATURE - GNA_MATURE) * exp(-age);
	double GCa = GCA_MATURE + (GCA_IMMATURE - GCA_MATURE) * exp(-age);
	double GCaK = GCAK_MATURE + (GCAK_IMMATURE - GCAK_MATURE) * exp(-age);
	double GsL = GSL_MATURE + (GSL_IMMATURE - GSL_MATURE) * exp(-age);
	double GdL = GDL_MATURE + (GDL_IMMATURE - GDL_MATURE) * exp(-age);
	
	double Ad = AD_MATURE + (AD_IMMATURE - AD_MATURE) * exp(-age);
	double Rc = RC_MATURE + (RC_IMMATURE - RC_MATURE) * exp(-age);
	
	std::cout << "E_GABA = " << gaba_potential << "\n"
			  << "E_REST = " << rest_potential << "\n"
			  << "G_K = "    << Gk             << "\n"
			  << "G_NA = "   << GNa            << "\n"
			  << "G_CA = "   << GCa            << "\n"
			  << "G_CAK = "  << GCaK           << "\n"
			  << "G_SL = "   << GsL            << "\n"
			  << "G_DL = "   << GdL            << "\n"
			  << "Ad = "     << Ad             << "\n"
			  << "Rc = "     << Rc             << std::endl;
			  
	
	std::vector<HH2_test> source_neurons(num_source_neurons);
	std::vector<HH2_test> target_neurons(num_target_neurons);
	HHI_buffer interneuron;
	
	unsigned seed = 1991;
	Poisson_noise noise_generator;
	
	noise_generator.set_seed(seed);
	
	// set parameters for source neurons
	for (int i = 0; i < num_source_neurons; i++)
	{
		source_neurons[i].set_dynamics(TIMESTEP);
		source_neurons[i].set_noise_generator(&noise_generator);
		source_neurons[i].set_white_noise(white_noise_mean_soma, white_noise_std_soma, 
									   white_noise_mean_dend, white_noise_std_dend);
		
		source_neurons[i].set_Ei(E_GABA_MATURE);
		source_neurons[i].set_Erest(E_REST_MATURE);
		source_neurons[i].set_Gk(GK_MATURE);
		source_neurons[i].set_GNa(GNA_MATURE);
		source_neurons[i].set_GCa(GCA_MATURE);
		source_neurons[i].set_GCaK(GCAK_MATURE);
		source_neurons[i].set_GsL(GSL_MATURE);
		source_neurons[i].set_GdL(GDL_MATURE);
		source_neurons[i].set_Ad(AD_MATURE);
		source_neurons[i].set_Rc(RC_MATURE);
	}
	
	;
	
	// set neuron properties in 2nd layer
	for (int i = 0; i < num_target_neurons; i++)
	{
		target_neurons[i].set_dynamics(TIMESTEP);
		
		target_neurons[i].set_Ei(gaba_potential);
		target_neurons[i].set_Erest(rest_potential);
		target_neurons[i].set_Gk(Gk);
		target_neurons[i].set_GNa(GNa);
		target_neurons[i].set_GCa(GCa);
		target_neurons[i].set_GCaK(GCaK);
		target_neurons[i].set_GsL(GsL);
		target_neurons[i].set_GdL(GdL);
		target_neurons[i].set_Ad(Ad);
		target_neurons[i].set_Rc(Rc);
		
		target_neurons[i].set_noise_generator(&noise_generator);
		target_neurons[i].set_white_noise(white_noise_mean_soma, white_noise_std_soma, 
									   white_noise_mean_dend, white_noise_std_dend);
	
	}
	
	// set parameters for interneuron
	interneuron.set_dynamics(TIMESTEP);
	
	double weight_ee = 0; // strenght of excitatory weight from first to second layer neurons
	double dweight_ee = 10;
	int num_steps = 15; 
	
	std::vector<double> weights(num_steps);
	
	std::vector<int> total_num_somatic_spikes_in_trials(num_steps);
	std::vector<int> total_num_dendritic_spikes_in_trials(num_steps);
	
	std::vector<int> num_soma_spiked_in_trials(num_steps);
	std::vector<int> num_dend_spiked_in_trials(num_steps);
	
	std::vector<double> first_soma_spike_delays_in_trials(num_steps);
	std::vector<double> dend_spike_delays_in_trials(num_steps);
	
	std::vector<double> mean_first_soma_spike_delays_in_trials(num_steps);
	std::vector<double> std_first_soma_spike_delays_in_trials(num_steps);
	std::vector<double> mean_dend_spike_delays_in_trials(num_steps);
	std::vector<double> std_dend_spike_delays_in_trials(num_steps);
	
	for (int k = 0; k < num_steps; k++)
	{
		weights[k] = weight_ee;
		
		std::string filename = "/home/eugene/Output/neuronTest/inhAndExcInputsResponse/RA" + std::to_string(k) + ".bin";
	
		struct stat buf;
	
		if ( stat(filename.c_str(), &buf) == 0 )
			std::remove(filename.c_str());
	
		target_neurons[0].set_recording_full(filename.c_str());
		
		bool training_excited = false;
		bool training_fired = false;
		
		double training_first_spike_time = T_KICK;
		
		std::vector<int> num_somatic_spikes(num_target_neurons);
		std::vector<int> num_dendritic_spikes(num_target_neurons);
		
		std::vector<double> first_soma_spike_delays;
		std::vector<double> dend_spike_delays;
		
		std::vector<bool> b_spike_indicators(num_target_neurons);
		std::fill(b_spike_indicators.begin(), b_spike_indicators.end(), false);
		
		// trial
		for (int t = 0; t < static_cast<int>(TRIAL_DURATION / TIMESTEP); t++)
		{
			if ( ( static_cast<double>(t) * TIMESTEP > T_KICK ) && ( !training_excited ) )
			{
				for (int i = 0; i < num_source_neurons; i++)
					source_neurons[i].raiseE(G_TRAINING_KICK);
				training_excited = true;
			}
				
			
			for (int i = 0; i < num_source_neurons; i++)
			{
				source_neurons[i].Debraband_step_no_target_update();
				
				
				if ( source_neurons[i].get_fired_soma() )
				{
					if ( !training_fired ){
						training_first_spike_time = static_cast<double>(t) * TIMESTEP;
						training_fired = true;
					}
					
					//std::cout << "source neuron " << i << " somatic spike at " << static_cast<double>(t) * TIMESTEP << std::endl;
					
					for ( int j = 0; j < num_target_neurons; j++ )
						target_neurons[j].raiseExcWeight(weight_ee);
						
					if ( i == 0 )
						interneuron.raiseE(Gei);
				}
			}
			
			for (int i = 0; i < num_target_neurons; i++)
			{
				target_neurons[i].Debraband_step_no_target_update();
				
				if ( target_neurons[i].get_fired_soma() )
				{
					if ( ( !b_spike_indicators[i] ) && (training_fired) )
					{
						//if ( training_fired )
							first_soma_spike_delays.push_back(static_cast<double>(t) * TIMESTEP - training_first_spike_time);
						//else{
							//std::cout << "Training neuron didn't fire! Time = " << static_cast<double>(t) * TIMESTEP << "\n" << std::endl;
							//first_soma_spike_delays.push_back(static_cast<double>(t) * TIMESTEP - T_KICK);
						//}
						
						b_spike_indicators[i] = true;
					}
					
					num_somatic_spikes[i] += 1;
					
					//std::cout << "target neuron " << i << " somatic spike at " << static_cast<double>(t) * TIMESTEP << std::endl;
					
				}
				if ( target_neurons[i].get_fired_dend() )
				{
					if ( training_fired )
						dend_spike_delays.push_back(static_cast<double>(t) * TIMESTEP - training_first_spike_time);
					else{
						std::cout << "Training neuron didn't fire! Time = " << static_cast<double>(t) * TIMESTEP << "\n" << std::endl;
						dend_spike_delays.push_back(static_cast<double>(t) * TIMESTEP - T_KICK);
					}
					num_dendritic_spikes[i] += 1;
				}	
			}
			
			interneuron.DP8_step_no_target_update();
			
			if ( interneuron.get_fired() )
			{
				for (int i = 0; i < num_target_neurons; i++)
					target_neurons[i].raiseI(Gie);
					
				std::cout << "interneuron spike at " << static_cast<double>(t) * TIMESTEP << std::endl;
			}
		}
		
		int num_spiked_neurons = static_cast<int>(first_soma_spike_delays.size());
		
		if ( num_spiked_neurons > 0 )
		{
			double mean = std::accumulate(first_soma_spike_delays.begin(), first_soma_spike_delays.end(), 0.0) / static_cast<double>(num_spiked_neurons);
			mean_first_soma_spike_delays_in_trials[k] = mean;
			
			if ( num_spiked_neurons > 1 )
			{
				double std = 0;
				std::for_each(first_soma_spike_delays.begin(), first_soma_spike_delays.end(),
									[&](double &x){std += (x - mean)*(x - mean);});
									
				std_first_soma_spike_delays_in_trials[k] = sqrt(std / (num_spiked_neurons - 1));
			}
			else
				std_first_soma_spike_delays_in_trials[k] = -1;
		}
		else
			mean_first_soma_spike_delays_in_trials[k] = -1;
			
		int num_bursted_neurons = static_cast<int>(dend_spike_delays.size());
		
		if ( num_bursted_neurons > 0 )
		{
			double mean = std::accumulate(dend_spike_delays.begin(), dend_spike_delays.end(), 0.0) / static_cast<double>(num_bursted_neurons);
			mean_dend_spike_delays_in_trials[k] = mean;
			
			if ( num_bursted_neurons > 1 )
			{
				double std = 0;
				std::for_each(dend_spike_delays.begin(), dend_spike_delays.end(),
									[&](double &x){std += (x - mean)*(x - mean);});
									
				std_dend_spike_delays_in_trials[k] = sqrt(std / (num_bursted_neurons - 1));
			}
			else
				std_dend_spike_delays_in_trials[k] = -1;
		}
		else
			mean_dend_spike_delays_in_trials[k] = -1;	
			
		total_num_somatic_spikes_in_trials[k] = std::accumulate(num_somatic_spikes.begin(), num_somatic_spikes.end(), 0);
		total_num_dendritic_spikes_in_trials[k] = std::accumulate(num_dendritic_spikes.begin(), num_dendritic_spikes.end(), 0);
		
		std::for_each(num_somatic_spikes.begin(), num_somatic_spikes.end(), 
							[&num_soma_spiked_in_trials, k](const int &x){if (x > 0) num_soma_spiked_in_trials[k]++;});
		
		std::for_each(num_dendritic_spikes.begin(), num_dendritic_spikes.end(), 
							[&num_dend_spiked_in_trials, k](const int &x){if (x > 0) num_dend_spiked_in_trials[k]++;});
		
		
		weight_ee += dweight_ee;
		
		for (int i = 0; i < num_source_neurons; i++)
			source_neurons[i].reset_time();
			
		for (int i = 0; i < num_target_neurons; i++)
			target_neurons[i].reset_time();
	
		interneuron.reset_time();
	}
	
	std::cout << "Weights_ee:\n";
	for (int i = 0; i < num_steps; i++)
		std::cout << weights[i] << ", ";
	std::cout << std::endl;
	
	std::cout << "Num neurons with somatic spikes:\n";
	for (int i = 0; i < num_steps; i++)
		std::cout << num_soma_spiked_in_trials[i] << ", ";
	std::cout << std::endl;
	
	std::cout << "Probability to produce somatic spike:\n";
	for (int i = 0; i < num_steps; i++)
		std::cout << static_cast<double>(num_soma_spiked_in_trials[i]) / static_cast<double>(num_target_neurons) << ", ";
	std::cout << std::endl;
		
	std::cout << "Total number of somatic spikes:\n";
	for (int i = 0; i < num_steps; i++)
		std::cout << total_num_somatic_spikes_in_trials[i] << ", ";
	std::cout << std::endl;
	
	
	std::cout << "Mean somatic spike delays:\n";
	for (int i = 0; i < num_steps; i++)
		std::cout << mean_first_soma_spike_delays_in_trials[i] << ", ";
	std::cout << std::endl;
	
	std::cout << "Std of somatic spike delays:\n";
	for (int i = 0; i < num_steps; i++)
		std::cout << std_first_soma_spike_delays_in_trials[i] << ", ";
	std::cout << std::endl;
	
	std::cout << "Num neurons with dendritic spikes:\n";
	for (int i = 0; i < num_steps; i++)
		std::cout << num_dend_spiked_in_trials[i] << ", ";
	std::cout << std::endl;
	
	std::cout << "Probability to produce dendritic spike:\n";
	for (int i = 0; i < num_steps; i++)
		std::cout << static_cast<double>(num_dend_spiked_in_trials[i]) / static_cast<double>(num_target_neurons) << ", ";
	std::cout << std::endl;
		
	std::cout << "Total number of dendritic spikes:\n";
	for (int i = 0; i < num_steps; i++)
		std::cout << total_num_dendritic_spikes_in_trials[i] << ", ";
	std::cout << std::endl;
	
	std::cout << "Mean dendritic spike delays:\n";
	for (int i = 0; i < num_steps; i++)
		std::cout << mean_dend_spike_delays_in_trials[i] << ", ";
	std::cout << std::endl;
	
	std::cout << "Std of dendritic spike delays:\n";
	for (int i = 0; i < num_steps; i++)
		std::cout << std_dend_spike_delays_in_trials[i] << ", ";
	std::cout << std::endl;
	
	
	return 0;
}
