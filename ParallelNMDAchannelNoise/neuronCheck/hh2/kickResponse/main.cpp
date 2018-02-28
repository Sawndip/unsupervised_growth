#include "HH2_buffer.h"
#include "poisson_noise.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <sstream>
#include <sys/stat.h>
#include <iomanip>

int main(int argc, char** argv)
{
	
	int num_groups = 20; // number of groups with different kick strength
	int num_neurons_in_group = 10; // number of neurons in groups
	int num_neurons = num_groups * num_neurons_in_group; // total number of HVC-RA neurons
	
	std::vector<std::vector<double>> burst_delays(num_groups); // time between inhibitory kick and burst
	std::vector<std::vector<double>> spike_delays(num_groups); // time between inhibitory kick and first spike time
	std::vector<std::vector<double>> firing_probability(num_groups); // probability of neuron to fire
	
	std::vector<double> mean_spike_delays(num_groups); // average spike time measured from the kick time
	std::vector<double> std_spike_delays(num_groups); // std of spike times measured from the kick time
		
	std::vector<HH2_buffer> neurons(num_neurons); // HVC-RA neurons
	std::vector<double> G_kick(num_neurons); // vector with strengths of inhibitory kicks
	std::vector<bool> b_spiked(num_neurons); // indicators that neuron already spiked
	
	// increment in strength of inhibitory kick
	double dG_kick = 5.0;
	double G_start = 0.0;
	
	// populate arrays with kicks
	for (int i = 0; i < num_groups; i++)
		for (int j = 0; j < num_neurons_in_group; j++)
			G_kick[i*num_neurons_in_group + j] = G_start + dG_kick * static_cast<double>(i);
	
	double trial_duration = 200; // trial duration in ms
	double timestep = 0.02; // timestep in neuron dynamics
	double wait_time = 100; // time to wait before the kick
	int kick_time = static_cast<int>(wait_time / timestep); // time when inhibitory kick is delivered to neurons
	
	// noise parameters
	double white_noise_mean_soma = 0;
	double white_noise_std_soma = 0.1; // 0.01
	double white_noise_mean_dend = 0;
	double white_noise_std_dend = 0.15; // 0.025
	
	Poisson_noise noise_generator;
	unsigned seed = 1991;
	
	noise_generator.set_seed(seed);
	
	double E_GABA_MATURE = -80.000000;
	double E_GABA_IMMATURE = -55.000000; // was -55.0

	double E_REST_MATURE = -80.000000;
	double E_REST_IMMATURE = -80.000000; // was -65.0

	double AD_MATURE = 10000.000000;
	double AD_IMMATURE = 1000.000000;

	double GK_MATURE = 8.000000;
	double GK_IMMATURE = 8.000000; // was 16.0

	double GNA_MATURE = 60.000000;
	double GNA_IMMATURE = 60.000000; // was 40.0

	double RC_MATURE = 55.000000;
	double RC_IMMATURE = 5.500000; // was 1.0

	double GCA_MATURE = 55.0;
	double GCA_IMMATURE = 0.0;

	double GCAK_MATURE = 150.0;
	double GCAK_IMMATURE = 0.0;

	double GSL_MATURE = 0.1;
	double GSL_IMMATURE = 0.1;

	double GDL_MATURE = 0.1;
	double GDL_IMMATURE = 0.1;

	double age = 0.0;

	double gaba_potential = E_GABA_MATURE + (E_GABA_IMMATURE - E_GABA_MATURE) * exp(-age);
	double rest_potential = E_REST_MATURE + (E_REST_IMMATURE - E_REST_MATURE) * exp(-age);
	double Gk = GK_MATURE + (GK_IMMATURE - GK_MATURE) * exp(-age);
	double GNa = GNA_MATURE + (GNA_IMMATURE - GNA_MATURE) * exp(-age);
	double Ad = AD_MATURE + (AD_IMMATURE - AD_MATURE) * exp(-age);
	double Rc = RC_MATURE + (RC_IMMATURE - RC_MATURE) * exp(-age);
	double GCa = GCA_MATURE + (GCA_IMMATURE - GCA_MATURE) * exp(-age);
	double GCaK = GCAK_MATURE + (GCAK_IMMATURE - GCAK_MATURE) * exp(-age);
	double GsL = GSL_MATURE + (GSL_IMMATURE - GSL_MATURE) * exp(-age);
	double GdL = GDL_MATURE + (GDL_IMMATURE - GDL_MATURE) * exp(-age);
	
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
	
	
	// prepare neurons for simulations
	for (size_t i = 0; i < neurons.size(); i++)
	{
		neurons[i].set_noise_generator(&noise_generator);
		neurons[i].set_white_noise(white_noise_mean_soma, white_noise_std_soma,
								   white_noise_mean_dend, white_noise_std_dend);
								   
		neurons[i].set_dynamics(timestep);
		neurons[i].set_Ei(gaba_potential);
		neurons[i].set_Erest(rest_potential);
		neurons[i].set_Gk(Gk);
		neurons[i].set_GNa(GNa);
		neurons[i].set_GCa(GCa);
		neurons[i].set_GCaK(GCaK);
		neurons[i].set_GsL(GsL);
		neurons[i].set_GdL(GdL);
		neurons[i].set_Ad(Ad);
		neurons[i].set_Rc(Rc);
	}
	
	std::string dirname = "/home/eugene/Output/neuronTest/kickResponse/";
	
	
	
	int num_trials = 10;
	
	
	
	for (int i = 0; i < num_trials; i++)
	{
		std::cout << "Trial " << i << std::endl;
		
		if (i == num_trials - 1)
		{
			for (int j = 0; j < num_neurons; j++)
			{
				// record from all neurons
				int neuron_id = j;
				int inside_group_id = j % num_neurons_in_group;
				
				std::stringstream sstream;
				sstream << std::fixed << std::setprecision(2) << G_kick[neuron_id];
				
				std::string filename = dirname + "Ginh_" + sstream.str() + "_" + std::to_string(inside_group_id) + ".bin";
				
				// if file already exists, delete it
				struct stat buf;
				
				if ( stat(filename.c_str(), &buf) == 0 )
					std::remove(filename.c_str());
					
				neurons[neuron_id].set_recording_full(filename);	
			}
			
			//~ for (int j = 0; j < num_groups; j++)
			//~ {
				//~ // record from the first neuron in each group
				//~ int neuron_id = j * num_neurons_in_group;
				//~ 
				//~ std::stringstream sstream;
				//~ sstream << std::fixed << std::setprecision(2) << G_kick[neuron_id];
				//~ 
				//~ std::string filename = dirname + "Ginh_" + sstream.str() + ".bin";
				//~ 
				//~ // if file already exists, delete it
				//~ struct stat buf;
				//~ 
				//~ if ( stat(filename.c_str(), &buf) == 0 )
					//~ std::remove(filename.c_str());
					//~ 
				//~ neurons[neuron_id].set_recording_full(filename);	
			//~ }
		}
		
		std::fill(b_spiked.begin(), b_spiked.end(), false);
		
		for (int t = 0; t < static_cast<int>(trial_duration / timestep); t++)
		{
			for (int j = 0; j < num_neurons; j++)
			{
				neurons[j].Debraband_step_no_target_update();
				
				
				if ( t == kick_time	)
					neurons[j].raiseI(G_kick[j]);
				
				if ( neurons[j].get_fired_dend() )
				{
					//std::cout << "Neuron " << j << " bursted at " << static_cast<double>(t) * timestep << "\n";
					burst_delays[j / num_neurons_in_group].push_back(static_cast<double>(t) * timestep - wait_time);
				}
				
				if ( ( neurons[j].get_fired_soma() ) && ( !b_spiked[j] ) )
				{
					//std::cout << "Neuron " << j << " first spiked at " << static_cast<double>(t) * timestep << "\n";
					spike_delays[j / num_neurons_in_group].push_back(static_cast<double>(t) * timestep - wait_time);
					b_spiked[j] = true;
				}
			}	
		}
		
		// reset neurons
		for (int j = 0; j < num_neurons; j++)
			neurons[j].set_to_rest();
	}
	
	std::fill(mean_spike_delays.begin(), mean_spike_delays.end(), -1.0);
	std::fill(std_spike_delays.begin(), std_spike_delays.end(), -1.0);
	
	// analyze bursts
	for (int i = 0; i < num_groups; i++)
	{
		std::cout << "Inhibitory kick Gi = " << G_kick[i*num_neurons_in_group] << " number of bursts: " << burst_delays[i].size() 
				  << " out of " << num_neurons_in_group * num_trials << "\n"
				  << "number of spikes: " << spike_delays[i].size() << std::endl;
		
		if ( !burst_delays[i].empty() )
		{
			double mean_burst_delay = std::accumulate(burst_delays[i].begin(), burst_delays[i].end(), 0.0) / static_cast<double>(burst_delays[i].size());
			
			std::cout << "Mean burst delay = " << mean_burst_delay << "\n" << std::endl;
			
			// calculate standard deviation of burst delays
			if ( burst_delays[i].size() > 1 )
			{
				double std_burst_delay = 0;
				
				for (size_t j = 0; j < burst_delays[i].size(); j++)
					std_burst_delay += ( burst_delays[i][j] - mean_burst_delay ) * ( burst_delays[i][j] - mean_burst_delay);
					
				std_burst_delay = sqrt( std_burst_delay / static_cast<double>(burst_delays[i].size() - 1));
				
				std::cout << "Std burst delay = " << std_burst_delay << "\n" << std::endl;
			}
			
		}
		
		// analyze spikes
		if ( !spike_delays[i].empty() )
		{
			mean_spike_delays[i] = std::accumulate(spike_delays[i].begin(), spike_delays[i].end(), 0.0) / static_cast<double>(spike_delays[i].size());
			
			std::cout << "Mean spike delay = " << mean_spike_delays[i] << "\n" << std::endl;
			
			// calculate standard deviation of burst delays
			if ( spike_delays[i].size() > 1 )
			{
				double std_spike_delay = 0;
				
				for (size_t j = 0; j < spike_delays[i].size(); j++)
					std_spike_delay += ( spike_delays[i][j] - mean_spike_delays[i] ) * ( spike_delays[i][j] - mean_spike_delays[i]);
					
				std_spike_delays[i] = sqrt( std_spike_delay / static_cast<double>(spike_delays[i].size() - 1));
				
				std::cout << "Std spike delay = " << std_spike_delays[i] << "\n" << std::endl;
			}
		}
		
	}
	
	// print results
	
	std::cout << "\nInhibitory kicks = ";
	
	for (int i = 0; i < num_groups-1; i++)
		std::cout << G_kick[i*num_neurons_in_group] << ", ";
	std::cout << G_kick[(num_groups-1)*num_neurons_in_group] << std::endl;;
	
	std::cout << "\nNum spikes = ";
	
	for (int i = 0; i < num_groups-1; i++)
		std::cout << static_cast<double>(spike_delays[i].size()) / static_cast<double>(num_trials * num_neurons_in_group) << ", ";
	std::cout << static_cast<double>(spike_delays[num_groups-1].size()) / static_cast<double>(num_trials * num_neurons_in_group) << std::endl;;
	
	std::cout << "\nMean spike delay = ";
	
	for (int i = 0; i < num_groups-1; i++)
		std::cout << mean_spike_delays[i] << ", ";
	std::cout << mean_spike_delays[num_groups-1] << std::endl;;
	
	std::cout << "\nStd of spike delays = ";
	
	for (int i = 0; i < num_groups-1; i++)
		std::cout << std_spike_delays[i] << ", ";
	std::cout << std_spike_delays[num_groups-1] << std::endl;;
	
	
	return 0;
}
