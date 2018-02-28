#include "HH2_buffer.h"
#include "poisson_noise.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <sys/stat.h>
#include <sstream>
#include <iomanip>
#include <algorithm>

int main(int argc, char** argv)
{
	double trial_duration = 250; // trial duration in ms
	double timestep = 0.02; // timestep in neuron dynamics
	
	double wait_time = 50; // time to wait before the inhibition input
	double inh_duration = 300; // duration of saturated inhibition
	
	std::string dirname = "/home/eugene/Output/neuronTest/saturatedInhibitionResponse/";
	
	int inh_start = static_cast<int>(wait_time / timestep); // time when inhibition is delivered to neuron
	int inh_end = static_cast<int>( (wait_time + inh_duration) / timestep); // time when inhibition ends
	
	// noise parameters
	double white_noise_mean_soma = 0;
	double white_noise_std_soma = 0.0; // was 0.15; last: 0.025
	double white_noise_mean_dend = 0;
	double white_noise_std_dend = 0.0; // was 0.25; last: 0.05
	
	double E_GABA_MATURE = -80.000000;
	double E_GABA_IMMATURE = -51.50000; // was -55.0

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
	
	Poisson_noise noise_generator;
	unsigned seed = 1988; // 1990
	
	noise_generator.set_seed(seed);
	
	int num_neurons = 10;
	
	std::vector<HH2_buffer> neurons(num_neurons); // HVC-RA neurons
	
	// prepare neuron for simulations
	for (int i = 0; i < num_neurons; i++)
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
	
	double G_inh_step = 0.25;
	int num_inh_steps = 100;
	
	double G_inh = 0.0;
	
	std::vector<double> dend_spike_times;
	std::vector<double> soma_spike_times;
	
	std::vector<double> first_soma_spike_times;
	
	std::vector<double> mean_first_soma_spike_times(num_inh_steps);
	std::vector<double> std_first_soma_spike_times(num_inh_steps);
	
	std::vector<bool> b_spiked(num_neurons);
	
	std::vector<double> mean_soma_spike_times(num_inh_steps);
	std::vector<double> std_soma_spike_times(num_inh_steps);
	
	std::vector<double> mean_dend_spike_times(num_inh_steps);
	std::vector<double> std_dend_spike_times(num_inh_steps);
	
	std::vector<double> fraction_soma_spikes(num_inh_steps);
	std::vector<double> fraction_dend_spikes(num_inh_steps);
	
	for (int i = 0; i < num_inh_steps; i++)
	{
		for (int j = 0; j < num_neurons; j++)
			neurons[j].set_to_rest();
	
		std::stringstream sstream;
		sstream << std::fixed << std::setprecision(2) << G_inh;
	
		std::cout << "Inhibition " << sstream.str() << "\n";
		
		// record from the first neuron
		
		std::string filename = dirname + "Ginh_" + sstream.str() + ".bin";
		
		// if file already exists, delete it
		struct stat buf;
		
		if ( stat(filename.c_str(), &buf) == 0 )
			std::remove(filename.c_str());
			
		neurons[0].set_recording_full(filename);
		
		std::fill(b_spiked.begin(), b_spiked.end(), false);
		
		for (int t = 0; t < static_cast<int>(trial_duration / timestep); t++)
		{
			for (int j = 0; j < num_neurons; j++)
			{
				if ( ( t >= inh_start ) && ( t <= inh_end ) )
				{
					neurons[j].set_inh_conductance(G_inh);
					//neuron.set_c(0);
				}	
				
				neurons[j].Debraband_step_no_target_update();
					
					
				if ( neurons[j].get_fired_dend() )
					//std::cout << "Dend spike time at " << static_cast<double>(t) * timestep << "\n";
					dend_spike_times.push_back(static_cast<double>(t) * timestep - wait_time);
			
			
				if ( neurons[j].get_fired_soma() ){
					//std::cout << "Soma spike time at " << static_cast<double>(t) * timestep << "\n";
					soma_spike_times.push_back(static_cast<double>(t) * timestep - wait_time);
					
					if ( !b_spiked[j] ){
						first_soma_spike_times.push_back(static_cast<double>(t) * timestep - wait_time);
						b_spiked[j] = true;
					}
						
				}
			}
		}	
		
		G_inh += G_inh_step;
		
		// analyze somatic and dendritic spike times
		std::cout << "Number of somatic spikes: " << soma_spike_times.size() << std::endl;
		
		fraction_soma_spikes[i] = static_cast<double>(soma_spike_times.size()) / static_cast<double>(num_neurons);
		
		if ( !soma_spike_times.empty() )
		{
			double mean_soma_spike_time = std::accumulate(soma_spike_times.begin(), soma_spike_times.end(), 0.0) / static_cast<double>(soma_spike_times.size());
			std::cout << "Mean somatic spike time = " << mean_soma_spike_time << std::endl;
			
			mean_soma_spike_times[i] = mean_soma_spike_time;
			
			// calculate standard deviation of somatic spike times
			if ( soma_spike_times.size() > 1 )
			{
				double std = 0.0;
				std::for_each(soma_spike_times.begin(), soma_spike_times.end(), 
						[&](double t){std += (t-mean_soma_spike_time)*(t-mean_soma_spike_time);});
						
				std = sqrt(std / static_cast<double>(soma_spike_times.size()-1));
			
				std::cout << "Std of somatic spike times = " << std << std::endl;
				std_soma_spike_times[i] = std;
			}
			else
				std_soma_spike_times[i] = -1.0;
		}
		else
		{
			mean_soma_spike_times[i] = -1.0;
			std_soma_spike_times[i] = -1.0;
		}
		
		fraction_soma_spikes[i] = static_cast<double>(soma_spike_times.size()) / static_cast<double>(num_neurons);
		
		if ( !first_soma_spike_times.empty() )
		{
			double mean_first_soma_spike_time = std::accumulate(first_soma_spike_times.begin(), first_soma_spike_times.end(), 0.0) / static_cast<double>(first_soma_spike_times.size());
			std::cout << "Mean first somatic spike time = " << mean_first_soma_spike_time << std::endl;
			
			mean_first_soma_spike_times[i] = mean_first_soma_spike_time;
			
			// calculate standard deviation of somatic spike times
			if ( first_soma_spike_times.size() > 1 )
			{
				double std = 0.0;
				std::for_each(first_soma_spike_times.begin(), first_soma_spike_times.end(), 
						[&](double t){std += (t-mean_first_soma_spike_time)*(t-mean_first_soma_spike_time);});
						
				std = sqrt(std / static_cast<double>(first_soma_spike_times.size()-1));
			
				std::cout << "Std of first somatic spike times = " << std << std::endl;
				std_first_soma_spike_times[i] = std;
			}
			else
				std_first_soma_spike_times[i] = -1.0;
		}
		else
		{
			mean_first_soma_spike_times[i] = -1.0;
			std_first_soma_spike_times[i] = -1.0;
		}
		
		std::cout << "Number of dendritic spikes: " << dend_spike_times.size() << std::endl;
		
		fraction_dend_spikes[i] = static_cast<double>(dend_spike_times.size()) / static_cast<double>(num_neurons);
		
		if ( !dend_spike_times.empty() )
		{
			double mean_dend_spike_time = std::accumulate(dend_spike_times.begin(), dend_spike_times.end(), 0.0) / static_cast<double>(dend_spike_times.size());
			std::cout << "Mean dendritic spike time = " << mean_dend_spike_time << std::endl;
			mean_dend_spike_times[i] = mean_dend_spike_time;
			
			// calculate standard deviation of dendritic spike times
			if ( dend_spike_times.size() > 1 )
			{
				double std = 0.0;
				std::for_each(dend_spike_times.begin(), dend_spike_times.end(), 
						[&](double t){std += (t-mean_dend_spike_time)*(t-mean_dend_spike_time);});
						
				std = sqrt(std / static_cast<double>(dend_spike_times.size()-1));
				
				
				std::cout << "Std of dendritic spike times = " << std << std::endl;
				std_dend_spike_times[i] = std;
			}
			else
				std_dend_spike_times[i] = -1.0;
			
			
		}
		else
		{
			mean_dend_spike_times[i] = -1.0;
			std_dend_spike_times[i] = -1.0;
		}
		
		soma_spike_times.clear();
		dend_spike_times.clear();
		
		first_soma_spike_times.clear();
	}
	
	// somatic spike info
	std::cout << "\nFraction somatic spikes: " << "\n";
	
	for (int i = 0; i < num_inh_steps-1; i++)
		std::cout << fraction_soma_spikes[i] << ", ";
	std::cout << fraction_soma_spikes[num_inh_steps - 1];
	
	std::cout << "\nMean somatic spike times: " << "\n";
	
	for (int i = 0; i < num_inh_steps-1; i++)
		std::cout << mean_soma_spike_times[i] << ", ";
	std::cout << mean_soma_spike_times[num_inh_steps - 1];
	
	std::cout << "\nStd somatic spike times: " << "\n";
	
	for (int i = 0; i < num_inh_steps-1; i++)
		std::cout << std_soma_spike_times[i] << ", ";
	std::cout << std_soma_spike_times[num_inh_steps - 1];
	
	std::cout << "\nMean first somatic spike times: " << "\n";
	
	for (int i = 0; i < num_inh_steps-1; i++)
		std::cout << mean_first_soma_spike_times[i] << ", ";
	std::cout << mean_first_soma_spike_times[num_inh_steps - 1];
	
	std::cout << "\nStd first somatic spike times: " << "\n";
	
	for (int i = 0; i < num_inh_steps-1; i++)
		std::cout << std_first_soma_spike_times[i] << ", ";
	std::cout << std_first_soma_spike_times[num_inh_steps - 1];
	
	// dendritic spike info	
	std::cout << "\nFraction dendritic spikes: " << "\n";
	
	for (int i = 0; i < num_inh_steps - 1; i++)
		std::cout << fraction_dend_spikes[i] << ", ";
	std::cout << fraction_dend_spikes[num_inh_steps - 1];
	
	std::cout << "\nMean dendritic spike times: " << "\n";
	
	for (int i = 0; i < num_inh_steps-1; i++)
		std::cout << mean_dend_spike_times[i] << ", ";
	std::cout << mean_dend_spike_times[num_inh_steps - 1];
	
	std::cout << "\nStd dendritic spike times: " << "\n";
	
	for (int i = 0; i < num_inh_steps-1; i++)
		std::cout << std_dend_spike_times[i] << ", ";
	std::cout << std_dend_spike_times[num_inh_steps - 1];
	
	return 0;
}
