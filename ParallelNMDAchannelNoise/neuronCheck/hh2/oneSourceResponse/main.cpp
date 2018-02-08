#include "HH2_buffer.h"
#include "poisson_noise.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <sys/stat.h>

int main(int argc, char** argv)
{
	std::string filename;
	
	if ( argc > 1 )
		filename = argv[1];
	else
	{
		std::cerr << "Usage: <filename>!" << std::endl;
		return -1;
	}
	
	// if file already exists, delete it
	struct stat buf;
	
	if ( stat(filename.c_str(), &buf) == 0 )
		std::remove(filename.c_str());
	
	// excitatory weight from source to target
	double w = 65.0;
	
	double trial_duration = 500; // trial duration in ms
	double timestep = 0.02; // timestep in neuron dynamics
	double wait_time = 200; // time to wait before the kick
	
	int source_kick_time = static_cast<int>(wait_time / timestep); // time when source neuron is excited
	double source_kick = 3.0;
	
	// noise parameters
	double white_noise_mean_soma = 0;
	double white_noise_std_soma = 0.015;
	double white_noise_mean_dend = 0;
	double white_noise_std_dend = 0.035;
	
	Poisson_noise noise_generator;
	unsigned seed = 1170; // 1170 for kick = 3.0; 1176 for kick = 10.0
	
	noise_generator.set_seed(seed);
	
	// target neuron model parameters
	double E_GABA_MATURE = -80.000000;
	double E_GABA_IMMATURE = -55.000000;

	double E_REST_MATURE = -80.000000;
	double E_REST_IMMATURE = -65.000000;

	double AD_MATURE = 10000.000000;
	double AD_IMMATURE = 1000.000000;

	double GK_MATURE = 8.000000;
	double GK_IMMATURE = 16.000000;

	double GNA_MATURE = 60.000000;
	double GNA_IMMATURE = 40.000000;

	double RC_MATURE = 55.000000;
	double RC_IMMATURE = 1.000000;

	double age = 0.0;
	
	double gaba_potential = E_GABA_MATURE + (E_GABA_IMMATURE - E_GABA_MATURE) * exp(-age);
	double rest_potential = E_REST_MATURE + (E_REST_IMMATURE - E_REST_MATURE) * exp(-age);
	double Gk = GK_MATURE + (GK_IMMATURE - GK_MATURE) * exp(-age);
	double GNa = GNA_MATURE + (GNA_IMMATURE - GNA_MATURE) * exp(-age);
	double Ad = AD_MATURE + (AD_IMMATURE - AD_MATURE) * exp(-age);
	double Rc = RC_MATURE + (RC_IMMATURE - RC_MATURE) * exp(-age);
	
	std::cout << "E_GABA = " << gaba_potential << "\n"
			  << "E_REST = " << rest_potential << "\n"
			  << "G_K = "    << Gk             << "\n"
			  << "G_NA = "   << GNa            << "\n"
			  << "Ad = "     << Ad             << "\n"
			  << "Rc = "     << Rc             << std::endl;
	
	HH2_buffer source, target;
	
	// set parameters for target neuron
	target.set_recording_full(filename);
	target.set_noise_generator(&noise_generator);
	target.set_white_noise(white_noise_mean_soma, white_noise_std_soma,
						   white_noise_mean_dend, white_noise_std_dend);
							   
	target.set_dynamics(timestep);
	
	target.set_Ei(gaba_potential);
	target.set_Erest(rest_potential);
	target.set_Gk(Gk);
	target.set_GNa(GNa);
	target.set_Ad(Ad);
	target.set_Rc(Rc);
	
	// set parameters for source neuron
	source.set_noise_generator(&noise_generator);
	source.set_white_noise(white_noise_mean_soma, white_noise_std_soma,
						   white_noise_mean_dend, white_noise_std_dend);
							   
	source.set_dynamics(timestep);
	
	int num_spikes = 0;

	for (int t = 0; t < static_cast<int>(trial_duration / timestep); t++)
	{
		source.Debraband_step_no_target_update();
		target.Debraband_step_no_target_update();
		
		if ( t == source_kick_time	)
			source.raiseE(source_kick);
		
		if ( source.get_fired_dend() )
			std::cout << "Source neuron bursted at " << static_cast<double>(t) * timestep << "\n";
			
		if ( source.get_fired_soma() )
		{
			std::cout << "Source neuron spiked at " << static_cast<double>(t) * timestep << "\n";
			target.raiseExcWeight(w);
		}
		
		if ( target.get_fired_dend() )
			std::cout << "Target neuron bursted at " << static_cast<double>(t) * timestep << "\n";
			
		if ( target.get_fired_soma() )
		{
			std::cout << "Target neuron spiked at " << static_cast<double>(t) * timestep << "\n";
			num_spikes++;
		}
	}
	
	std::cout << "Spike frequency of target neuron = " << static_cast<double>(num_spikes) / (trial_duration / 1000.0) << " Hz" << std::endl;
	
	return 0;
}
