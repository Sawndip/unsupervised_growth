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
	
	// inhibitory kick
	double G_kick = 0.0;
	
	double trial_duration = 10000; // trial duration in ms
	double timestep = 0.02; // timestep in neuron dynamics
	double wait_time = 50; // time to wait before the kick
	int kick_time = static_cast<int>(wait_time / timestep); // time when inhibitory kick is delivered to neurons
	
	// noise parameters
	double white_noise_mean_soma = 0;
	double white_noise_std_soma = 0.025;
	double white_noise_mean_dend = 0;
	double white_noise_std_dend = 0.05;
	
	Poisson_noise noise_generator;
	unsigned seed = 1170; // 1170 for kick = 3.0; 1176 for kick = 10.0
	
	noise_generator.set_seed(seed);
	
	double E_GABA = -50; // GABA reversal potential
	
	HH2_buffer neuron;
	
	// prepare neuron for simulations
	neuron.set_recording_full(filename);
	neuron.set_noise_generator(&noise_generator);
	neuron.set_white_noise(white_noise_mean_soma, white_noise_std_soma,
						   white_noise_mean_dend, white_noise_std_dend);
							   
	neuron.set_dynamics(timestep);
	neuron.set_Ei(E_GABA);

	int num_spikes = 0;

	for (int t = 0; t < static_cast<int>(trial_duration / timestep); t++)
	{
		neuron.Debraband_step_no_target_update();
		
		if ( t == kick_time	)
			neuron.raiseI(G_kick);
		
		if ( neuron.get_fired_dend() )
			std::cout << "Neuron bursted at " << static_cast<double>(t) * timestep << "\n";
			
		if ( neuron.get_fired_soma() )
		{
			std::cout << "Neuron spiked at " << static_cast<double>(t) * timestep << "\n";
			num_spikes++;
		}
	}
	
	std::cout << "Spike frequency = " << static_cast<double>(num_spikes) / (trial_duration / 1000.0) << " Hz" << std::endl;
	
	return 0;
}
