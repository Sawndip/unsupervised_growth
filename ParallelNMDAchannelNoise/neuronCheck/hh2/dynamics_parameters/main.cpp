#include "HH2_test.h"
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
	
	double trial_duration = 100000; // trial duration in ms
	double timestep = 0.02; // timestep in neuron dynamics
	double wait_time = 50; // time to wait before the kick
	int kick_time = static_cast<int>(wait_time / timestep); // time when inhibitory kick is delivered to neurons
	
	// noise parameters
	double white_noise_mean_soma = 0;
	double white_noise_std_soma = 0.1;
	double white_noise_mean_dend = 0;
	double white_noise_std_dend = 0.2;
	
	double age = 0.3;
	
	
	double E_GABA = -80.000000;
	
	double E_REST_MATURE = -80.000000;
	double E_REST_IMMATURE = -55.000000;

	double AD = 10000.000000;

	double GK = 8.000000;
	double GNA = 60.000000;

	double RC = 55.000000;

	double GCA_MATURE = 55.0; 
	double GCA_IMMATURE = 0.0;
	
	double GCAK = 150.0;
	double GSL = 0.1;
	double GDL = 0.1;

	// exponential schedule	
//	double GCa = GCA_MATURE + (GCA_IMMATURE - GCA_MATURE) * exp(-age);
//	double E_rest = E_REST_MATURE + (E_REST_IMMATURE - E_REST_MATURE) * exp(-age);
	
	// linear schedule
	double MATURATION_SCALE = 1.0;
	double GCa = GCA_IMMATURE + (GCA_MATURE - GCA_IMMATURE) * age / MATURATION_SCALE;
	double E_rest = E_REST_IMMATURE + (E_REST_MATURE - E_REST_IMMATURE) * age / MATURATION_SCALE;

	Poisson_noise noise_generator;
	unsigned seed = 1170; // 1170 for kick = 3.0; 1176 for kick = 10.0
	
	noise_generator.set_seed(seed);
	
	
	HH2_test neuron;
	
	// prepare neuron for simulations
	//neuron.set_recording_full(filename);
	neuron.set_noise_generator(&noise_generator);
	neuron.set_white_noise(white_noise_mean_soma, white_noise_std_soma,
						   white_noise_mean_dend, white_noise_std_dend);
							   
	neuron.set_dynamics(timestep);
	
	neuron.set_Ei(E_GABA);	
	neuron.set_Erest(E_rest);
	
	neuron.set_Gk(GK);
	neuron.set_GNa(GNA);
	neuron.set_GCa(GCa);
	neuron.set_GCaK(GCAK);
	neuron.set_GsL(GSL);
	neuron.set_GdL(GDL);
	
	neuron.set_Ad(AD);
	neuron.set_Rc(RC);

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
