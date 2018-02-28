#include <iostream>
#include "HH2_buffer.h"
#include "poisson_noise.h"
#include <fstream>
#include <vector>
#include <sys/stat.h>


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


int main(){
	double timestep = 0.02;
	double sim_duration = 300;
	double kick_time = 100.0;
	bool excited;
	bool spiked;
	
	std::string dirname = "/home/eugene/Output/neuronTest/modelStability/";
	
	unsigned seed = 1991;
	
	Poisson_noise noise_generator;
	
	noise_generator.set_seed(seed);
	
	
	int num_iter = static_cast<int>(sim_duration / timestep);
	
	int N_TR = 10;
	
	std::vector<HH2_buffer> training(N_TR);
	
	for (int i = 0; i < N_TR; i++){
		training[i].set_dynamics(timestep);
		training[i].set_noise_generator(&noise_generator);
		training[i].set_white_noise(0.0, 0.010, 0.0, 0.030);
	}
	
	int N = 1000; 
	
	std::vector<double> burst_onset_times(N);
	std::vector<double> spike_times(N);
	
	std::fill(burst_onset_times.begin(), burst_onset_times.end(), -1.0);
	std::fill(spike_times.begin(), spike_times.end(), -1.0);
	
	double W = 20.0;
	double W_step = 0.0;
	double G_kick = 3.0;
	
	
	
	
	for (int i = 0; i < N; i++)
	{
		HH2_buffer *neuron = new HH2_buffer;
	
		neuron->set_dynamics(timestep);
		
		neuron->set_Ei(gaba_potential);
		neuron->set_Erest(rest_potential);
		neuron->set_Gk(Gk);
		neuron->set_GNa(GNa);
		neuron->set_Ad(Ad);
		neuron->set_Rc(Rc);
		
		
		neuron->set_noise_generator(&noise_generator);
		neuron->set_white_noise(0.0, 0.010, 0.0, 0.030);
	
		excited = false;
		spiked = false;
	
		std::string filename = dirname + "RA" + std::to_string(i) + ".bin";
		
		struct stat buf;
	
		if ( stat(filename.c_str(), &buf) == 0 )
			std::remove(filename.c_str());
		
		
		//neuron->set_recording_full(filename);
		
	
	
		for (int j = 0; j < num_iter; j++)
		{

			neuron->Debraband_step_no_target_update();	
			
			if ( neuron->check_bad_values() < 0 )
				std::cout << "Bad values for weight = " << W << std::endl;

			if ( ( static_cast<double>(j) * timestep > kick_time ) && ( !excited ) ){
				for (int k = 0; k < N_TR; k++){
					
					training[k].raiseE(G_kick);
					excited = true;
				}
			}


			
			for (int k = 0; k < N_TR; k++){
				
				training[k].Debraband_step_no_target_update();	
			
				if ( training[k].get_fired_soma() )
					neuron->raiseExcWeight(W);
			}
			
						
			//~ if ( neuron->get_fired_dend() )
				//~ burst_onset_times[i] = static_cast<double>(j) * timestep - kick_time;
				//~ 
			//~ if ( ( neuron->get_fired_soma() ) && ( !spiked ) )
			//~ {
				//~ spike_times[i] = static_cast<double>(j) * timestep - kick_time;
				//~ spiked = true;
			//~ }
		}
		
		
		delete neuron;
		
		for (int k = 0; k < N_TR; k++)
			training[k].reset_time();
		
		
		W += W_step;
	}
	
	return 0;
}
