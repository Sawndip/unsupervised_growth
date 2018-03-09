#include "HHI_buffer.h" 
#include "HH2_clopath.h"
#include <iostream>
#include "poisson_noise.h"
#include <vector>
#include <algorithm>
#include <sys/stat.h>
#include <string>
#include <fstream>

void write_clopath_sim_results(std::vector<double>& time, 
							   std::vector<double>& vd, 
							   std::vector<double>& u_minus,
							   std::vector<double>& u_plus,
							   std::vector<double>& x,
							   std::vector<double>& w,
							   const char *filename){
							
	std::ofstream out;
	
	out.open(filename, std::ios::out | std::ios::binary);
	
	int num_datapoints = static_cast<int>(vd.size());
	
	out.write(reinterpret_cast<char*>(&num_datapoints), sizeof(int));
	
	for (int i = 0; i < num_datapoints; i++){
		out.write(reinterpret_cast<char*>(&time[i]), sizeof(double));
		out.write(reinterpret_cast<char*>(&vd[i]), sizeof(double));
		out.write(reinterpret_cast<char*>(&u_minus[i]), sizeof(double));
		out.write(reinterpret_cast<char*>(&u_plus[i]), sizeof(double));
		out.write(reinterpret_cast<char*>(&x[i]), sizeof(double));
		out.write(reinterpret_cast<char*>(&w[i]), sizeof(double));
	}
	
	out.close();			
}

int main()
{
	
	// trial parameters
	double TIMESTEP = 0.02;
	
	double E_GABA_MATURE = -80.000000;
	double E_GABA_IMMATURE = -55.000000;

	double E_REST_MATURE = -80.000000;
	double E_REST_IMMATURE = -80.000000;

	double AD_MATURE = 10000.000000;
	double AD_IMMATURE = 1000.000000;

	double GK_MATURE = 8.000000;
	double GK_IMMATURE = 8.000000;

	double GNA_MATURE = 60.000000;
	double GNA_IMMATURE = 60.000000;

	double RC_MATURE = 55.000000;
	double RC_IMMATURE = 5.50000;
	
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
			  
	
	HH2_clopath source_neuron, target_neuron;
	
	Poisson_noise noise_generator;
	unsigned seed = 1991;
	
	noise_generator.set_seed(seed);
	
	// set parameters for neurons
	source_neuron.set_dynamics(TIMESTEP);
	source_neuron.set_noise_generator(&noise_generator);
	
	source_neuron.set_Ei(gaba_potential);
	source_neuron.set_Erest(rest_potential);
	source_neuron.set_Gk(Gk);
	source_neuron.set_GNa(GNa);
	source_neuron.set_GCa(GCa);
	source_neuron.set_GCaK(GCaK);
	source_neuron.set_GsL(GsL);
	source_neuron.set_GdL(GdL);
	source_neuron.set_Ad(Ad);
	source_neuron.set_Rc(Rc);
	
	target_neuron.set_dynamics(TIMESTEP);
	target_neuron.set_noise_generator(&noise_generator);
	
	target_neuron.set_Ei(gaba_potential);
	target_neuron.set_Erest(rest_potential);
	target_neuron.set_Gk(Gk);
	target_neuron.set_GNa(GNa);
	target_neuron.set_GCa(GCa);
	target_neuron.set_GCaK(GCaK);
	target_neuron.set_GsL(GsL);
	target_neuron.set_GdL(GdL);
	target_neuron.set_Ad(Ad);
	target_neuron.set_Rc(Rc);
		
	
	std::string dirname = "/home/eugene/Output/clopath/frequency/preSpikePostSpike/";
	
	double W_TRAINING_KICK = 1000.0; // strength of excitatory kick delivered to training neurons	 
	double T_TRAINING_KICK = 150; // time when training neurons receive excitatory kick
	
	double W_TARGET_KICK = 1000;
	
	int num_pairs = 60;
	double frequency = 20.0; // in Hz
	
	double TRIAL_DURATION = T_TRAINING_KICK + 1000 * num_pairs / frequency + 50.0;
	
	
	int num_steps = 80;
	double time_kick_step = 0.25;
	double dt_min = - static_cast<double>(num_steps / 2) * time_kick_step;
	
	int num_timepoints = static_cast<int>(TRIAL_DURATION / TIMESTEP);
	
	std::vector<double> x(num_timepoints); // filtered glutamate bounded to postsynaptic receptor
	std::vector<double> vd(num_timepoints);
	std::vector<double> u_minus(num_timepoints);
	std::vector<double> u_plus(num_timepoints);
	std::vector<double> w(num_timepoints);
	std::vector<double> time(num_timepoints);
	
	double DX = 1;
	
	double THETA_MINUS = -70.0;
	double THETA_PLUS = -40.0;
	
	double A_LTD = 0.5e0; // was 14e0
	double A_LTP = 0.5e0; // was 8e0
	double TAU_X = 15.0;
	
	double start_weight = 100.0;
	
	std::vector<double> weights(num_steps);
	std::vector<double> dt(num_steps);
	
	
	for (int i = 0; i < num_timepoints; i++)
		time[i] = TIMESTEP*i;
	
	
	std::fill(weights.begin(), weights.end(), start_weight);
	
	for (int k = 0; k < num_steps; k++)
	{	
		//~ std::string filename = "/mnt/hodgkin_home/eugene/Output/neuronTest/inhAndExcInputsResponse/RA" + std::to_string(k) + ".bin";
	//~ 
		//~ struct stat buf;
	//~ 
		//~ if ( stat(filename.c_str(), &buf) == 0 )
			//~ std::remove(filename.c_str());
	//~ 
		//~ target_neuron.set_recording_full(filename.c_str());
		//~ 
		target_neuron.set_to_rest();
		source_neuron.set_to_rest();
		
		
		std::fill(x.begin(), x.end(), -1000.0);
		std::fill(vd.begin(), vd.end(), -1000.0);
		std::fill(u_minus.begin(), u_minus.end(), -1000.0);
		std::fill(u_plus.begin(), u_plus.end(), -1000.0);
		std::fill(w.begin(), w.end(), -1000.0);
		
		x[0] = 0;
		vd[0] = target_neuron.get_Vd();
		u_minus[0] = target_neuron.get_u_minus();
		u_plus[0] = target_neuron.get_u_plus();
		w[0] = start_weight;
		
		
		bool source_spiked = false;
		bool target_spiked = false;
		
		double first_spike_time_source = -1.0;
		double first_spike_time_target = -1.0;
	
		double t_training_kick = T_TRAINING_KICK;
		double t_target_kick = t_training_kick + dt_min + time_kick_step * static_cast<double>(k);
	
		
		// trial
		for (int t = 1; t < static_cast<int>(TRIAL_DURATION / TIMESTEP); t++)
		{
			if ( static_cast<double>(t) * TIMESTEP > t_training_kick )
			{
				//source_neuron.raiseE(G_TRAINING_KICK);
				source_neuron.raiseExcWeight(W_TRAINING_KICK);
				t_training_kick += 1000 / frequency;
			}
				
			if ( static_cast<double>(t) * TIMESTEP > t_target_kick )
			{
				target_neuron.raiseExcWeight(W_TARGET_KICK);
				t_target_kick += 1000 / frequency;
			}
		
			source_neuron.Debraband_step_no_target_update();
			target_neuron.Debraband_step_no_target_update();
			
			x[t] = x[t-1] * exp(-TIMESTEP / TAU_X);
			w[t] = w[t-1];
			u_minus[t] = target_neuron.get_u_minus();
			u_plus[t] = target_neuron.get_u_plus();	
			vd[t] = target_neuron.get_Vd();
			
			// if source neuron fired, update glutamate and apply LTD
			if ( source_neuron.get_fired_soma() )
			{
				x[t] += DX;
				std::cout << "source neuron somatic spike at " << static_cast<double>(t) * TIMESTEP << std::endl;
				
				
				std::cout << "u_minus = " << u_minus[t] << std::endl;
				
				if ( u_minus[t] >= THETA_MINUS ){
					weights[k] -= A_LTD;
					w[t] -= A_LTD;
				}
					
				if ( !source_spiked){
					first_spike_time_source = static_cast<double>(t) * TIMESTEP;
					source_spiked = true;
				}
					
			}
			
			
		
			
			//std::cout << "u_plus = " << u_plus << std::endl;
			//std::cout << "vd = " << vd << std::endl;
			
			if ( ( u_plus[t] >= THETA_MINUS ) && ( vd[t] >= THETA_PLUS ) ){
				weights[k] += TIMESTEP * A_LTP * x[t];
				w[t] += TIMESTEP * A_LTP * x[t];
			}
			
			if ( target_neuron.get_fired_soma() )
			{
				if ( !target_spiked){
					first_spike_time_target = static_cast<double>(t) * TIMESTEP;
					target_spiked = true;
				}
				
				std::cout << "target neuron somatic spike at " << static_cast<double>(t) * TIMESTEP << std::endl;
			}
			
		}
		
		
		
		t_target_kick += time_kick_step;
		
		if ( ( first_spike_time_target > 0 ) && ( first_spike_time_source > 0 ) ){
			dt[k] = first_spike_time_target - first_spike_time_source;
			std::string filename = dirname + std::to_string(k) + ".bin";
		
			write_clopath_sim_results(time, vd, u_minus, u_plus, x, w, filename.c_str());
		}
	}
	
	std::cout << "Time difference in first spike times:\n";
	for (int k = 0; k < num_steps; k++)
		std::cout << dt[k] << ", ";
	std::cout << "\n" << std::endl;
	
	std::cout << "Weights:\n";
	for (int k = 0; k < num_steps; k++)
		std::cout << weights[k] << ", ";
	std::cout << "\n" << std::endl;
	
	return 0;
}
