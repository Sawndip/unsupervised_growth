#include "NetworkGrowthSimulator.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <functional>
#include "training_current.h"
#include <sys/stat.h>
#include <set>
#include <iomanip>

const double NetworkGrowthSimulator::STDP_WINDOW = 100.0;
const double NetworkGrowthSimulator::WAITING_TIME = 100.0;
const double NetworkGrowthSimulator::TIMESTEP = 0.02; // timestep for solving neuronal dynamics in ms
const double NetworkGrowthSimulator::NETWORK_UPDATE_FREQUENCY = 0.25; // frequency of communication between processes in ms
const double NetworkGrowthSimulator::TRIAL_DURATION = 500.0; // duration of simulation trial in ms

const double NetworkGrowthSimulator::EXPERIMENTAL_INTERNEURON_DISTANCE = 40.0; // experimental distance between interneurons in microns
const double NetworkGrowthSimulator::MIN_INTERSOMATIC_DISTANCE = 10.0; // minimum allowed distance between somas of two neurons in microns

const double NetworkGrowthSimulator::G_TRAINING_KICK = 3.0; // strength of excitatory conductance delivered to training neurons
const double NetworkGrowthSimulator::MATURATION_SCALE_SPONTANEOUS = 100000; // timescale of neuron properties of spontaneous maturation
const double NetworkGrowthSimulator::MATURATION_SCALE_DRIVEN = 1000; // timescale of neuron properties of driven maturation
const double NetworkGrowthSimulator::EXPERIMENTAL_AXONAL_DELAY_PER_MICRON = 0.01; // experimental value of axonal time delay in ms for 1 micron distance
const double NetworkGrowthSimulator::MATURE_SYNAPSE_SCALING = 15.0; // constant to scale excitatory synapses to mature neurons

//const double NetworkGrowthSimulator::INHIBITION_SYNAPSE_SCALING = 10.0; // constant to scale excitatory synapses to mature neurons


static double distance2d(double x1, double y1, double x2, double y2)
{
	return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

static double distance_on_sphere(double x1, double y1, double z1, double x2, double y2, double z2)
{
	double tmp = x1*x2 + y1*y2 + z1*z2;
	
	if (tmp > 1.0)
		return 0.0;
	else
		return acos(tmp);
}

using namespace std::placeholders;

NetworkGrowthSimulator::NetworkGrowthSimulator()
{

    MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
    
    // initialize internal time
    internal_time = 0.0;
    network_time = 0.0;
    trial_number = 0;

    this->initialize_generator();

	// set all parameters
	//~ 
	//~ struct SynapticParameters synaptic_params = cfgGrowth.get_synaptic_parameters();
    //~ struct GabaParameters gaba_params = cfgGrowth.get_gaba_parameters();
    //~ struct TimeParameters time_params = cfgGrowth.get_time_parameters();
    //~ struct NoiseParameters noise_params = cfgGrowth.get_noise_parameters();
    //~ 
    //~ network_params = cfgGrowth.get_network_parameters();
    //~ 
     //~ // set inhibitory parameters for network generator
    //~ double Gie_mean, Gie_std; // inhibitory strength parameters
    //~ 
    //~ Gie_mean = network_params.Gie_mean;
    //~ Gie_std = network_params.Gie_std;
    //~ 
    //~ networkGen.set_inhibitory_strength(Gie_mean, Gie_std);
    //~ 
    //~ this->set_synaptic_parameters(synaptic_params);
    //~ this->set_gaba_parameters(gaba_params);
    //~ this->set_time_parameters(time_params);
    //~ this->set_noise_parameters(noise_params);
}

void NetworkGrowthSimulator::generate_synthetic_chain(int num_layers, int num_group, double probability, double Gee)
{
	N_TR = num_group;
	
	training_neurons.resize(N_TR);

	std::iota(training_neurons.begin(), training_neurons.end(), 0);
		
	for (int i = 0; i < num_layers-1; i++)
	{
		for (int j = 0; j < num_group; j++)
		{
			for (int k = 0; k < num_group; k++)
			{
				if ( noise_generator.random(1.0) < probability )
				{
					active_synapses_global[i*num_group+j].push_back((i+1)*num_group+k);
					supersynapses_global[i*num_group+j].push_back((i+1)*num_group+k);
					weights_RA_RA_global[i*num_group+j][(i+1)*num_group+k] = this->sample_G(Gee);
					
				}
			}
		}
	}
}

void NetworkGrowthSimulator::generate_network_topology(int N_ra, int N_i, int N_tr, 
								const struct TopologyParameters &top_params, std::string networkDir)
{
	if (MPI_rank == 0)
	{
		N_RA = N_ra;
		N_I = N_i;
		N_TR = N_tr;
		
		topology_params = top_params;
	
		xx_RA.resize(N_RA);
		yy_RA.resize(N_RA);
		
		xx_I.resize(N_I);
		yy_I.resize(N_I);
		
		training_neurons.resize(N_TR);
		
		if ( topology_params.arrangement == "sphere" )
		{
			zz_RA.resize(N_RA);
			zz_I.resize(N_I);
		}
		
		weights_RA_I_global.resize(N_RA);
		syn_ID_RA_I_global.resize(N_RA);
		syn_lengths_RA_I_global.resize(N_RA);
		axonal_delays_RA_I_global.resize(N_RA);
		
		weights_I_RA_global.resize(N_I);
		syn_ID_I_RA_global.resize(N_I);
		syn_lengths_I_RA_global.resize(N_I);
		axonal_delays_I_RA_global.resize(N_I);
		
		
		this->sample_coordinates();
		this->sample_connections();
		
		
		this->write_different_training_arrangements(networkDir);
		
		this->write_network_topology(networkDir);
		this->write_pajek_topology((networkDir + "network_topology.net").c_str());
	}
}

void NetworkGrowthSimulator::sample_coordinates_for_replaced(std::vector<int>& neurons_to_replace)
{
	std::vector<int> neurons_with_fixed_coordinates;
	
	std::set<int> set_neurons_to_replace(neurons_to_replace.begin(), neurons_to_replace.end());
	
	for (int i = 0; i < N_RA; i++)
		if ( set_neurons_to_replace.find(i) == set_neurons_to_replace.end() )
			neurons_with_fixed_coordinates.push_back(i);
	
	const double pi = 3.14159265358979323846; // value of constant pi
		
	double too_close_distance = MODEL_INTERNEURON_DISTANCE * MIN_INTERSOMATIC_DISTANCE / EXPERIMENTAL_INTERNEURON_DISTANCE;
		
	// HVC(RA) neurons are sampled from uniform distribution on sphere
	for (size_t i = 0; i < neurons_to_replace.size(); i++)
	{
		double tmp_x, tmp_y, tmp_z;
		
		bool too_close = true;
		
		// check if neuron is too close to other interneurons or previously sampled HVC-RA
		while ( too_close )
		{
			too_close = false;
		
			double z = 2*noise_generator.random(1.0) - 1; // generate z in range(-1, 1)
			double phi = noise_generator.random(2.0 * pi); // generate phi in range(0, 2*pi)
		
			tmp_x = sin(acos(z)) * cos(phi);
			tmp_y = sin(acos(z)) * sin(phi);
			tmp_z = z;		
			
			for (int j = 0; j < N_I; j++)
			{
				if ( distance_on_sphere(tmp_x, tmp_y, tmp_z, xx_I[j], yy_I[j], zz_I[j]) < too_close_distance )
				{
					too_close = true;
					break;
				}
			}
			
			if ( !too_close )
			{
				for (size_t j = 0; j < neurons_with_fixed_coordinates.size(); j++)
				{
					int neuron_id = neurons_with_fixed_coordinates[j];
					
					if ( distance_on_sphere(tmp_x, tmp_y, tmp_z, xx_RA[neuron_id], yy_RA[neuron_id], zz_RA[neuron_id]) < too_close_distance )
					{
						too_close = true;
						break;
					}
				}
			}
		}
		
		xx_RA[neurons_to_replace[i]] = tmp_x;
		yy_RA[neurons_to_replace[i]] = tmp_y;
		zz_RA[neurons_to_replace[i]] = tmp_z;
		
		neurons_with_fixed_coordinates.push_back(neurons_to_replace[i]);
	}
}

void NetworkGrowthSimulator::sample_coordinates()
{   
	// check arrangement
	
	double noise_shift_in_interneurons = 0.15; // level of noisyness in ideal pattern of interneurons
	
	if ( topology_params.arrangement == "square" )
	{
		// set coordinates for HVC(I) neurons
		MODEL_INTERNEURON_DISTANCE = 1 / ( sqrt(N_I) + 1 );
	
		std::cout << "Model interneuron distance = " << MODEL_INTERNEURON_DISTANCE << std::endl;
	
		double too_close_distance = MODEL_INTERNEURON_DISTANCE * MIN_INTERSOMATIC_DISTANCE / EXPERIMENTAL_INTERNEURON_DISTANCE;
		
		// HVC(I) neurons form nearly ideal lattice
		
		
		for (int i = 0; i < (int) sqrt(N_I); i++)
		{
			for (int k = 0; k < (int) sqrt(N_I); k++)
			{
				xx_I[i * ((int) sqrt(N_I)) + k] = (double) (i+1) / (sqrt(N_I)+1) - noise_shift_in_interneurons * MODEL_INTERNEURON_DISTANCE + noise_generator.random(2*noise_shift_in_interneurons * MODEL_INTERNEURON_DISTANCE);
				yy_I[i * ((int) sqrt(N_I)) + k] = (double) (k+1) / (sqrt(N_I)+1) - noise_shift_in_interneurons * MODEL_INTERNEURON_DISTANCE + noise_generator.random(2*noise_shift_in_interneurons * MODEL_INTERNEURON_DISTANCE);
			}
		}
		
		// set coordinates for HVC(RA) neurons
		// HVC(RA) neurons are sampled unuformly from the square, but not too close to other neurons
		
		for (int i = 0; i < N_RA; i++)
		{
			bool too_close = true;
		
			double tmp_x, tmp_y; 
		
			while ( too_close )
			{
				too_close = false;
				
				tmp_x = 0.5 / (sqrt(N_I)+1) + noise_generator.random(1.0 - 1.0 / (sqrt(N_I)+1));
				tmp_y = 0.5 / (sqrt(N_I)+1) + noise_generator.random(1.0 - 1.0 / (sqrt(N_I)+1));	
				
				// check if neuron is too close to other interneurons or previously sampled HVC-RA
				
				for (int j = 0; j < N_I; j++)
				{
					if ( distance2d(tmp_x, tmp_y, xx_I[j], yy_I[j]) < too_close_distance )
					{
						too_close = true;
						break;
					}
				}
				
				if ( !too_close )
				{
					for (int j = 0; j < i; j++)
					{
						if ( distance2d(tmp_x, tmp_y, xx_RA[j], yy_RA[j]) < too_close_distance )
						{
							too_close = true;
							break;
						}
					}
				}
			}
			
			xx_RA[i] = tmp_x;
			yy_RA[i] = tmp_y;
		}
	}
		
	if ( topology_params.arrangement == "sphere" )
	{
		// set coordinates for HVC(I) neurons
		// HVC(I) neurons are sampled equally distantly on sphere
		const double pi = 3.14159265358979323846; // value of constant pi
		
		double dphi = pi * (3.0 - sqrt(5.0)); // angle phi increment value
		double dz = 2 / static_cast<double>(N_I); // z increment value
		
		double phi = 0; // current value of angle phi
		double z = 1 - dz/2.0; // current value of z  
		
		double tmp_x1 = cos(phi) * sqrt(1.0 - z*z);
		double tmp_y1 = sin(phi) * sqrt(1.0 - z*z);
		
		double tmp_x2 = cos(phi+dphi) * sqrt(1.0 - (z-dz)*(z-dz));
		double tmp_y2 = sin(phi+dphi) * sqrt(1.0 - (z-dz)*(z-dz));
		
		MODEL_INTERNEURON_DISTANCE = distance_on_sphere(tmp_x1, tmp_y1, z, tmp_x2, tmp_y2, z-dz);
		
		std::cout << "Model interneuron distance = " << MODEL_INTERNEURON_DISTANCE << std::endl;
		
		double too_close_distance = MODEL_INTERNEURON_DISTANCE * MIN_INTERSOMATIC_DISTANCE / EXPERIMENTAL_INTERNEURON_DISTANCE;
		
		for (int i = 0; i < N_I; i++)
		{
			double tmp_x, tmp_y, tmp_z;
		
			bool too_close = true;
			
			// check if neuron is too close to previously sampled interneurons
			while ( too_close )
			{
				too_close = false;
				
				double theta = acos(z);
			
				double theta_noise = -noise_shift_in_interneurons * MODEL_INTERNEURON_DISTANCE + noise_generator.random(2*noise_shift_in_interneurons * MODEL_INTERNEURON_DISTANCE);
				double phi_noise = ( -noise_shift_in_interneurons * MODEL_INTERNEURON_DISTANCE + noise_generator.random(2*noise_shift_in_interneurons * MODEL_INTERNEURON_DISTANCE) ) / sin(theta);
			
				tmp_x = cos(phi + phi_noise) * sin(theta + theta_noise);
				tmp_y = sin(phi + phi_noise) * sin(theta + theta_noise);
				tmp_z = cos(theta + theta_noise);
				
				for (int j = 0; j < i; j++)
				{
					if ( distance_on_sphere(tmp_x, tmp_y, tmp_z, xx_I[j], yy_I[j], zz_I[j]) < too_close_distance )
					{
						too_close = true;
						break;
					}
				}
			}
				
			xx_I[i] = tmp_x;
			yy_I[i] = tmp_y;
			zz_I[i] = tmp_z;
			
			z = z - dz;
			phi = phi + dphi;
		}

		// set coordinates for HVC(RA) neurons
		// HVC(RA) neurons are sample from uniform distribution on sphere
		for (int i = 0; i < N_RA; i++)
		{
			double tmp_x, tmp_y, tmp_z;
			
			bool too_close = true;
			
			// check if neuron is too close to other interneurons or previously sampled HVC-RA
			while ( too_close )
			{
				too_close = false;
			
				z = 2*noise_generator.random(1.0) - 1; // generate z in range(-1, 1)
				phi = noise_generator.random(2.0 * pi); // generate phi in range(0, 2*pi)
			
				tmp_x = sin(acos(z)) * cos(phi);
				tmp_y = sin(acos(z)) * sin(phi);
				tmp_z = z;		
				
				for (int j = 0; j < N_I; j++)
				{
					if ( distance_on_sphere(tmp_x, tmp_y, tmp_z, xx_I[j], yy_I[j], zz_I[j]) < too_close_distance )
					{
						too_close = true;
						break;
					}
				}
				
				if ( !too_close )
				{
					for (int j = 0; j < i; j++)
					{
						if ( distance_on_sphere(tmp_x, tmp_y, tmp_z, xx_RA[j], yy_RA[j], zz_RA[j]) < too_close_distance )
						{
							too_close = true;
							break;
						}
					}
				}
			}
			
			xx_RA[i] = tmp_x;
			yy_RA[i] = tmp_y;
			zz_RA[i] = tmp_z;
		}
	}
}


void NetworkGrowthSimulator::write_different_training_arrangements(std::string networkDir)
{
	///////////////////////////////////
	// Sample clustered training set
	///////////////////////////////////
	
	// find neurons that are close to neuron 0 to form a training set cluster
	// select neuron 0 as the first training neuron and find distances to all other neurons in the pool
	double xx = xx_RA[0];
	double yy = yy_RA[0];
	double zz = zz_RA[0];

	std::vector<double> distances_to_pool_neurons(xx_RA.size()-1); 

	for (size_t i = 1; i < xx_RA.size(); i++)
		distances_to_pool_neurons[i-1] = distance_on_sphere(xx, yy, zz, xx_RA[i], yy_RA[i], zz_RA[i]);

	// sort distances
	std::vector<size_t> idx(distances_to_pool_neurons.size()); // vector with sorted indices of distances

	std::iota(idx.begin(), idx.end(), 0);

	std::sort(idx.begin(), idx.end(), [&distances_to_pool_neurons](size_t i1, size_t i2)
										{return distances_to_pool_neurons[i1] < distances_to_pool_neurons[i2];});

	// select neuron 0 and N_TR - 1 closest neurons to be training set
	
	training_neurons[0] = 0;
	
	for (int i = 1; i < N_TR; i++) 
		training_neurons[i] = idx[i-1] + 1;
		
	this->write_training_neurons((networkDir + "training_neurons_clustered.bin").c_str());
	
	///////////////////////////////////////////////////
	// Sample training neurons randomly from pool
	///////////////////////////////////////////////////
	std::vector<int> unsampled_neurons(N_RA);
	std::iota(unsampled_neurons.begin(), unsampled_neurons.end(), 0);
	
	for (int i = 0; i < N_TR; i++)
	{
		int random_ind = noise_generator.sample_integer(0, N_RA - i - 1);
		
		training_neurons[i] = unsampled_neurons[random_ind];
		unsampled_neurons[random_ind] = unsampled_neurons.back();
		unsampled_neurons.pop_back();
	}
	
	this->write_training_neurons((networkDir + "training_neurons_random.bin").c_str());
	
}

void NetworkGrowthSimulator::set_delays_RA2RA(double delay)
{
	for (int i = 0; i < N_RA; i++)
		std::fill(axonal_delays_RA_RA_global[i].begin(), axonal_delays_RA_RA_global[i].end(), delay);
	
}

void NetworkGrowthSimulator::sample_axonal_delays()
{
	//axonal_delays_RA_RA_global.resize(N_RA);
	//syn_lengths_RA_RA_global.resize(N_RA);
	
	///////////////////////////////////////////////////////////////
	// Sample HVC-RA -> HVC-I  and HVC-RA -> HVC-RA axonal delays
	///////////////////////////////////////////////////////////////
	for (int i = 0; i < N_RA; i++)
	{
		
		for (size_t j = 0; j < syn_ID_RA_I_global[i].size(); j++)
		{
			double delay = syn_lengths_RA_I_global[i][j] * EXPERIMENTAL_AXONAL_DELAY_PER_MICRON * connection_params.delay_constant;
													
			axonal_delays_RA_I_global[i][j] = delay;
			
			//std::cout << "length, delay = " << syn_lengths_RA_I_global[i][j] << ", " << delay << std::endl;
		}
		
		//axonal_delays_RA_RA_global[i].resize(N_RA);
		//syn_lengths_RA_RA_global[i].resize(N_RA);
		
		for (int j = 0; j < N_RA; j++)
		{
			if ( j != i )
				syn_lengths_RA_RA_global[i][j] = distance_on_sphere(xx_RA[i], yy_RA[i], zz_RA[i], 
											  xx_RA[j], yy_RA[j], zz_RA[j]) * EXPERIMENTAL_INTERNEURON_DISTANCE / MODEL_INTERNEURON_DISTANCE;
			else
				syn_lengths_RA_RA_global[i][j] = 0.0; 
										  
			double delay = syn_lengths_RA_RA_global[i][j] * EXPERIMENTAL_AXONAL_DELAY_PER_MICRON * connection_params.delay_constant;
													
			axonal_delays_RA_RA_global[i][j] = delay;
		}
	}
	
	std::cout << "Delay constant = " << connection_params.delay_constant << std::endl;
	
	//~ std::cout << "Axonal delays between HVC-RA neurons:\n";
	//~ for (int i = 0; i < N_RA; i++)
	//~ {
		//~ for (int j = 0; j < N_RA; j++)
			//~ std::cout << axonal_delays_RA_RA_global[i][j] << " ";
		//~ std::cout << "\n";
	//~ }
	//~ std::cout << std::endl;
	
	///////////////////////////////////////////////////////////////
	// Sample HVC-I -> HVC-RA axonal delays
	///////////////////////////////////////////////////////////////
	for (int i = 0; i < N_I; i++)
	{	
		for (size_t j = 0; j < syn_ID_I_RA_global[i].size(); j++)
		{
			double delay = syn_lengths_I_RA_global[i][j] * EXPERIMENTAL_AXONAL_DELAY_PER_MICRON * connection_params.delay_constant;
													
			axonal_delays_I_RA_global[i][j] = delay;
		}
	}
}

void NetworkGrowthSimulator::sample_connectionsAndDelays_for_replaced(std::vector<int>& neurons_to_replace)
{
	// connections for HVC(RA) neurons
	for (size_t i = 0; i < neurons_to_replace.size(); i++)
	{
		int neuron_id = neurons_to_replace[i];
		
		for (int j = 0; j < N_I; j++)
		{
			double d = distance_on_sphere(xx_RA[neuron_id], yy_RA[neuron_id], zz_RA[neuron_id], xx_I[j], yy_I[j], zz_I[j]) * EXPERIMENTAL_INTERNEURON_DISTANCE / MODEL_INTERNEURON_DISTANCE; // distance between HVC(RA) and HVC(I)
		
			if (noise_generator.random(1) < p_RA2I(d))
			{
				double G = this->sample_G(connection_params.Gei_max);

				weights_RA_I_global[neuron_id].push_back(G);
				syn_ID_RA_I_global[neuron_id].push_back(j);
				syn_lengths_RA_I_global[neuron_id].push_back(d);
			
				double delay = d * EXPERIMENTAL_AXONAL_DELAY_PER_MICRON * connection_params.delay_constant;
													
				axonal_delays_RA_I_global[neuron_id].push_back(delay);
			}
		}
	 }
	// connections for HVC(I) neurons

	for (int i = 0; i < N_I; i++)
	{
		for (size_t j = 0; j < neurons_to_replace.size(); j++)
		{
			int neuron_id = neurons_to_replace[j];
			
			double d = distance_on_sphere(xx_I[i], yy_I[i], zz_I[i], xx_RA[neuron_id], yy_RA[neuron_id], zz_RA[neuron_id]) * EXPERIMENTAL_INTERNEURON_DISTANCE / MODEL_INTERNEURON_DISTANCE; // distance between HVC(I) and HVC(RA)
		
			if (noise_generator.random(1) < p_I2RA(d))
			{
				double G = this->sample_G(connection_params.Gie_max);

				weights_I_RA_global[i].push_back(G);
				syn_ID_I_RA_global[i].push_back(neuron_id);
				syn_lengths_I_RA_global[i].push_back(d);
				
				double delay = d * EXPERIMENTAL_AXONAL_DELAY_PER_MICRON * connection_params.delay_constant;
													
				axonal_delays_I_RA_global[i].push_back(delay);
			}
		}
	}

	// update axonal delays for HVC-RA -> HVC-RA connections to and from replaced neurons
	
	// input connections
	for (int i = 0; i < N_RA; i++)
	{
		for (size_t j = 0; j < neurons_to_replace.size(); j++)
		{
			int target_id = neurons_to_replace[j];
			
			if ( target_id != i )
				syn_lengths_RA_RA_global[i][target_id] = distance_on_sphere(xx_RA[i], yy_RA[i], zz_RA[i], 
											  xx_RA[target_id], yy_RA[target_id], zz_RA[target_id]) * EXPERIMENTAL_INTERNEURON_DISTANCE / MODEL_INTERNEURON_DISTANCE;
			else
				syn_lengths_RA_RA_global[i][target_id] = 0.0; 
										  
			double delay = syn_lengths_RA_RA_global[i][target_id] * EXPERIMENTAL_AXONAL_DELAY_PER_MICRON * connection_params.delay_constant;
													
			axonal_delays_RA_RA_global[i][target_id] = delay;
		}
	}
	
	// output connections
	for (size_t i = 0; i < neurons_to_replace.size(); i++)
	{
		for (int j = 0; j < N_RA; j++)
		{
			int source_id = neurons_to_replace[i];
			
			if ( j != source_id )
				syn_lengths_RA_RA_global[source_id][j] = distance_on_sphere(xx_RA[source_id], yy_RA[source_id], zz_RA[source_id], 
											  xx_RA[j], yy_RA[j], zz_RA[j]) * EXPERIMENTAL_INTERNEURON_DISTANCE / MODEL_INTERNEURON_DISTANCE;
			else
				syn_lengths_RA_RA_global[source_id][j] = 0.0; 
										  
			double delay = syn_lengths_RA_RA_global[source_id][j] * EXPERIMENTAL_AXONAL_DELAY_PER_MICRON * connection_params.delay_constant;
													
			axonal_delays_RA_RA_global[source_id][j] = delay;
		}
	}
	
}

void NetworkGrowthSimulator::sample_connections()
{
	if ( topology_params.arrangement == "square" )
	{	
		// connections for HVC(RA) neurons
		for (int i = 0; i < N_RA; i++)
		{
			for (int j = 0; j < N_I; j++)
			{
				double d = distance2d(xx_RA[i], yy_RA[i], xx_I[j], yy_I[j]) * EXPERIMENTAL_INTERNEURON_DISTANCE / MODEL_INTERNEURON_DISTANCE; // distance between HVC(RA) and HVC(I)
			
				if (noise_generator.random(1) < p_RA2I(d))
				{
					double G = this->sample_G(0.0);

					weights_RA_I_global[i].push_back(G);
					syn_ID_RA_I_global[i].push_back(j);
					syn_lengths_RA_I_global[i].push_back(d);
					axonal_delays_RA_I_global[i].push_back(0.0);
				}
			}

		 }
		// connections for HVC(I) neurons

		for (int i = 0; i < N_I; i++)
		{
			for (int j = 0; j < N_RA; j++)
			{
				double d = distance2d(xx_I[i], yy_I[i], xx_RA[j], yy_RA[j]) * EXPERIMENTAL_INTERNEURON_DISTANCE / MODEL_INTERNEURON_DISTANCE; // distance between HVC(I) and HVC(RA)
			
				if (noise_generator.random(1) < p_I2RA(d))
				{
					double G = this->sample_G(0.0);

					weights_I_RA_global[i].push_back(G);
					syn_ID_I_RA_global[i].push_back(j);
					syn_lengths_I_RA_global[i].push_back(d);
					axonal_delays_I_RA_global[i].push_back(0.0);
				}
			}	
		}	
	}
		
	if ( topology_params.arrangement == "sphere" )
	{
		// connections for HVC(RA) neurons
		for (int i = 0; i < N_RA; i++)
		{
			for (int j = 0; j < N_I; j++)
			{
				double d = distance_on_sphere(xx_RA[i], yy_RA[i], zz_RA[i], xx_I[j], yy_I[j], zz_I[j]) * EXPERIMENTAL_INTERNEURON_DISTANCE / MODEL_INTERNEURON_DISTANCE; // distance between HVC(RA) and HVC(I)
			
				if (noise_generator.random(1) < p_RA2I(d))
				{
					double G = this->sample_G(0.0);

					weights_RA_I_global[i].push_back(G);
					syn_ID_RA_I_global[i].push_back(j);
					syn_lengths_RA_I_global[i].push_back(d);
					axonal_delays_RA_I_global[i].push_back(0.0);
				}
			}

		 }
		// connections for HVC(I) neurons

		for (int i = 0; i < N_I; i++)
		{
			for (int j = 0; j < N_RA; j++)
			{
				double d = distance_on_sphere(xx_I[i], yy_I[i], zz_I[i], xx_RA[j], yy_RA[j], zz_RA[j]) * EXPERIMENTAL_INTERNEURON_DISTANCE / MODEL_INTERNEURON_DISTANCE; // distance between HVC(I) and HVC(RA)
			
				if (noise_generator.random(1) < p_I2RA(d))
				{
					double G = this->sample_G(0.0);

					weights_I_RA_global[i].push_back(G);
					syn_ID_I_RA_global[i].push_back(j);
					syn_lengths_I_RA_global[i].push_back(d);
					axonal_delays_I_RA_global[i].push_back(0.0);
				}
			}
		}
	}
}

double NetworkGrowthSimulator::p_RA2I(double d)
{
	double sigma_square = topology_params.SIGMA_RA2I * topology_params.SIGMA_RA2I;
	
	return topology_params.A_RA2I * exp(-d*d / (2 * sigma_square));
}

double NetworkGrowthSimulator::p_I2RA(double d)
{
	double sigma_square = topology_params.SIGMA_I2RA * topology_params.SIGMA_I2RA;
	
	return topology_params.B_I2RA * exp(-d*d / (2 * sigma_square));		
}

void NetworkGrowthSimulator::resample_weights()
{
	if (MPI_rank == 0)
	{
		// Resample HVC-I -> HVC-RA connections
		for (size_t i = 0; i < weights_I_RA_global.size(); i++)
			for (size_t j = 0; j < weights_I_RA_global[i].size(); j++)
				weights_I_RA_global[i][j] = this->sample_G(connection_params.Gie_max);
		
		// Resample HVC-RA -> HVC-I connections
		for (size_t i = 0; i < weights_RA_I_global.size(); i++)
			for (size_t j = 0; j < weights_RA_I_global[i].size(); j++)
				weights_RA_I_global[i][j] = this->sample_G(connection_params.Gei_max);
		
	}
	
}

double NetworkGrowthSimulator::sample_G(double G_max)
{
    return noise_generator.random(1.0) * G_max;
}


void NetworkGrowthSimulator::initialize_network()
{
	if (MPI_rank == 0)
	{
		// get number of neurons in the network
		std::cout << "N_RA = " << N_RA << std::endl;
		std::cout << "N_TR = " << N_TR << std::endl;
		std::cout << "N_I = " << N_I << std::endl;
		
	}

	// send number of neurons to all processes
	MPI_Bcast(&N_RA, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&N_TR, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&N_I, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if ( MPI_rank != 0 )
		training_neurons.resize(N_TR);
	
	// get connections between HVC(RA) and HVC(I) to master process
	if (MPI_rank == 0)
	{
		std::cout << "Training neurons: " << std::endl;
		for (size_t i = 0; i < training_neurons.size(); i++)
			std::cout << training_neurons[i] << "\t";
		std::cout << std::endl;
	}
	
	//this->print_invariable_connections();
	
	// send training neurons to all processes
	MPI_Bcast(&training_neurons[0], N_TR, MPI_INT, 0, MPI_COMM_WORLD);
	
	// send training spread to all processes
	training_spread_times.resize(N_TR);
	MPI_Bcast(&training_spread_times[0], N_TR, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	
	
	// distibute work among processes
	N_RA_sizes.resize(MPI_size); // array with number of RA neurons per process
	N_I_sizes.resize(MPI_size); // array with number of I neurons per process

	for (int i = 0; i < MPI_size; i++)
	{
		N_RA_sizes[i] = N_RA / MPI_size;
		N_I_sizes[i] = N_I / MPI_size;
	}
	int RA_remain = N_RA % MPI_size;
	int I_remain = N_I % MPI_size;
	int j = 0;

	// distribute RA neurons
	while (RA_remain > 0)
	{
		N_RA_sizes[j] += 1;
		RA_remain -= 1;
		j += 1;

		if (j >= MPI_size)
			j -= MPI_size;
	}

	// distribute I neurons
	j = 0;
	while (I_remain > 0)
	{
		
		N_I_sizes[j] += 1;
		I_remain -= 1;
		j += 1;

		if (j >= MPI_size)
			j -= MPI_size;
	}

	N_RA_local = N_RA_sizes[MPI_rank]; // assign number of RA neurons for each process
	N_I_local = N_I_sizes[MPI_rank]; // assign number of I neurons for each process

	printf("My rank is %d; N_RA_local = %d; N_I_local = %d\n", MPI_rank, N_RA_local, N_I_local);

    // resize all data arrays
    this->resize_arrays_for_all_processes();
	
    // assign real id to neurons
	for (int i = 0; i < N_RA_local; i++)
	{
		// assign real Id for RA neurons
		int N = 0; // number of neurons in the processes with lower rank
		
		for (int k = 0; k < MPI_rank; k++)
			N += N_RA_sizes[k];

		Id_RA_local[i] = N + i;
	}
	
    // assign real ID for I neurons
    for (int i = 0; i < N_I_local; i++)
	{
		int N = 0; // number of neurons in the processes with lower rank
		
		for (int k = 0; k < MPI_rank; k++)
			N += N_I_sizes[k];

		Id_I_local[i] = N + i;
	}

    // prepare global array with ids of all HVC(RA) neurons
    Id_RA_global.resize(N_RA); 
    

    for (int i = 0; i < N_RA; i++)
        Id_RA_global[i] = i;
	
	
	// make training neurons mature on master process
	//if (MPI_rank == 0)
	//{
		//for (size_t i = 0; i < training_neurons.size(); i++)
			//mature_global[training_neurons[i]] = 1;
	//}
	

	
	this->set_noise();
	this->set_dynamics();
	
	this->send_connections_RAandI();
	this->send_connections_RA2RA();
	
	this->send_active_synapses();
	this->send_super_synapses();
	this->send_activity_history();
	this->send_replacement_history();
	this->send_remodeled_indicators();
	this->send_maturation_properties();
	//this->send_mature_indicators();
	
	
	this->set_synapse_indicators();
	
	// set neuron model parameters
	//this->update_neuron_properties();
	//this->set_neuron_properties_sudden_maturation();
	this->set_neuron_properties();
}

void NetworkGrowthSimulator::test_inhAndExc_response(const struct ConnectionParameters &con_par,
													const struct NoiseParameters &n_par,
													std::string networkDirectory, std::string fileTraining,
													std::string outputDirectory)
{
	if (MPI_rank == 0)
	{
		connection_params = con_par;
		noise_params = n_par;
	}
	
	MPI_Bcast(&connection_params.Gei_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&connection_params.Gie_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&connection_params.delay_constant, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	MPI_Bcast(&noise_params.white_noise_mean_soma, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&noise_params.white_noise_std_soma, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&noise_params.white_noise_mean_dend, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&noise_params.white_noise_std_dend, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if (MPI_rank == 0)
	{
		// read initial network from directory
		this->read_network_topology(networkDirectory);
		
		// write topology parameters to output directory
		ConfigurationNetworkTopology cfg_topology;
		cfg_topology.set_topology_parameters(topology_params);

		this->sample_axonal_delays();
		this->resample_weights();
		
		this->read_training_neurons(fileTraining.c_str());
		
		// sort training neurons for convenience
		std::sort(training_neurons.begin(), training_neurons.end());
		
		this->create_second_layer();
	}
	
	this->initialize_network();
	
	//this->set_all_neurons_immature();
	//this->set_training_neurons_mature();
	
	// rescale inhibitory synapses to training neurons
	for (size_t i = 0; i < training_neurons.size(); i++)
		this->rescale_synapses_to_mature(training_neurons[i]);
		
	this->disable_RA2I_immature_outputs();
	
	
	double window = 25.0;
	int num_trials = 100;
	
	std::vector<std::vector<double>> spike_times;
	std::vector<std::vector<double>> relevant_spike_times;
	
	std::vector<std::vector<int>> num_somatic_spikes_in_trials;
	std::vector<std::vector<int>> num_relevant_spikes_in_trials;
	
	std::vector<std::vector<double>> Ginh;
	
	// set training kick time
	double training_kick_time = 100;
	
	if (MPI_rank == 0)
	{
		spike_times.resize(N_RA);
		relevant_spike_times.resize(N_RA);
		
		num_somatic_spikes_in_trials.resize(N_RA);
		num_relevant_spikes_in_trials.resize(N_RA);
		
		Ginh.resize(N_RA);
		
		for (int i = 0; i < N_RA; i++)
		{
			num_somatic_spikes_in_trials[i].resize(num_trials);
			num_relevant_spikes_in_trials[i].resize(num_trials);
			//Ginh[i].resize(static_cast<int>(TRIAL_DURATION / TIMESTEP));
		}
	}
	
	int num_steps = 5;
	double dG = 10.0;
	double G = 0.0;
	
	for (int t = 0; t < num_steps; t++){
	
		this->set_active_synapse_weights(G);
		
		for (int i = 0; i < num_trials; i++)
		{
			if (MPI_rank == 0)
				std::cout << "Trial " << i << std::endl;
			
			// set recording for the last trial
			//if ( i == num_trials - 1)
				//this->set_recording(RAtoWrite, ItoWrite, outputDirectory);
		
			
			//this->trial_no_stdp_inhAndExc(training_kick_time, &Ginh);
			this->trial_no_stdp(training_kick_time);
			
			
			int bad_values_indicator_local = this->check_bad_values();
			int bad_values_indicator_global;
			
			MPI_Allreduce(&bad_values_indicator_local, &bad_values_indicator_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		   
			if (bad_values_indicator_global < 0)
			{
				//std::endl << ""
				return;
			}  
			
		  
			this->gather_graph_state_data();

			if (MPI_rank == 0)
			{
				// spikes and bursts of HVC-RA neurons
				
				size_t num_HVCRA_spikes = 0; // total number of spikes of HVC-RA neurons during trial
				size_t num_HVCRA_silent = 0; // number of HVC-RA neurons that did not burst
				
				for (int j = 0; j < N_RA; j++)
				{
					if ( spikes_in_trial_soma_global[j].size() >= 1 ){
						std::copy(spikes_in_trial_soma_global[j].begin(), spikes_in_trial_soma_global[j].end(),
										std::back_inserter(spike_times[j]));
						
						auto it_up = std::upper_bound(spikes_in_trial_soma_global[j].begin(), spikes_in_trial_soma_global[j].end(),
																	training_kick_time + window);
						
						auto it_bottom = std::lower_bound(spikes_in_trial_soma_global[j].begin(), spikes_in_trial_soma_global[j].end(),
																	training_kick_time);
																	
						std::copy(it_bottom, it_up, std::back_inserter(relevant_spike_times[j]));
						
						num_relevant_spikes_in_trials[j][i] = std::distance(it_bottom, it_up);
						
						//printf("Average dendritic spike time = %f\n", average_spike_time);
						
					}
					else
						num_HVCRA_silent += 1;
						
					
					num_HVCRA_spikes += spikes_in_trial_soma_global[j].size();
					
					
					num_somatic_spikes_in_trials[j][i] = static_cast<int>(spikes_in_trial_soma_global[j].size());
			
				}
				
				std::cout << "Number of silent HVC-RA = " << num_HVCRA_silent << "\n" << std::endl;
				
				std::cout << "Average number of spikes of HVC-RA = " << static_cast<double>(num_HVCRA_spikes) / static_cast<double>(N_RA) << "\n" << std::endl;
				
				std::cout << "Average spike frequency of HVC-RA = " << static_cast<double>(num_HVCRA_spikes) * 1000.0 / (TRIAL_DURATION * static_cast<double>(N_RA)) << "\n" << std::endl;
				
				// spikes of HVC-I neurons
				size_t num_HVCI_spikes = 0; // total number of spikes of HVC-I neurons during trial
				
				for (int j = 0; j < N_I; j++)
					num_HVCI_spikes += spikes_in_trial_interneuron_global[j].size();
					
				std::cout << "Average number of spikes of HVC-I = " << static_cast<double>(num_HVCI_spikes) / static_cast<double>(N_I) << "\n" << std::endl;
				std::cout << "Average spike frequency of HVC-I = " << static_cast<double>(num_HVCI_spikes) * 1000.0 / (TRIAL_DURATION * static_cast<double>(N_I)) << "\n" << std::endl;
			
			}

			this->reset_after_chain_test();
		}
		
		if ( MPI_rank == 0 ){
			std::stringstream sstream;
			sstream << std::fixed << std::setprecision(2) << "Ginh_" << connection_params.Gie_max << "_Gee_" << G << ".bin";
	
		
			std::string filename = outputDirectory + sstream.str();
		
		
			this->calculate_and_write_inhAndExc(spike_times, relevant_spike_times,
									num_somatic_spikes_in_trials, num_relevant_spikes_in_trials,
									filename.c_str());
		
											
			for (int i = 0; i < N_RA; i++){
				std::fill(num_somatic_spikes_in_trials[i].begin(), num_somatic_spikes_in_trials[i].end(), 0.0);
				std::fill(num_relevant_spikes_in_trials[i].begin(), num_relevant_spikes_in_trials[i].end(), 0.0);
				
				spike_times[i].clear();
				relevant_spike_times[i].clear();
			}
		}
		// change excitatory weight
		G += dG;
		
		
	}
	
}

void NetworkGrowthSimulator::create_second_layer()
{
	// connect training neurons to all immature pool neurons
	for (size_t i = 0; i < training_neurons.size(); i++){
		for (int j = 0; j < N_RA; j++){
			auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), j);
									
			if ( p.first == p.second )
				active_synapses_global[training_neurons[i]].push_back(j);
		}	
	}
	
	
	
}

void NetworkGrowthSimulator::test_synthetic_chain(const ConfigurationNetworkGrowth &cfg, int num_layers, int num_group, double probability, double Gee,
										std::string networkDirectory, int num_trials, std::string outputDirectory)
{
	if (MPI_rank == 0)
	{
		connection_params = cfg.get_connection_parameters();
		noise_params = cfg.get_noise_parameters();
		maturation_params = cfg.get_maturation_parameters(); 
	}
	
	this->send_growth_parameters();
	
	if (MPI_rank == 0)
	{
		// read initial network from directory
		this->read_network_topology(networkDirectory);
		
		// write topology parameters to output directory
		ConfigurationNetworkTopology cfg_topology;
		cfg_topology.set_topology_parameters(topology_params);
		cfg_topology.write_configuration((outputDirectory + "network_parameters.cfg").c_str());
		
		this->sample_axonal_delays();
		this->resample_weights();
		
		this->set_delays_RA2RA(3.0);
		
		this->generate_synthetic_chain(num_layers, num_group, probability, Gee);
		
		// sort training neurons for convenience
		std::sort(training_neurons.begin(), training_neurons.end());
		
		// make all chain neurons mature
		for (int i = 0; i < num_layers; i++)
		{
			for (int j = 0; j < num_group; j++)
			{
				mature_global[i*num_group+j] = 1;
				GCa_global[i*num_group+j] = maturation_params.GCA_MATURE;
				rest_potential_global[i*num_group+j] = maturation_params.E_REST_MATURE;
			}
		}
		
		// write network to output directory both as initial
		std::string filename_I_coordinates = outputDirectory + "I_xy.bin";
		std::string filename_out_training = outputDirectory + "training_neurons.bin";
		std::string filename_training_spread = outputDirectory + "training_spread.bin";
		std::string filename_num_neurons = outputDirectory + "num_neurons.bin";
		
		// generate fixed spread times
		double spread = 5.0;
		
		training_spread_times.resize(N_TR);
		for (int i = 0; i < N_TR; i++)
			training_spread_times[i] = noise_generator.random(spread) - spread / 2.;
	
	
		this->write_coordinates(xx_I,  yy_I,  zz_I,  filename_I_coordinates.c_str());
		this->write_training_neurons(filename_out_training.c_str());
		this->write_training_spread(filename_training_spread.c_str());
		this->write_number_of_neurons(filename_num_neurons.c_str());
		this->write_alterable_network("_initial", outputDirectory);
		
		
		cfg.write_configuration((outputDirectory + "growth_parameters.cfg").c_str());
		
		std::string fileTraining = outputDirectory + "training_neurons.bin";
		
		this->write_full_network_state("_" + std::to_string(0), outputDirectory);
		this->write_training_neurons(fileTraining.c_str());
	}
	
	this->initialize_network();
	
		
	// start simulation
	std::vector<std::vector<double>> first_soma_spike_times;
	std::vector<std::vector<double>> first_dend_spike_times;
	
	std::vector<std::vector<int>> num_somatic_spikes_in_trials;
	std::vector<std::vector<int>> num_dendritic_spikes_in_trials;
	
	
	std::vector<double> spread_times(N_RA);
	
	if (MPI_rank == 0){
		for (int i = 0; i < N_TR; i++){
			spread_times[training_neurons[i]] = training_spread_times[i];
			std::cout << "Neuron " << training_neurons[i] << " spread time = " << spread_times[training_neurons[i]] << std::endl;
		}
	}
	
	MPI_Bcast(&spread_times[0], N_RA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&training_spread_times[0], N_TR, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if (MPI_rank == 0)
	{
		first_soma_spike_times.resize(N_RA);
		first_dend_spike_times.resize(N_RA);
		
		num_somatic_spikes_in_trials.resize(N_RA);
		num_dendritic_spikes_in_trials.resize(N_RA);
		
		for (int i = 0; i < N_RA; i++)
		{
			num_somatic_spikes_in_trials[i].resize(num_trials);
			num_dendritic_spikes_in_trials[i].resize(num_trials);
		}
	}
	
    for (int i = 0; i < num_trials; i++)
    {
        if (MPI_rank == 0)
            std::cout << "Trial " << i << std::endl;
		
		// set recording for the last trial
		//if ( i == num_trials - 1)
			//this->set_recording(RAtoWrite, ItoWrite, outputDirectory);
	
		
		//this->trial_no_stdp(training_kick_time);
	    this->trial_no_stdp_fixedSpread(spread_times);
	    
	    int bad_values_indicator_local = this->check_bad_values();
	    int bad_values_indicator_global;
	    
	    MPI_Allreduce(&bad_values_indicator_local, &bad_values_indicator_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	   
		if (bad_values_indicator_global < 0)
		{
			//std::endl << ""
			return;
        }  
	    
      
		this->gather_graph_state_data();

		this->write_chain_test(outputDirectory, i);
		
		if (MPI_rank == 0)
		{
			// spikes and bursts of HVC-RA neurons
			
			size_t num_HVCRA_spikes = 0; // total number of spikes of HVC-RA neurons during trial
			size_t num_HVCRA_bursts = 0; // total number of bursts of HVC-RA neurons during trial
			size_t num_HVCRA_silent = 0; // number of HVC-RA neurons that did not burst
			
			for (int j = 0; j < N_RA; j++)
			{
				if ( spikes_in_trial_soma_global[j].size() >= 1 )
					//printf("Average dendritic spike time = %f\n", average_spike_time);
					first_soma_spike_times[j].push_back(spikes_in_trial_soma_global[j][0]);
				else
					num_HVCRA_silent += 1;
					
				if ( spikes_in_trial_dend_global[j].size() >= 1 )
					//printf("Average dendritic spike time = %f\n", average_spike_time);
					first_dend_spike_times[j].push_back(spikes_in_trial_dend_global[j][0]);
				
				num_HVCRA_spikes += spikes_in_trial_soma_global[j].size();
				num_HVCRA_bursts += spikes_in_trial_dend_global[j].size();
				
				
                num_dendritic_spikes_in_trials[j][i] = static_cast<int>(spikes_in_trial_dend_global[j].size());
                num_somatic_spikes_in_trials[j][i] = static_cast<int>(spikes_in_trial_soma_global[j].size());
        
			}
			
			std::cout << "Number of silent HVC-RA = " << num_HVCRA_silent << "\n" << std::endl;
			
			std::cout << "Average number of bursts of HVC-RA = " << static_cast<double>(num_HVCRA_bursts) / static_cast<double>(N_RA) << "\n" << std::endl;
			std::cout << "Average number of spikes of HVC-RA = " << static_cast<double>(num_HVCRA_spikes) / static_cast<double>(N_RA) << "\n" << std::endl;
			
			std::cout << "Average burst frequency of HVC-RA = " << static_cast<double>(num_HVCRA_bursts) * 1000.0 / (TRIAL_DURATION * static_cast<double>(N_RA)) << "\n" << std::endl;
			std::cout << "Average spike frequency of HVC-RA = " << static_cast<double>(num_HVCRA_spikes) * 1000.0 / (TRIAL_DURATION * static_cast<double>(N_RA)) << "\n" << std::endl;
			
			// spikes of HVC-I neurons
			size_t num_HVCI_spikes = 0; // total number of spikes of HVC-I neurons during trial
			
			for (int j = 0; j < N_I; j++)
				num_HVCI_spikes += spikes_in_trial_interneuron_global[j].size();
				
			std::cout << "Average number of spikes of HVC-I = " << static_cast<double>(num_HVCI_spikes) / static_cast<double>(N_I) << "\n" << std::endl;
			std::cout << "Average spike frequency of HVC-I = " << static_cast<double>(num_HVCI_spikes) * 1000.0 / (TRIAL_DURATION * static_cast<double>(N_I)) << "\n" << std::endl;
		
		}

		this->reset_after_chain_test();
	}
	
	// process dendritic spikes
	if (MPI_rank == 0)
		this->calculate_and_write_jitter(num_trials, first_soma_spike_times, first_dend_spike_times,
									num_somatic_spikes_in_trials, num_dendritic_spikes_in_trials,
											 outputDirectory);                      
}



void NetworkGrowthSimulator::test_chain(std::string networkDirectory,  int starting_trial, int num_trials, std::string outputDirectory)
{
	trial_number = starting_trial;
	
	if (MPI_rank == 0)
		this->read_network_state(networkDirectory, starting_trial);
    
    
    this->send_growth_parameters();
    this->initialize_network();
    
    
    //this->set_time_for_neurons(trial_number * trial_duration);
    //this->print_super();
    
    // rescale inhibition to mature neurons because rescaled synapses were not saved
    // in global arrays
    //~ for (int i = 0; i < N_RA; i++){
		//~ if ( mature_global[i] == 1 )
			//~ this->rescale_inhibition_to_mature(i);
    //~ }
    //~ 
    this->gather_graph_state_data();
    this->gather_full_state_data();
    
    //if (MPI_rank == 0)
		//this->write_full_network_state("_" + std::to_string(trial_number) + "afterReading", outputDirectory);
	
	
	// start simulations
	
	std::vector<int> RAtoWrite{};
    //std::vector<int> ItoWrite{0, 1, 2, 3};
    std::vector<int> ItoWrite{};
    

    //std::vector<int> source{0, 1, 2, 3};
    
    std::vector<int> source{};
    std::vector<int> target{};

    for (int i = 0; i < N_RA; i++)
        RAtoWrite.push_back(i);

    // make all neurons able to make output connections
    //this->set_all_mature();

	
	std::vector<std::vector<double>> first_soma_spike_times;
	std::vector<std::vector<double>> first_dend_spike_times;
	
	std::vector<std::vector<int>> num_somatic_spikes_in_trials;
	std::vector<std::vector<int>> num_dendritic_spikes_in_trials;
	
	
	std::vector<double> spread_times(N_RA);
	
	if (MPI_rank == 0){
		for (int i = 0; i < N_TR; i++){
			spread_times[training_neurons[i]] = training_spread_times[i];
			std::cout << "Neuron " << training_neurons[i] << " spread time = " << spread_times[training_neurons[i]] << std::endl;
		}
	}
	
	MPI_Bcast(&spread_times[0], N_RA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&training_spread_times[0], N_TR, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if (MPI_rank == 0)
	{
		first_soma_spike_times.resize(N_RA);
		first_dend_spike_times.resize(N_RA);
		
		num_somatic_spikes_in_trials.resize(N_RA);
		num_dendritic_spikes_in_trials.resize(N_RA);
		
		for (int i = 0; i < N_RA; i++)
		{
			num_somatic_spikes_in_trials[i].resize(num_trials);
			num_dendritic_spikes_in_trials[i].resize(num_trials);
		}
	}
	
    for (int i = 0; i < num_trials; i++)
    {
        if (MPI_rank == 0)
            std::cout << "Trial " << i << std::endl;
		
		// set recording for the last trial
		
		this->set_recording(RAtoWrite, ItoWrite, i, outputDirectory);
	
		
		//this->trial_no_stdp(training_kick_time);
	    this->trial_no_stdp_fixedSpread(spread_times);
	    
	    int bad_values_indicator_local = this->check_bad_values();
	    int bad_values_indicator_global;
	    
	    MPI_Allreduce(&bad_values_indicator_local, &bad_values_indicator_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	   
		if (bad_values_indicator_global < 0)
		{
			//std::endl << ""
			return;
        }  
	    
      
		this->gather_graph_state_data();

		this->write_chain_test(outputDirectory, i);
		
		if (MPI_rank == 0)
		{
			// spikes and bursts of HVC-RA neurons
			
			size_t num_HVCRA_spikes = 0; // total number of spikes of HVC-RA neurons during trial
			size_t num_HVCRA_bursts = 0; // total number of bursts of HVC-RA neurons during trial
			size_t num_HVCRA_silent = 0; // number of HVC-RA neurons that did not burst
			
			for (int j = 0; j < N_RA; j++)
			{
				if ( spikes_in_trial_soma_global[j].size() >= 1 )
					//printf("Average dendritic spike time = %f\n", average_spike_time);
					first_soma_spike_times[j].push_back(spikes_in_trial_soma_global[j][0]);
				else
					num_HVCRA_silent += 1;
					
				if ( spikes_in_trial_dend_global[j].size() >= 1 )
					//printf("Average dendritic spike time = %f\n", average_spike_time);
					first_dend_spike_times[j].push_back(spikes_in_trial_dend_global[j][0]);
				
				num_HVCRA_spikes += spikes_in_trial_soma_global[j].size();
				num_HVCRA_bursts += spikes_in_trial_dend_global[j].size();
				
				
                num_dendritic_spikes_in_trials[j][i] = static_cast<int>(spikes_in_trial_dend_global[j].size());
                num_somatic_spikes_in_trials[j][i] = static_cast<int>(spikes_in_trial_soma_global[j].size());
        
			}
			
			std::cout << "Number of silent HVC-RA = " << num_HVCRA_silent << "\n" << std::endl;
			
			std::cout << "Average number of bursts of HVC-RA = " << static_cast<double>(num_HVCRA_bursts) / static_cast<double>(N_RA) << "\n" << std::endl;
			std::cout << "Average number of spikes of HVC-RA = " << static_cast<double>(num_HVCRA_spikes) / static_cast<double>(N_RA) << "\n" << std::endl;
			
			std::cout << "Average burst frequency of HVC-RA = " << static_cast<double>(num_HVCRA_bursts) * 1000.0 / (TRIAL_DURATION * static_cast<double>(N_RA)) << "\n" << std::endl;
			std::cout << "Average spike frequency of HVC-RA = " << static_cast<double>(num_HVCRA_spikes) * 1000.0 / (TRIAL_DURATION * static_cast<double>(N_RA)) << "\n" << std::endl;
			
			// spikes of HVC-I neurons
			size_t num_HVCI_spikes = 0; // total number of spikes of HVC-I neurons during trial
			
			for (int j = 0; j < N_I; j++)
				num_HVCI_spikes += spikes_in_trial_interneuron_global[j].size();
				
			std::cout << "Average number of spikes of HVC-I = " << static_cast<double>(num_HVCI_spikes) / static_cast<double>(N_I) << "\n" << std::endl;
			std::cout << "Average spike frequency of HVC-I = " << static_cast<double>(num_HVCI_spikes) * 1000.0 / (TRIAL_DURATION * static_cast<double>(N_I)) << "\n" << std::endl;
		
		}

		this->reset_after_chain_test();
	}
	
	// process dendritic spikes
	if (MPI_rank == 0)
		this->calculate_and_write_jitter(num_trials, first_soma_spike_times, first_dend_spike_times,
									num_somatic_spikes_in_trials, num_dendritic_spikes_in_trials,
											 outputDirectory);
                           
	
}

void NetworkGrowthSimulator::calculate_and_write_inhAndExc(const std::vector<std::vector<double>> &spike_times,
											  const std::vector<std::vector<double>> &relevant_spike_times,
											  const std::vector<std::vector<int>> &num_somatic_spikes_in_trials,
											  const std::vector<std::vector<int>> &num_relevant_spikes_in_trials,
											  const char *filename)
{
	int num_trials;
	
	if ( !num_somatic_spikes_in_trials.empty() )
		num_trials = num_somatic_spikes_in_trials[0].size();
	else{
		std::cerr << "Empty array in calculate_and_write_inhAndExc!\n" << std::endl;
		return;
	}
	
	std::vector<double> mean_relevant_spike_time; // average of first somatic spike time
	std::vector<double> std_relevant_spike_time; // standard deviation of first somatic spike time
	std::vector<double> mean_spike_time; // average of first dendritic spike time
	std::vector<double> std_spike_time; // standard deviation of first dendritic spike time
	std::vector<double> average_num_soma_spikes_in_trial; // average number of somatic spikes in trials
	std::vector<int> num_trials_with_soma_spikes; // number of trials in which neuron produced somatic spikes
	std::vector<int> num_trials_with_relevant_spikes; // number of trials in which neuron produced relevant somatic spikes
	std::vector<int> num_trials_with_nonrelevant_spikes; // number of trials in which neuron produced relevant somatic spikes

	mean_relevant_spike_time.resize(N_RA);
	std_relevant_spike_time.resize(N_RA);
	mean_spike_time.resize(N_RA);
	std_spike_time.resize(N_RA);
	average_num_soma_spikes_in_trial.resize(N_RA);
	num_trials_with_soma_spikes.resize(N_RA);
	num_trials_with_relevant_spikes.resize(N_RA);
	num_trials_with_nonrelevant_spikes.resize(N_RA);
	
	
	for (int i = 0; i < N_RA; i++)
	{
		average_num_soma_spikes_in_trial[i] = std::accumulate(num_somatic_spikes_in_trials[i].begin(), num_somatic_spikes_in_trials[i].end(), 0.0) 
											/ static_cast<double>(num_trials);
		
		for (int j = 0; j < num_trials; j++){
			if ( num_somatic_spikes_in_trials[i][j] > 0 )
				num_trials_with_soma_spikes[i] += 1;
				
			if ( num_relevant_spikes_in_trials[i][j] > 0 )
				num_trials_with_relevant_spikes[i] += 1;
				
			if ( num_somatic_spikes_in_trials[i][j] > num_relevant_spikes_in_trials[i][j] )
				num_trials_with_nonrelevant_spikes[i] += 1;
		}
		
		///////////////////////////
		// Mean somatic spike times
		///////////////////////////////
		if ( spike_times[i].size() > 0)
		{

			//for (int j = 0; j < (int) average_dendritic_spike_time[i].size(); j++)
			  //  printf("average_dendritic_spike_time[%d][%d] = %f\n", i, j, average_dendritic_spike_time[i][j]);
			mean_spike_time[i] = std::accumulate(spike_times[i].begin(), spike_times[i].end(), 0.0) / 
									static_cast<double>(spike_times[i].size());

		}

		else
			mean_spike_time[i] = -1;
		// calculate standard deviation of spike times

		double accum = 0;
		double mean = mean_spike_time[i];

		std::for_each(spike_times[i].begin(), spike_times[i].end(), [&accum, mean](const double t)
		{
			accum += (t - mean) * (t - mean);
		});

		if (static_cast<int>(spike_times[i].size() > 1))
			std_spike_time[i] = sqrt(accum / (static_cast<double>(spike_times[i].size()) - 1));
		else
			std_spike_time[i] = -1;
		
		
		///////////////////////////
		// Relevant somatic spike times
		///////////////////////////////
		if ( relevant_spike_times[i].size() > 0)
		{

			//for (int j = 0; j < (int) average_dendritic_spike_time[i].size(); j++)
			  //  printf("average_dendritic_spike_time[%d][%d] = %f\n", i, j, average_dendritic_spike_time[i][j]);
			mean_relevant_spike_time[i] = std::accumulate(relevant_spike_times[i].begin(), relevant_spike_times[i].end(), 0.0) / 
									static_cast<double>(relevant_spike_times[i].size());

		}

		else
			mean_relevant_spike_time[i] = -1;
		// calculate standard deviation of spike times

		accum = 0;
		mean = mean_relevant_spike_time[i];

		std::for_each(relevant_spike_times[i].begin(), relevant_spike_times[i].end(), [&accum, mean](const double t)
		{
			accum += (t - mean) * (t - mean);
		});

		if (static_cast<int>(relevant_spike_times[i].size() > 1))
			std_relevant_spike_time[i] = sqrt(accum / (static_cast<double>(relevant_spike_times[i].size()) - 1));
		else
			std_relevant_spike_time[i] = -1;
		
	}

	std::cout << "Mean relevant spike times: \n";
	for (int i = 0; i < N_RA; i++){
		std::cout << mean_relevant_spike_time[i] << " ";
	}
	std::cout << "\n" << std::endl;
	
	std::cout << "Probability for relevant spike: \n";
	for (int i = 0; i < N_RA; i++){
		std::cout << static_cast<double>(num_trials_with_relevant_spikes[i]) / static_cast<double>(num_trials) << " ";
	}
	std::cout << "\n" << std::endl;
	
	std::cout << "Std relevant spike times: \n";
	for (int i = 0; i < N_RA; i++){
		std::cout << std_relevant_spike_time[i] << " ";
	}
	std::cout << "\n" << std::endl;
	
	
	std::cout << "Mean spike times: \n";
	for (int i = 0; i < N_RA; i++){
		std::cout << mean_spike_time[i] << " ";
	}
	std::cout << "\n" << std::endl;
	
	std::cout << "Probability for non-relevant spike in trial: \n";
	for (int i = 0; i < N_RA; i++){
		std::cout << static_cast<double>(num_trials_with_nonrelevant_spikes[i]) / static_cast<double>(num_trials) << " ";
	}
	std::cout << "\n" << std::endl;
	
	std::cout << "Std spike times: \n";
	for (int i = 0; i < N_RA; i++){
		std::cout << std_spike_time[i] << " ";
	}
	std::cout << "\n" << std::endl;	
	
	
	this->write_inhAndExc(num_trials, num_trials_with_relevant_spikes, mean_relevant_spike_time, 
							std_relevant_spike_time, num_trials_with_nonrelevant_spikes, mean_spike_time, std_spike_time,
							 filename);
}


void NetworkGrowthSimulator::calculate_and_write_jitter(int num_trials,
												const std::vector<std::vector<double>> &first_soma_spike_times,
											  const std::vector<std::vector<double>> &first_dend_spike_times,
											  const std::vector<std::vector<int>> &num_somatic_spikes_in_trials,
											  const std::vector<std::vector<int>> &num_dendritic_spikes_in_trials,
											  std::string outputDirectory)
{
	//~ std::cout << "Somatic spikes of HVC-RA neurons\n";
	//~ for (int i = 0; i < N_RA; i++)
	//~ {
		//~ std::cout << "HVC-RA " << i << " ";
		//~ for (size_t j = 0; j < first_soma_spike_times[i].size(); j++)
			//~ std::cout << first_soma_spike_times[i][j] << " ";
		//~ std::cout << std::endl;
	//~ }
	//~ 
	//~ std::cout << "\nNumber of somatic spikes of HVC-RA neurons\n";
	//~ for (int i = 0; i < N_RA; i++)
	//~ {
		//~ std::cout << "HVC-RA " << i << " ";
		//~ for (size_t j = 0; j < num_somatic_spikes_in_trials[i].size(); j++)
			//~ std::cout << num_somatic_spikes_in_trials[i][j] << " ";
		//~ std::cout << std::endl;
		//~ 
	//~ }
	//~ std::cout << std::endl;
	
	std::vector<double> mean_first_soma_spike_time; // average of first somatic spike time
	std::vector<double> std_first_soma_spike_time; // standard deviation of first somatic spike time
	std::vector<double> mean_first_dend_spike_time; // average of first dendritic spike time
	std::vector<double> std_first_dend_spike_time; // standard deviation of first dendritic spike time
    std::vector<double> average_num_dend_spikes_in_trial; // average number of dendritic spikes in trial
    std::vector<double> average_num_soma_spikes_in_trial; // average number of somatic spikes in trials

    std::vector<int> num_trials_with_soma_spikes; // number of trials in which neuron produced somatic spikes
	std::vector<int> num_trials_with_dend_spikes; // number of trials in which neuron produced dendritic spikes


	if (MPI_rank == 0)
	{
		mean_first_soma_spike_time.resize(N_RA);
		std_first_soma_spike_time.resize(N_RA);
		mean_first_dend_spike_time.resize(N_RA);
		std_first_dend_spike_time.resize(N_RA);
		average_num_dend_spikes_in_trial.resize(N_RA);
		average_num_soma_spikes_in_trial.resize(N_RA);
		num_trials_with_soma_spikes.resize(N_RA);
		num_trials_with_dend_spikes.resize(N_RA);

		
		for (int i = 0; i < N_RA; i++)
		{
            average_num_dend_spikes_in_trial[i] = std::accumulate(num_dendritic_spikes_in_trials[i].begin(), num_dendritic_spikes_in_trials[i].end(), 0.0) 
                                                / static_cast<double>(num_trials);

            average_num_soma_spikes_in_trial[i] = std::accumulate(num_somatic_spikes_in_trials[i].begin(), num_somatic_spikes_in_trials[i].end(), 0.0) 
                                                / static_cast<double>(num_trials);
            
            num_trials_with_soma_spikes[i] = static_cast<int>(first_soma_spike_times[i].size());
			num_trials_with_dend_spikes[i] = static_cast<int>(first_dend_spike_times[i].size());
			
			///////////////////////////
			// First somatic spike jitter
			///////////////////////////////
            if (num_trials_with_soma_spikes[i] > 0)
            {

                //for (int j = 0; j < (int) average_dendritic_spike_time[i].size(); j++)
                  //  printf("average_dendritic_spike_time[%d][%d] = %f\n", i, j, average_dendritic_spike_time[i][j]);
				mean_first_soma_spike_time[i] = std::accumulate(first_soma_spike_times[i].begin(), first_soma_spike_times[i].end(), 0.0) / 
                                        static_cast<double>(num_trials_with_soma_spikes[i]);

            }

			else
				mean_first_soma_spike_time[i] = -1;
			// calculate standard deviation of burst times

			double accum = 0;
			double mean = mean_first_soma_spike_time[i];

			std::for_each(first_soma_spike_times[i].begin(), first_soma_spike_times[i].end(), [&accum, mean](const double t)
			{
				accum += (t - mean) * (t - mean);
			});

			if (static_cast<int>(first_soma_spike_times[i].size() > 1))
				std_first_soma_spike_time[i] = sqrt(accum / (static_cast<double>(first_soma_spike_times[i].size()) - 1));
			else
				std_first_soma_spike_time[i] = -1;
			
			///////////////////////////
			// First dendritic spike jitter
			///////////////////////////////

			if (num_trials_with_dend_spikes[i] > 0)
            {

                //for (int j = 0; j < (int) average_dendritic_spike_time[i].size(); j++)
                  //  printf("average_dendritic_spike_time[%d][%d] = %f\n", i, j, average_dendritic_spike_time[i][j]);
				mean_first_dend_spike_time[i] = std::accumulate(first_dend_spike_times[i].begin(), first_dend_spike_times[i].end(), 0.0) / 
                                        static_cast<double>(num_trials_with_dend_spikes[i]);

            }

			else
				mean_first_dend_spike_time[i] = -1;
			// calculate standard deviation of burst times

			accum = 0;
			mean = mean_first_dend_spike_time[i];

			std::for_each(first_dend_spike_times[i].begin(), first_dend_spike_times[i].end(), [&accum, mean](const double t)
			{
				accum += (t - mean) * (t - mean);
			});

			if (static_cast<int>(first_dend_spike_times[i].size() > 1))
				std_first_dend_spike_time[i] = sqrt(accum / (static_cast<double>(first_dend_spike_times[i].size()) - 1));
			else
				std_first_dend_spike_time[i] = -1;
		}
	}

	std::string fileJitter = outputDirectory + "jitter.bin";
	
	//~ std::cout << "Num trials with somatic spikes of HVC-RA neurons\n";
	//~ for (int i = 0; i < N_RA; i++)
	//~ {
		//~ std::cout << "HVC-RA " << i << " "
				  //~ << num_trials_with_soma_spikes[i] << std::endl;
	//~ }
	
	this->write_jitter(num_trials, 
						num_trials_with_soma_spikes, average_num_soma_spikes_in_trial, 
						mean_first_soma_spike_time, std_first_soma_spike_time, 
						num_trials_with_dend_spikes, average_num_dend_spikes_in_trial, 
						mean_first_dend_spike_time, std_first_dend_spike_time, 
						fileJitter.c_str());
	
	
}

void NetworkGrowthSimulator::new_chain_growth_with_inhibition_tracking(const ConfigurationNetworkGrowth &cfg, std::string networkDirectory, std::string fileTraining,  std::string outputDirectory, 
															bool training, int save_freq_long, int num_trials, double time_resolution_conductance)
{
	if (MPI_rank == 0)
	{
		connection_params = cfg.get_connection_parameters();
		noise_params = cfg.get_noise_parameters();
		maturation_params = cfg.get_maturation_parameters(); 
		synaptic_params = cfg.get_synaptic_parameters();
	}
	
	this->send_growth_parameters();
	
	if (MPI_rank == 0)
	{
		// read initial network from directory
		this->read_network_topology(networkDirectory);
		
		// write topology parameters to output directory
		ConfigurationNetworkTopology cfg_topology;
		cfg_topology.set_topology_parameters(topology_params);
		cfg_topology.write_configuration((outputDirectory + "network_parameters.cfg").c_str());
		
		this->sample_axonal_delays();
		this->resample_weights();
		
		this->read_training_neurons(fileTraining.c_str());
		
		// sort training neurons for convenience
		std::sort(training_neurons.begin(), training_neurons.end());
		
		// make training neurons mature
		for (size_t i = 0; i < training_neurons.size(); i++)
		{
			mature_global[training_neurons[i]] = 1;
			GCa_global[training_neurons[i]] = maturation_params.GCA_MATURE;
			rest_potential_global[training_neurons[i]] = maturation_params.E_REST_MATURE;
		}
		
		// initialize global inhibitory conductances
		this->initialize_global_inhibitory_conductances(time_resolution_conductance);
		
		// write network to output directory both as initial
		std::string filename_I_coordinates = outputDirectory + "I_xy.bin";
		std::string filename_out_training = outputDirectory + "training_neurons.bin";
		std::string filename_training_spread = outputDirectory + "training_spread.bin";
		std::string filename_num_neurons = outputDirectory + "num_neurons.bin";
		
		// generate fixed spread times
		double spread = 0.0;
		
		training_spread_times.resize(N_TR);
		for (int i = 0; i < N_TR; i++)
			training_spread_times[i] = noise_generator.random(spread) - spread / 2.;
	
	
		this->write_coordinates(xx_I,  yy_I,  zz_I,  filename_I_coordinates.c_str());
		this->write_training_neurons(filename_out_training.c_str());
		this->write_training_spread(filename_training_spread.c_str());
		this->write_number_of_neurons(filename_num_neurons.c_str());
		this->write_alterable_network("_initial", outputDirectory);
		
		
		cfg.write_configuration((outputDirectory + "growth_parameters.cfg").c_str());
	}
	
	this->initialize_network();
	this->initialize_local_inhibitory_conductances(time_resolution_conductance);
	
	//this->set_all_neurons_immature();
	//this->set_training_neurons_mature();
	
	// rescale inhibitory synapses to training neurons
	//for (size_t i = 0; i < training_neurons.size(); i++)
	//	this->rescale_synapses_to_mature(training_neurons[i]);
		
	//this->disable_RA2I_immature_outputs();
		
	this->chain_growth_with_inhibition_tracking(training, save_freq_long, num_trials, time_resolution_conductance, outputDirectory);
}

void NetworkGrowthSimulator::new_chain_growth(const ConfigurationNetworkGrowth &cfg, std::string networkDirectory, std::string fileTraining,  std::string outputDirectory, 
															bool training, int save_freq_short, int save_freq_long)
{
	if (MPI_rank == 0)
	{
		connection_params = cfg.get_connection_parameters();
		noise_params = cfg.get_noise_parameters();
		maturation_params = cfg.get_maturation_parameters(); 
		synaptic_params = cfg.get_synaptic_parameters();
	}
	
	this->send_growth_parameters();
	
	if (MPI_rank == 0)
	{
		// read initial network from directory
		this->read_network_topology(networkDirectory);
		
		// write topology parameters to output directory
		ConfigurationNetworkTopology cfg_topology;
		cfg_topology.set_topology_parameters(topology_params);
		cfg_topology.write_configuration((outputDirectory + "network_parameters.cfg").c_str());
		
		this->sample_axonal_delays();
		this->resample_weights();
		
		this->read_training_neurons(fileTraining.c_str());
		
		// sort training neurons for convenience
		std::sort(training_neurons.begin(), training_neurons.end());
		
		// make training neurons mature
		for (size_t i = 0; i < training_neurons.size(); i++)
		{
			mature_global[training_neurons[i]] = 1;
			GCa_global[training_neurons[i]] = maturation_params.GCA_MATURE;
			rest_potential_global[training_neurons[i]] = maturation_params.E_REST_MATURE;
	
		}
		
		// write network to output directory both as initial
		std::string filename_I_coordinates = outputDirectory + "I_xy.bin";
		std::string filename_out_training = outputDirectory + "training_neurons.bin";
		std::string filename_training_spread = outputDirectory + "training_spread.bin";
		std::string filename_num_neurons = outputDirectory + "num_neurons.bin";
		
		// generate fixed spread times
		double spread = 0.0;
		
		training_spread_times.resize(N_TR);
		for (int i = 0; i < N_TR; i++)
			training_spread_times[i] = noise_generator.random(spread) - spread / 2.;
	
	
		this->write_coordinates(xx_I,  yy_I,  zz_I,  filename_I_coordinates.c_str());
		this->write_training_neurons(filename_out_training.c_str());
		this->write_training_spread(filename_training_spread.c_str());
		this->write_number_of_neurons(filename_num_neurons.c_str());
		this->write_alterable_network("_initial", outputDirectory);
		
		
		cfg.write_configuration((outputDirectory + "growth_parameters.cfg").c_str());
	}
	
	this->initialize_network();
	
	//this->set_all_neurons_immature();
	//this->set_training_neurons_mature();
	
	// rescale inhibitory synapses to training neurons
	//for (size_t i = 0; i < training_neurons.size(); i++)
	//	this->rescale_synapses_to_mature(training_neurons[i]);
		
	//this->disable_RA2I_immature_outputs();
		
	this->chain_growth(training, save_freq_short, save_freq_long, outputDirectory);
}

void NetworkGrowthSimulator::read_network_state(std::string dataDir, int starting_trial)
{
	std::string fileEnding; // extension for file name, which takes into account trial number
	
	if (starting_trial > 0)
		fileEnding = "_" + std::to_string(starting_trial);
	else
		fileEnding = "";
	
	std::string fileRAxy = dataDir + "RA_xy" + fileEnding + ".bin";
	std::string fileIxy = dataDir + "I_xy.bin";
	
	std::string fileRA2I = dataDir + "RA_I_connections" + fileEnding + ".bin";
	std::string fileI2RA = dataDir + "I_RA_connections" + fileEnding + ".bin";
    
    std::string fileActive = dataDir + "RA_RA_active_connections" + fileEnding + ".bin";
    std::string fileSuper = dataDir + "RA_RA_super_connections" + fileEnding + ".bin";
    std::string fileMaturationProperties = dataDir + "maturation_properties" + fileEnding + ".bin";
	std::string fileActivityHistory = dataDir + "activity_history" + fileEnding + ".bin";
	//std::string fileLastBurstTimes = dataDir + "last_dendritic_spike_times" + trial_extension + ".bin";
	std::string fileReplacementHistory = dataDir + "replacement_history" + fileEnding + ".bin";
    std::string fileWeights = dataDir + "weights"  + fileEnding + ".bin";
    
    std::string fileAxonalDelaysRA2RA = dataDir + "axonal_delays_RA2RA" + fileEnding + ".bin";
    std::string fileAxonalDelaysRA2I = dataDir + "axonal_delays_RA2I" + fileEnding + ".bin";
    std::string fileAxonalDelaysI2RA = dataDir + "axonal_delays_I2RA" + fileEnding + ".bin";
    
    std::string fileRemodeledIndicators = dataDir + "remodeled_indicators" + fileEnding + ".bin";
    
    std::string fileTrainingSpread = dataDir + "training_spread.bin";
    std::string fileTraining = dataDir + "training_neurons.bin";
    std::string fileNumNeurons = dataDir + "num_neurons.bin";
   
   
	std::string fileTopologyConfiguration = dataDir + "network_parameters.cfg";
	std::string fileGrowthConfiguration = dataDir + "growth_parameters.cfg";
	
	// read configuration files
	
	ConfigurationNetworkGrowth cfg_growth;
	ConfigurationNetworkTopology cfg_topology;
	
	cfg_topology.read_configuration(fileTopologyConfiguration.c_str());
	cfg_growth.read_configuration(fileGrowthConfiguration.c_str());
   
	topology_params = cfg_topology.get_topology_parameters();
	
	connection_params = cfg_growth.get_connection_parameters();
	noise_params = cfg_growth.get_noise_parameters();
	maturation_params = cfg_growth.get_maturation_parameters(); 
	synaptic_params = cfg_growth.get_synaptic_parameters();
   
	this->read_number_of_neurons(fileNumNeurons.c_str());
	this->read_training_neurons(fileTraining.c_str());
	this->read_training_spread(fileTrainingSpread.c_str());
	
		
	// sort training neurons for convenience
	std::sort(training_neurons.begin(), training_neurons.end());
	
	this->resize_arrays_for_master_process();
	
	this->read_all_coordinates(fileRAxy.c_str(), fileIxy.c_str());
	this->read_connections_RAandI(fileRA2I.c_str(), fileI2RA.c_str());
	
	this->read_super_synapses(fileSuper.c_str());
    this->read_active_synapses(fileActive.c_str());
    this->read_weights(fileWeights.c_str());
    
    this->read_axonal_delays(axonal_delays_RA_RA_global, fileAxonalDelaysRA2RA.c_str());
    this->read_axonal_delays(axonal_delays_RA_I_global, fileAxonalDelaysRA2I.c_str());
    this->read_axonal_delays(axonal_delays_I_RA_global, fileAxonalDelaysI2RA.c_str());
    
    this->read_activity_history(fileActivityHistory.c_str());
    this->read_replacement_history(fileReplacementHistory.c_str());
    this->read_remodeled_indicators(fileRemodeledIndicators.c_str());
    this->read_maturation_properties(fileMaturationProperties.c_str());
}

//~ void NetworkGrowthSimulator::test_chain_recovery(std::string dataDir, int starting_trial, double fraction, bool training, int save_freq_short, int save_freq_long)
//~ {
	//~ this->read_network_state(dataDir, starting_trial);
    //~ 
    //~ 
    //~ this->set_time_for_neurons(trial_number * trial_duration);
	//~ this->make_lesion(fraction);
	//~ 
	//~ this->gather_data();
    //~ trial_number++;
	//~ //this->chain_growth(training, save_freq_short, save_freq_long);
//~ }

void NetworkGrowthSimulator::continue_chain_growth(std::string dataDir, std::string outDir, int starting_trial, bool training, int save_freq_short, int save_freq_long)
{
	trial_number = starting_trial;
	
	if (MPI_rank == 0)
		this->read_network_state(dataDir, starting_trial);
    
    this->send_growth_parameters();
    this->initialize_network();
    
    //this->set_time_for_neurons(trial_number * trial_duration);
    //this->print_super();
   
    
    
    this->gather_graph_state_data();
    this->gather_full_state_data();
    
    //if (MPI_rank == 0)
		//this->write_full_network_state("_" + std::to_string(trial_number) + "afterReading", dataDir);
	
	this->chain_growth(training, save_freq_short, save_freq_long, outDir);
	
}


void NetworkGrowthSimulator::set_noise()
{
	for (int i = 0; i < N_RA_local; i++)
	{	
		HVCRA_local[i].set_noise_generator(&noise_generator);
		HVCRA_local[i].set_white_noise(noise_params.white_noise_mean_soma, noise_params.white_noise_std_soma, 
									   noise_params.white_noise_mean_dend, noise_params.white_noise_std_dend);
	}

    for (int i = 0; i < N_I_local; i++)
	{
        HVCI_local[i].set_noise_generator(&noise_generator);
        HVCI_local[i].set_poisson_noise();
    }
}

//~ void NetworkGrowthSimulator::set_time_for_neurons(double t)
//~ {
	//~ for (int i = 0; i < N_RA_local; i++)
		//~ HVCRA_local[i].set_time(t);
//~ 
	//~ 
    //~ for (int i = 0; i < N_I_local; i++)
		//~ HVCI_local[i].set_time(t);
//~ }

void NetworkGrowthSimulator::set_dynamics()
{
	for (int i = 0; i < N_RA_local; i++)
		HVCRA_local[i].set_dynamics(TIMESTEP);

	
    for (int i = 0; i < N_I_local; i++)
		HVCI_local[i].set_dynamics(TIMESTEP);
}

void NetworkGrowthSimulator::set_neuron_properties_sudden_maturation()
{
	for (int i = 0; i < N_RA_local; i++)
	{
		//if ( mature_global[Id_RA_local[i]] == 1 )
		if ( mature_local[i] == 1 )
			this->set_neuron_mature(i);
		else
			this->set_neuron_immature(i);
	}
}

void NetworkGrowthSimulator::update_replaced_neurons(std::vector<int>& neurons_to_replace)
{
	for (size_t i = 0; i < neurons_to_replace.size(); i++)
	{
		int neuron_id = neurons_to_replace[i];
		int rank, shift;
		
		this->get_neuronRA_location(neuron_id, &rank, &shift);
		
		if (MPI_rank == rank)
		{
			rest_potential_local[shift] = maturation_params.E_REST_IMMATURE;
			GCa_local[shift] = maturation_params.GCA_IMMATURE;
			
			mature_local[shift] = 0;
			maturation_scale_local[shift] = MATURATION_SCALE_SPONTANEOUS;
			remodeled_local[shift] = 0;
			
			num_trials_after_replacement_local[shift] = 0;
			
			for (int j = 0; j < maturation_params.RATE_WINDOW_LONG; j++)
				num_spikes_in_recent_trials_local[shift].push_back(0);
	
		}
	}
}

void NetworkGrowthSimulator::set_neuron_properties()
{
	for (int i = 0; i < N_RA_local; i++)
	{
		//HVCRA_local[i].set_Ei(gaba_potential_local[i]);
		HVCRA_local[i].set_Erest(rest_potential_local[i]);
		HVCRA_local[i].set_GCa(GCa_local[i]);
		
		//HVCRA_local[i].set_Gk(Gk_local[i]);
		//HVCRA_local[i].set_GNa(GNa_local[i]);
		//HVCRA_local[i].set_Ad(Ad_local[i]);
		//HVCRA_local[i].set_Rc(Rc_local[i]);
	}
}


void NetworkGrowthSimulator::update_neuron_properties_sameDendrite_diffMaturationRate()
{
	for (int i = 0; i < N_RA_local; i++)
	{
		// update only properties of the neurons that are not mature
		if ( mature_local[i] != 1){
			
			//rest_potential_local[i] = rest_potential_local[i] * exp(1.0 / maturation_scale_local[i]);
			//GCa_local[i] = GCa_local[i] * exp(-1.0 / maturation_scale_local[i]);
			
			rest_potential_local[i] = rest_potential_local[i]  + (rest_potential_local[i] - maturation_params.E_REST_MATURE) 
														* ( exp(-1.0 / maturation_scale_local[i]) - 1);
			
			GCa_local[i] = GCa_local[i]  + (GCa_local[i] - maturation_params.GCA_MATURE) 
														* ( exp(-1.0 / maturation_scale_local[i]) - 1);
			
			double potential_threshold = 0.5;
			double conductance_threshold = 0.5;
			
			if ( ( GCa_local[i] >= maturation_params.GCA_MATURE - conductance_threshold) || 
					(rest_potential_local[i] <= maturation_params.E_REST_MATURE - potential_threshold) ){
				GCa_local[i] = maturation_params.GCA_MATURE;
				rest_potential_local[i] = maturation_params.E_REST_MATURE;
				
				mature_local[i] = 1;
			}
			
			HVCRA_local[i].set_Erest(rest_potential_local[i]);
			HVCRA_local[i].set_GCa(GCa_local[i]);
		}
	}
}


void NetworkGrowthSimulator::update_neuron_properties_sameDendrite()
{
	for (int i = 0; i < N_RA_local; i++)
	{
		if ( mature_global[Id_RA_local[i]] != 1){
			// update only properties of the neurons that are not mature
			double age = static_cast<int>(num_trials_after_replacement_local[i]);
			

			double rest_potential = maturation_params.E_REST_MATURE + 
										(maturation_params.E_REST_IMMATURE - maturation_params.E_REST_MATURE) * 
																	exp(-age / MATURATION_SCALE_SPONTANEOUS);
			
			double GCa = maturation_params.GCA_MATURE + (maturation_params.GCA_IMMATURE - maturation_params.GCA_MATURE) * exp(-age / MATURATION_SCALE_SPONTANEOUS);
			
			
			HVCRA_local[i].set_Erest(rest_potential);
			HVCRA_local[i].set_GCa(GCa);
		}
	}
}

void NetworkGrowthSimulator::update_neuron_properties()
{
	for (int i = 0; i < N_RA_local; i++)
	{
		// update only properties of the neurons that are not mature
		if ( mature_local[i] == 0 )
		{
		
			double age = static_cast<int>(num_trials_after_replacement_local[i]);
			
			gaba_potential_local[i] = maturation_params.E_GABA_MATURE + 
										(maturation_params.E_GABA_IMMATURE - maturation_params.E_GABA_MATURE) * 
																	exp(-age / MATURATION_SCALE_SPONTANEOUS);
			
			rest_potential_local[i] = maturation_params.E_REST_MATURE + 
										(maturation_params.E_REST_IMMATURE - maturation_params.E_REST_MATURE) * 
																	exp(-age / MATURATION_SCALE_SPONTANEOUS);
			
			Gk_local[i] = maturation_params.GK_MATURE + 
						(maturation_params.GK_IMMATURE - maturation_params.GK_MATURE) * 
																	exp(-age / MATURATION_SCALE_SPONTANEOUS);														
			
			GNa_local[i] = maturation_params.GNA_MATURE + 
						(maturation_params.GNA_IMMATURE - maturation_params.GNA_MATURE) * 
																	exp(-age / MATURATION_SCALE_SPONTANEOUS);
																	
			Ad_local[i] = maturation_params.AD_MATURE + 
						(maturation_params.AD_IMMATURE - maturation_params.AD_MATURE) * 
																	exp(-age / MATURATION_SCALE_SPONTANEOUS);	
																														
			Rc_local[i] = maturation_params.RC_MATURE + 
						(maturation_params.RC_IMMATURE - maturation_params.RC_MATURE) * 
																		exp(-age / MATURATION_SCALE_SPONTANEOUS);
		}
	}
}

void NetworkGrowthSimulator::set_synapse_indicators()
{
	for (int i = 0; i < N_RA_local; i++)
	{
		// set active synapse indicators
		for (size_t j = 0; j < active_synapses_local[i].size(); j++)
		{
			int target_id = active_synapses_local[i][j];
			active_indicators_local[i][target_id] = 1;
		}
		
		// set super synapses
		for (size_t j = 0; j < supersynapses_local[i].size(); j++)
		{
			int target_id = supersynapses_local[i][j];
			supersynapses_indicators_local[i][target_id] = 1;
		}	
	}
}

void NetworkGrowthSimulator::set_neuron_immature(int local_id)
{
	//~ gaba_potential_local[local_id] = maturation_params.E_GABA_IMMATURE;
	//~ rest_potential_local[local_id] = maturation_params.E_REST_IMMATURE;
	//~ Gk_local[local_id] = maturation_params.GK_IMMATURE;
	//~ GNa_local[local_id] = maturation_params.GNA_IMMATURE;
	//~ Ad_local[local_id] = maturation_params.AD_IMMATURE;
	//~ Rc_local[local_id] = maturation_params.RC_IMMATURE;
	//~ 
	//~ HVCRA_local[local_id].set_Ei(gaba_potential_local[local_id]);
	//~ HVCRA_local[local_id].set_Erest(rest_potential_local[local_id]);
	//~ HVCRA_local[local_id].set_Gk(Gk_local[local_id]);
	//~ HVCRA_local[local_id].set_GNa(GNa_local[local_id]);
	//~ HVCRA_local[local_id].set_Ad(Ad_local[local_id]);
	//~ HVCRA_local[local_id].set_Rc(Rc_local[local_id]);
	rest_potential_local[local_id] = maturation_params.E_REST_IMMATURE;
	GCa_local[local_id] = maturation_params.GCA_IMMATURE;
	
	
	HVCRA_local[local_id].set_Ei(maturation_params.E_GABA_IMMATURE);
	HVCRA_local[local_id].set_Erest(maturation_params.E_REST_IMMATURE);
	HVCRA_local[local_id].set_Gk(maturation_params.GK_IMMATURE);
	HVCRA_local[local_id].set_GNa(maturation_params.GNA_IMMATURE);
	HVCRA_local[local_id].set_GCa(maturation_params.GCA_IMMATURE);
	HVCRA_local[local_id].set_GCaK(maturation_params.GCAK_IMMATURE);
	HVCRA_local[local_id].set_GsL(maturation_params.GSL_IMMATURE);
	HVCRA_local[local_id].set_GdL(maturation_params.GDL_IMMATURE);
	HVCRA_local[local_id].set_Ad(maturation_params.AD_IMMATURE);
	HVCRA_local[local_id].set_Rc(maturation_params.RC_IMMATURE);
}

void NetworkGrowthSimulator::set_neuron_mature(int local_id)
{
	//~ gaba_potential_local[local_id] = maturation_params.E_GABA_MATURE;
	//~ rest_potential_local[local_id] = maturation_params.E_REST_MATURE;
	//~ Gk_local[local_id] = maturation_params.GK_MATURE;
	//~ GNa_local[local_id] = maturation_params.GNA_MATURE;
	//~ Ad_local[local_id] = maturation_params.AD_MATURE;
	//~ Rc_local[local_id] = maturation_params.RC_MATURE;
	
	rest_potential_local[local_id] = maturation_params.E_REST_MATURE;
	GCa_local[local_id] = maturation_params.GCA_MATURE;
	
	//~ HVCRA_local[local_id].set_Ei(gaba_potential_local[local_id]);
	//~ HVCRA_local[local_id].set_Erest(rest_potential_local[local_id]);
	//~ HVCRA_local[local_id].set_Gk(Gk_local[local_id]);
	//~ HVCRA_local[local_id].set_GNa(GNa_local[local_id]);
	//~ HVCRA_local[local_id].set_Ad(Ad_local[local_id]);
	//~ HVCRA_local[local_id].set_Rc(Rc_local[local_id]);
	//~ 
	HVCRA_local[local_id].set_Ei(maturation_params.E_GABA_MATURE);
	HVCRA_local[local_id].set_Erest(maturation_params.E_REST_MATURE);
	HVCRA_local[local_id].set_Gk(maturation_params.GK_MATURE);
	HVCRA_local[local_id].set_GNa(maturation_params.GNA_MATURE);
	HVCRA_local[local_id].set_GCa(maturation_params.GCA_MATURE);
	HVCRA_local[local_id].set_GCaK(maturation_params.GCAK_MATURE);
	HVCRA_local[local_id].set_GsL(maturation_params.GSL_MATURE);
	HVCRA_local[local_id].set_GdL(maturation_params.GDL_MATURE);
	HVCRA_local[local_id].set_Ad(maturation_params.AD_MATURE);
	HVCRA_local[local_id].set_Rc(maturation_params.RC_MATURE);
}

void NetworkGrowthSimulator::set_active_synapse_weights(double G)
{
	for (int i = 0; i < N_RA_local; i++)
		for (int j = 0; j < N_RA; j++ )
			if ( active_indicators_local[i][j] == 1 )
				weights_RA_RA_local[i][j] = G;
	
}

void NetworkGrowthSimulator::set_all_neurons_immature()
{
	// set all to immature
	for (int i = 0; i < N_RA_local; i++)
	{
		gaba_potential_local[i] = maturation_params.E_GABA_IMMATURE;
		rest_potential_local[i] = maturation_params.E_REST_IMMATURE;
		Gk_local[i] = maturation_params.GK_IMMATURE;
		GNa_local[i] = maturation_params.GNA_IMMATURE;
		Ad_local[i] = maturation_params.AD_IMMATURE;
		Rc_local[i] = maturation_params.RC_IMMATURE;
	}
	
	std::fill(mature_local.begin(), mature_local.end(), 0);
}

void NetworkGrowthSimulator::set_recording(const std::vector<int> &RA_neurons, const std::vector<int> &I_neurons, int trial_number,
																std::string outputDirectory)
{
	std::string RAdir = outputDirectory + "RA/";
	std::string Idir = outputDirectory + "I/";
	
	
	// check if files already exist: then delete them
	if (MPI_rank == 0)
	{
		struct stat buf;
	
		for (size_t i = 0; i < RA_neurons.size(); i++)
		{
			std::string filename = RAdir + "RA" + std::to_string(RA_neurons[i]) + ".bin";
			
			if ( stat(filename.c_str(), &buf) == 0 )
				std::remove(filename.c_str());
		}
		
		for (size_t i = 0; i < I_neurons.size(); i++)
		{
			std::string filename = Idir + "I" + std::to_string(I_neurons[i]) + ".bin";
			
			if ( stat(filename.c_str(), &buf) == 0 )
				std::remove(filename.c_str());
		}
	}
	
	int rank;
	int shift;
	
	for (size_t i = 0; i < RA_neurons.size(); i++)
	{
		this->get_neuronRA_location(RA_neurons[i], &rank, &shift);
		
		if (MPI_rank == rank)
		{	
			std::string filename = RAdir + "RA" + std::to_string(RA_neurons[i]) + "_trial" + std::to_string(trial_number) + ".bin";
			//HVCRA_local[shift].set_recording_full(filename.c_str());
			HVCRA_local[shift].set_recording(filename.c_str());
		}		
	}
	
	for (size_t i = 0; i < I_neurons.size(); i++)
	{
		this->get_neuronI_location(I_neurons[i], &rank, &shift);
		
		if (MPI_rank == rank)
		{	
			std::string filename = Idir + "I" + std::to_string(I_neurons[i]) + "_trial" + std::to_string(trial_number) + ".bin";
			HVCI_local[shift].set_recording(filename.c_str());
		}		
	}
	
}

void NetworkGrowthSimulator::set_training_neurons_mature()
{
	// set training neurons to be mature
	for (int i = 0; i < N_TR; i++)
	{
		mature_global[training_neurons[i]] = 1;
		
		int rank;
		int shift;
		
		this->get_neuronRA_location(training_neurons[i], &rank, &shift);
		
		if (MPI_rank == rank)
		{	
			gaba_potential_local[shift] = maturation_params.E_GABA_MATURE;
			rest_potential_local[shift] = maturation_params.E_REST_MATURE;
			Gk_local[shift] = maturation_params.GK_MATURE;
			GNa_local[shift] = maturation_params.GNA_MATURE;
			Ad_local[shift] = maturation_params.AD_MATURE;
			Rc_local[shift] = maturation_params.RC_MATURE;
			
			mature_local[shift] = 1;
		}		
	}	
}

void NetworkGrowthSimulator::initialize_generator()
{
    // prepare seeds and initialize each generator with different seed
    srand((unsigned int) time(NULL));

    //unsigned seed = rand();
    //unsigned seed = rand() + MPI_rank;

    unsigned seed = 15000 + MPI_rank*1000;

    noise_generator.set_seed(seed);
}

void NetworkGrowthSimulator::print_simulation_parameters()
{

	printf("\nGrowth parameters: \n");

	printf("\nNoise:\n");
	
	printf("\nwhite_noise_mean_soma = %f\n", noise_params.white_noise_mean_soma);
	printf("white_noise_std_soma = %f\n", noise_params.white_noise_std_soma);
	printf("white_noise_mean_dend = %f\n", noise_params.white_noise_mean_dend);
	printf("white_noise_std_dend = %f\n", noise_params.white_noise_std_dend);

	


	printf("\nSynaptic parameters\n");
	
	printf("\nNss = %d\n", synaptic_params.Nss);
	
	printf("\nA_P = %f\n", synaptic_params.A_P);
	printf("A_D = %f\n", synaptic_params.A_D);
	printf("T_0 = %f\n", synaptic_params.T_0);
	printf("T_P = %f\n", synaptic_params.T_P);
	printf("TAU_P = %f\n", synaptic_params.TAU_P);
	printf("T_D = %f\n", synaptic_params.T_D);
	printf("TAU_D = %f\n", synaptic_params.TAU_D);
	
	printf("\nBETA = %f\n", synaptic_params.BETA);
	printf("BETA_SUPERSYNAPSE = %f\n", synaptic_params.BETA_SUPERSYNAPSE);
	
	printf("\nACTIVATION_THRESHOLD = %f\n", synaptic_params.ACTIVATION_THRESHOLD);
	printf("SUPERSYNAPSE_THRESHOLD = %f\n", synaptic_params.SUPERSYNAPSE_THRESHOLD);
	printf("WEIGHT_MAX = %f\n", synaptic_params.WEIGHT_MAX);
	
	
	printf("\nMaturation parameters\n");
	
	printf("\nE_GABA_MATURE = %f\n", maturation_params.E_GABA_MATURE);
	printf("E_GABA_IMMATURE = %f\n", maturation_params.E_GABA_IMMATURE);
	
	printf("DEATH_RATE_THRESHOLD = %f\n", maturation_params.DEATH_RATE_THRESHOLD);
	
	printf("RATE_WINDOW_LONG = %d\n", maturation_params.RATE_WINDOW_LONG);

}

//~ void NetworkGrowthSimulator::read_maturation_info(const char* filename)
//~ {
    //~ std::ifstream inp;
//~ 
    //~ inp.open(filename, std::ios::in | std::ios::binary);
//~ 
    //~ // read number of HVC(RA) neurons
//~ 
  //~ int N;
//~ 
    //~ inp.read(reinterpret_cast<char *>(&N), sizeof(N));
//~ 
    //~ if (N != N_RA)
        //~ std::cerr << "Number of HVC(RA) neurons read from file with active synapses: N = " << N << "is different from N_RA = " << N_RA << std::endl;
//~ 
    //~ inp.read(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));
//~ 
    //~ if (MPI_rank == 0)
        //~ std::cout << "Trial number read from maturation file: trial_number = " << trial_number << std::endl;
    //~ 
    //~ int counter = 0; // counter of neuron id in the process
   //~ 
    //~ for (int i = 0; i < N_RA; i++)
    //~ {
        //~ bool neuron_found = false; // indicator that we found a neuron
//~ 
        //~ if (counter < N_RA_local)
            //~ if (Id_RA_local[counter] == i)
            //~ {
                //~ neuron_found = true;
            //~ }
//~ 
        //~ double gaba; // reverse GABA potential
        //~ double firing_rate_short; // neuronal firing rate in small window
        //~ double firing_rate_long; // neuronal firing rate in large window
        //~ int remodeled; // indicator of neuronal axon remodeling
      //~ 
        //~ inp.read(reinterpret_cast<char *>(&gaba), sizeof(gaba));
        //~ inp.read(reinterpret_cast<char *>(&firing_rate_short), sizeof(firing_rate_short));
        //~ inp.read(reinterpret_cast<char *>(&firing_rate_long), sizeof(firing_rate_long));
        //~ inp.read(reinterpret_cast<char *>(&remodeled), sizeof(remodeled));
       //~ 
//~ 
        //~ if (neuron_found)
        //~ {
            //~ gaba_potential_local[counter] = gaba;
            //~ firing_rate_short_local[counter] = firing_rate_short;
            //~ firing_rate_long_local[counter] = firing_rate_long;
            //~ remodeled_local[counter] = remodeled;
            //~ 
        //~ }
//~ 
        //~ // increase counter
        //~ if (neuron_found)
            //~ counter++;
//~ 
    //~ }
    //~ inp.close();
//~ }

void NetworkGrowthSimulator::read_training_spread(const char* filename)
{
	struct stat buf;
	
	if ( stat(filename, &buf) == 1 )
	{
		std::cerr << " File " << filename << " doesn't exist!\n" << std::endl;
		return;
	}
	
	std::ifstream in;
  
	in.open(filename, std::ios::in | std::ios::binary );
	
	int N;
	
	in.read(reinterpret_cast<char *>(&N), sizeof(int));

	if (N != N_TR)
		std::cerr << "Number of neurons read from file with training spread: N = " << N << "is different from N_TR = " << N_TR << std::endl;
	

	training_spread_times.resize(N_TR);

	for (int i = 0; i < N_TR; i++)
		in.read(reinterpret_cast<char *>(&training_spread_times[i]), sizeof(double));
		
	in.close();
}

void NetworkGrowthSimulator::read_training_neurons(const char* filename)
{
	struct stat buf;
	
	if ( stat(filename, &buf) == 1 )
	{
		std::cerr << " File " << filename << " doesn't exist!\n" << std::endl;
		return;
	}
	
	std::ifstream in;
  
	in.open(filename, std::ios::in | std::ios::binary );
	
	in.read(reinterpret_cast<char *>(&N_TR), sizeof(int));

	training_neurons.resize(N_TR);

	for (int i = 0; i < N_TR; i++)
		in.read(reinterpret_cast<char *>(&training_neurons[i]), sizeof(training_neurons[i]));
		
	in.close();
}

void NetworkGrowthSimulator::read_number_of_neurons(const char* filename)
{
	struct stat buf;
	
	if ( stat(filename, &buf) == 1 )
	{
		std::cerr << " File " << filename << " doesn't exist!\n" << std::endl;
		return;
	}
	
	std::ifstream in;
  
	in.open(filename, std::ios::in | std::ios::binary );
	
	in.read(reinterpret_cast<char *>(&N_RA), sizeof(int));
	in.read(reinterpret_cast<char *>(&N_I) , sizeof(int));

	in.close();
}

void NetworkGrowthSimulator::read_network_topology(std::string networkDir)
{
	std::string filename_I_coordinates = networkDir + "I_xy.bin";
	std::string filename_RA_coordinates = networkDir + "RA_xy.bin";
	
	std::string filename_RA_I = networkDir + "RA_I_connections.bin";
	std::string filename_I_RA = networkDir + "I_RA_connections.bin";
	
	std::string filename_num_neurons = networkDir + "num_neurons.bin";
	
	std::string filename_topology_configuration = networkDir + "network_parameters.cfg";
	
	ConfigurationNetworkTopology cfg;
	cfg.read_configuration(filename_topology_configuration.c_str());
	
	topology_params = cfg.get_topology_parameters();
	
	this->read_number_of_neurons(filename_num_neurons.c_str());
	this->resize_arrays_for_master_process();
	
	this->read_all_coordinates(filename_RA_coordinates.c_str(), filename_I_coordinates.c_str());
	this->read_connections_RAandI(filename_RA_I.c_str(), filename_I_RA.c_str());
}

void NetworkGrowthSimulator::read_all_coordinates(const char *RA_xy, const char *I_xy)
{
	struct stat buf;
			
	if ( stat(RA_xy, &buf) == 1 )
	{
		std::cerr << " File " << RA_xy << " doesn't exist!\n" << std::endl;
		return;
	}
	
	if ( stat(I_xy, &buf) == 1 )
	{
		std::cerr << " File " << I_xy << " doesn't exist!\n" << std::endl;
		return;
	}
	
	std::ifstream inp_RA, inp_I;

	// open files
	inp_RA.open(RA_xy, std::ios::binary | std::ios::in);
	inp_I.open(I_xy, std::ios::binary | std::ios::in);

	// read number of neurons
	
	int N_ra, N_i;
	
	inp_RA.read(reinterpret_cast<char*>(&N_ra), sizeof(int));
	inp_I.read(reinterpret_cast<char*>(&N_i), sizeof(int));
	
	if (N_ra != N_RA)
		std::cerr << "Number of HVC(RA) neurons read from file with fixed synapses: N_ra = " << N_ra << "is different from N_RA = " << N_RA << std::endl;
	
	if (N_i != N_I)
		std::cerr << "Number of HVC(I) neurons read from file with fixed synapses: N_i = " << N_i << "is different from N_I = " << N_I << std::endl;
	
	int dim_RA, dim_I; 
	
	inp_RA.read(reinterpret_cast<char*>(&dim_RA), sizeof(int));
	inp_I.read(reinterpret_cast<char*>(&dim_I), sizeof(int));
	
	if (dim_RA != dim_I)
		std::cerr << "Dimensionality of HVC-RA coord = " << dim_RA << " is different from dimensionality for HVC-I coord = " << dim_I << std::endl;
	
	double model_distance;
	
	inp_RA.read(reinterpret_cast<char*>(&MODEL_INTERNEURON_DISTANCE), sizeof(double));
	inp_I.read(reinterpret_cast<char*>(&model_distance), sizeof(double));
	
	if (MODEL_INTERNEURON_DISTANCE != model_distance)
		std::cerr << "Model interneuron distance of HVC-RA coord = " << MODEL_INTERNEURON_DISTANCE 
					<< " is different from model interneuron distance for HVC-I coord = " << model_distance << std::endl;
	
	// resize arrays
	xx_RA.resize(N_RA);
	yy_RA.resize(N_RA);
	
	xx_I.resize(N_I);
	yy_I.resize(N_I);
	
	if ( dim_RA == 3 )
	{
		zz_RA.resize(N_RA);
		zz_I.resize(N_I);
	}
	
	for (int i = 0; i < N_RA; i++)
	{
		inp_RA.read(reinterpret_cast<char*>(&xx_RA[i]), sizeof(xx_RA[i]));
		inp_RA.read(reinterpret_cast<char*>(&yy_RA[i]), sizeof(yy_RA[i]));
		
		if (dim_RA == 3)
			inp_RA.read(reinterpret_cast<char*>(&zz_RA[i]), sizeof(zz_RA[i]));
	}

	// read I coordinates
	for (int i = 0; i < N_I; i++)
	{
		inp_I.read(reinterpret_cast<char*>(&xx_I[i]), sizeof(xx_I[i]));
		inp_I.read(reinterpret_cast<char*>(&yy_I[i]), sizeof(yy_I[i]));
		
		if (dim_I == 3)
			inp_I.read(reinterpret_cast<char*>(&zz_I[i]), sizeof(zz_I[i]));
	}
	
	// close files
	inp_RA.close();	
	inp_I.close();	
	
}

void NetworkGrowthSimulator::read_connections_RAandI(const char* RA_I, const char* I_RA)
{
	struct stat buf;
	
	if ( stat(RA_I, &buf) == 1 )
	{
		std::cerr << " File " << RA_I << " doesn't exist!\n" << std::endl;
		return;
	}
	
	if ( stat(I_RA, &buf) == 1 )
	{
		std::cerr << " File " << I_RA << " doesn't exist!\n" << std::endl;
		return;
	}
	
	std::ifstream inp_RA_I, inp_I_RA;
        
	// open files
	inp_RA_I.open(RA_I, std::ios::binary | std::ios::in);
	inp_I_RA.open(I_RA, std::ios::binary | std::ios::in);

	int N_ra;
	int N_i;

	inp_RA_I.read(reinterpret_cast<char *>(&N_ra), sizeof(int));
	inp_I_RA.read(reinterpret_cast<char *>(&N_i), sizeof(int));

	if (N_ra != N_RA)
		std::cerr << "Number of HVC(RA) neurons read from file with fixed synapses: N_ra = " << N_ra << " is different from N_RA = " << N_RA << std::endl;
	
	if (N_i != N_I)
		std::cerr << "Number of HVC(I) neurons read from file with fixed synapses: N_i = " << N_i << " is different from N_I = " << N_I << std::endl;
	
	// read connections from RA to I
	for (int i = 0; i < N_RA; i++)
	{
		int n_id; // neuronal id
		int size; // number of outgoing connections

		inp_RA_I.read(reinterpret_cast<char *>(&n_id), sizeof(int));
		inp_RA_I.read(reinterpret_cast<char *>(&size), sizeof(int)); // write neuron's ID

		syn_ID_RA_I_global[i].resize(size);
		weights_RA_I_global[i].resize(size);
		syn_lengths_RA_I_global[i].resize(size);
		axonal_delays_RA_I_global[i].resize(size);
		
		
		for (int j = 0; j < size; j++)
		{
			inp_RA_I.read(reinterpret_cast<char *>(&syn_ID_RA_I_global[i][j]), sizeof(int));
			inp_RA_I.read(reinterpret_cast<char *>(&weights_RA_I_global[i][j]), sizeof(double));
			inp_RA_I.read(reinterpret_cast<char *>(&syn_lengths_RA_I_global[i][j]), sizeof(double));
			inp_RA_I.read(reinterpret_cast<char *>(&axonal_delays_RA_I_global[i][j]), sizeof(double));
		}

	}

	// read connections from I to RA
	for (int i = 0; i < N_I; i++)
	{
		int n_id; // neuronal id
		int size; // number of outgoing connections

		inp_I_RA.read(reinterpret_cast<char *>(&n_id), sizeof(n_id));
		inp_I_RA.read(reinterpret_cast<char *>(&size), sizeof(size)); // write neuron's ID
		
		syn_ID_I_RA_global[i].resize(size);
		weights_I_RA_global[i].resize(size);
		syn_lengths_I_RA_global[i].resize(size);
		axonal_delays_I_RA_global[i].resize(size);
		
		for (int j = 0; j < size; j++)
		{
			inp_I_RA.read(reinterpret_cast<char *>(&syn_ID_I_RA_global[i][j]), sizeof(int)); // write targets ID
			inp_I_RA.read(reinterpret_cast<char *>(&weights_I_RA_global[i][j]), sizeof(double)); // write targets conductance
			inp_I_RA.read(reinterpret_cast<char *>(&syn_lengths_I_RA_global[i][j]), sizeof(double));
			inp_I_RA.read(reinterpret_cast<char *>(&axonal_delays_I_RA_global[i][j]), sizeof(double));
		}
	}
	// close files
	inp_I_RA.close();
	inp_RA_I.close();
}

//~ void NetworkGrowthSimulator::read_replacement_history(const char* filename)
//~ {
	//~ std::ifstream inp;
//~ 
	//~ inp.open(filename, std::ios::in | std::ios::binary);
//~ 
	//~ // read number of HVC(RA) neurons
//~ 
	//~ int N;
//~ 
	//~ inp.read(reinterpret_cast<char *>(&N), sizeof(N));
//~ 
	//~ if (N != N_RA)
		//~ std::cerr << "Number of HVC(RA) neurons read from file with replacement history: N = " << N << "is different from N_RA = " << N_RA << std::endl;
//~ 
	//~ inp.read(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));
//~ 
	//~ if (MPI_rank == 0)
		//~ std::cout << "Trial number read from replacement history file: trial_number = " << trial_number << std::endl;
	//~ 
	//~ int counter = 0; // counter of neuron id in the process
   //~ 
	//~ for (int i = 0; i < N_RA; i++)
	//~ {
		//~ bool neuron_found = false; // indicator that we found a neuron
//~ 
		//~ if (counter < N_RA_local)
			//~ if (Id_RA_local[counter] == i)
			//~ {
				//~ neuron_found = true;
			//~ }
//~ 
	//~ 
		//~ int trials_since_replacement; // number of trials since previous neuron replacement
	  //~ 
		//~ inp.read(reinterpret_cast<char *>(&trials_since_replacement), sizeof(trials_since_replacement));
	  //~ 
//~ 
		//~ if (neuron_found)
		//~ {
			//~ num_trials_after_replacement_local[counter] = trials_since_replacement;
			//~ counter++;
		//~ }
//~ 
	//~ }
	//~ 
	//~ 
	//~ inp.close();
//~ }
//~ 
//~ void NetworkGrowthSimulator::read_active_synapses(const char *filename)
//~ {
    //~ std::ifstream inp;
//~ 
    //~ inp.open(filename, std::ios::in | std::ios::binary);
//~ 
    //~ // read number of HVC(RA) neurons
//~ 
  //~ int N;
//~ 
    //~ inp.read(reinterpret_cast<char *>(&N), sizeof(N));
//~ 
    //~ if (N != N_RA)
        //~ std::cerr << "Number of HVC(RA) neurons read from file with active synapses: N = " << N << "is different from N_RA = " << N_RA << std::endl;
//~ 
    //~ int counter = 0; // counter of neuron id in the process
//~ 
    //~ for (int i = 0; i < N_RA; i++)
    //~ {
        //~ // read neuron's ID number
        //~ int n_id;
        //~ 
        //~ inp.read(reinterpret_cast<char *>(&n_id), sizeof(n_id));
        //~ 
        //~ // read number of connections to RA neurons
//~ 
      //~ int num_active;
        //~ bool source_neuron_found = false; // indicator that we found a source neuron
//~ 
        //~ inp.read(reinterpret_cast<char *>(&num_active), sizeof(num_active));
        //~ 
        //~ if (counter < N_RA_local)
            //~ if (Id_RA_local[counter] == i)
            //~ {
                //~ active_synapses_local[counter].resize(num_active);
                //~ source_neuron_found = true;
            //~ }
//~ 
        //~ // read active synapses
        //~ for (int j = 0; j < num_active; j++)
        //~ {
            //~ int target_id; // id of active target
            //~ double synapse_weight; // weight of synapse
//~ 
            //~ inp.read(reinterpret_cast<char *>(&target_id), sizeof(target_id));
            //~ inp.read(reinterpret_cast<char *>(&synapse_weight), sizeof(synapse_weight));
//~ 
            //~ if (source_neuron_found)
            //~ {
                //~ active_synapses_local[counter][j] = target_id;
                //~ active_indicators_local[counter][target_id] = 1;
            //~ }
//~ 
        //~ }
        //~ // increase counter
        //~ if (source_neuron_found)
            //~ counter++;
//~ 
    //~ }
    //~ inp.close();
//~ }

//~ void NetworkGrowthSimulator::print_active()
//~ {
    //~ for (int i = 0; i < N_RA_local; i++)
        //~ for (int j = 0; j < N_RA; j++)
            //~ if (active_indicators_local[i][j] == 1)
                //~ std::cout << "Active synapse " << Id_RA_local[i] << " -> " << j << std::endl; 
//~ }
//~ 
//~ void NetworkGrowthSimulator::print_super()
//~ {
    //~ for (int i = 0; i < N_RA_local; i++)
        //~ for (int j = 0; j < N_RA; j++)
            //~ if (supersynapses_indicators_local[i][j] == 1)
                //~ std::cout << "Active supersynapse " << Id_RA_local[i] << " -> " << j << std::endl; 
//~ }

//~ void NetworkGrowthSimulator::read_last_dendritic_spike_times(const char* filename)
//~ {
	//~ std::ifstream inp;
//~ 
    //~ inp.open(filename, std::ios::in | std::ios::binary);
//~ 
    //~ // read number of HVC(RA) neurons
//~ 
	//~ int N;
//~ 
    //~ inp.read(reinterpret_cast<char *>(&N), sizeof(N));
//~ 
    //~ if (N != N_RA)
        //~ std::cerr << "Number of HVC(RA) neurons read from file with last dendritic spike times: N = " << N 
				//~ << "is different from N_RA = " << N_RA << std::endl;
//~ 
	//~ 
//~ 
    //~ int num_trial;
    //~ 
    //~ inp.read(reinterpret_cast<char *>(&num_trial), sizeof(num_trial));
    //~ 
    //~ if (MPI_rank == 0)
        //~ std::cout << "Trial number read from file with last dendritic spike times: " << num_trial << std::endl;
//~ 
    //~ int counter = 0; // counter of neuron id in the process
	//~ double spike_time; // last dendritic spike time
//~ 
    //~ for (int i = 0; i < N_RA; i++)
    //~ {
        //~ bool source_neuron_found = false; // indicator that source neuron is found
        //~ 
        //~ if (counter < N_RA_local)
            //~ if (Id_RA_local[counter] == i)
            //~ {
//~ 
              //~ source_neuron_found = true;
            //~ }
       //~ 
		//~ inp.read(reinterpret_cast<char *>(&spike_time), sizeof(spike_time));
//~ 
		//~ if (source_neuron_found)
		//~ {
			//~ last_dend_spike_time_local[counter] = spike_time;
			//~ counter++;
//~ 
		//~ }
//~ 
         //~ 
    //~ }
    //~ inp.close();
	//~ 
//~ }
//~ 
//~ void NetworkGrowthSimulator::read_num_bursts_in_recent_trials(const char* filename)
//~ {
    //~ std::ifstream inp;
//~ 
    //~ inp.open(filename, std::ios::in | std::ios::binary);
//~ 
    //~ // read number of HVC(RA) neurons
//~ 
	//~ int N;
//~ 
    //~ inp.read(reinterpret_cast<char *>(&N), sizeof(N));
//~ 
    //~ if (N != N_RA)
        //~ std::cerr << "Number of HVC(RA) neurons read from file with num_bursts_in_recent_trials: N = " << N 
				//~ << "is different from N_RA = " << N_RA << std::endl;
//~ 
	//~ int rate_window_long;
	//~ inp.read(reinterpret_cast<char *>(&rate_window_long), sizeof(rate_window_long));
//~ 
//~ 
	//~ if (rate_window_long != RATE_WINDOW_LONG)
        //~ std::cerr << "rate_window_long read from file with num_bursts_in_recent_trials: rate_window_long = " << rate_window_long 
				//~ << "is different from RATE_WINDOW_LONG = " << RATE_WINDOW_LONG << std::endl;
//~ 
//~ 
    //~ int num_trial;
    //~ 
    //~ inp.read(reinterpret_cast<char *>(&num_trial), sizeof(num_trial));
    //~ 
    //~ if (MPI_rank == 0)
        //~ std::cout << "Trial number read from file with num_bursts_in_recent_trials: " << num_trial << std::endl;
//~ 
    //~ int counter = 0; // counter of neuron id in the process
	//~ int num_bursts; // number of bursts produced in previoud trial j
//~ 
    //~ for (int i = 0; i < N_RA; i++)
    //~ {
        //~ bool source_neuron_found = false; // indicator that source neuron is found
        //~ 
        //~ if (counter < N_RA_local)
            //~ if (Id_RA_local[counter] == i)
            //~ {
//~ 
              //~ source_neuron_found = true;
            //~ }
       //~ 
        //~ for (int j = 0; j < rate_window_long; j++)
        //~ {            
            //~ inp.read(reinterpret_cast<char *>(&num_bursts), sizeof(num_bursts));
//~ 
            //~ if (source_neuron_found)
            //~ {
                //~ num_bursts_in_recent_trials[counter][j] = num_bursts;
            //~ }
//~ 
        //~ }
        //~ // increase counter
        //~ if (source_neuron_found)
            //~ counter++;
//~ 
    //~ }
    //~ inp.close();
//~ }

void NetworkGrowthSimulator::read_replacement_history(const char* filename)
{
	struct stat buf;
	
	if ( stat(filename, &buf) == 1 )
	{
		std::cerr << " File " << filename << " doesn't exist!\n" << std::endl;
		return;
	}
	
	std::ifstream in;

	in.open(filename, std::ios::in | std::ios::binary );
	
	int N, num_trial;
	
	in.read(reinterpret_cast<char *>(&N), sizeof(N));
	in.read(reinterpret_cast<char *>(&num_trial), sizeof(num_trial));
	
	if (N != N_RA)
		std::cerr << "Number of HVC-RA neurons in file with replacement history N = " << N 
				  << " is different from N_RA = " << N_RA << std::endl; 
	
			
	std::cout << "Trial number read from file with replacement history: " << num_trial << std::endl;
		  
	for (int i = 0; i < N; i++)
		in.read(reinterpret_cast<char *>(&num_trials_after_replacement_global[i]), sizeof(int));
		   
	in.close();
}

void NetworkGrowthSimulator::read_maturation_properties(const char* filename)
{
	struct stat buf;
			
	if ( stat(filename, &buf) == 1 )
	{
		std::cerr << " File " << filename << " doesn't exist!\n" << std::endl;
		return;
	}
	
	std::ifstream in;

	in.open(filename, std::ios::in | std::ios::binary);

	int N, num_trial;
	
	in.read(reinterpret_cast<char *>(&N), sizeof(int));
	in.read(reinterpret_cast<char *>(&num_trial), sizeof(int));

	if (N != N_RA)
		std::cerr << "Number of HVC-RA neurons in file with maturation properties N = " << N 
				  << " is different from N_RA = " << N_RA << std::endl; 
			
	std::cout << "Trial number read from file with maturation properties: " << num_trial << std::endl;	  
	
	for (int i = 0; i < N; i++){
		in.read(reinterpret_cast<char *>(&mature_global[i]), sizeof(int));
		in.read(reinterpret_cast<char *>(&maturation_scale_global[i]), sizeof(int));
		in.read(reinterpret_cast<char *>(&rest_potential_global[i]), sizeof(double));
		in.read(reinterpret_cast<char *>(&GCa_global[i]), sizeof(double));
	}
	
	in.close();
}

void NetworkGrowthSimulator::read_mature_indicators(const char* filename)
{
	struct stat buf;
			
	if ( stat(filename, &buf) == 1 )
	{
		std::cerr << " File " << filename << " doesn't exist!\n" << std::endl;
		return;
	}
	
	std::ifstream in;

	in.open(filename, std::ios::in | std::ios::binary);

	int N, num_trial;
	
	in.read(reinterpret_cast<char *>(&N), sizeof(int));
	in.read(reinterpret_cast<char *>(&num_trial), sizeof(int));

	if (N != N_RA)
		std::cerr << "Number of HVC-RA neurons in file with mature indicators N = " << N 
				  << " is different from N_RA = " << N_RA << std::endl; 
			
	std::cout << "Trial number read from file with mature indicators: " << num_trial << std::endl;	  
	
	for (int i = 0; i < N; i++)
		in.read(reinterpret_cast<char *>(&mature_global[i]), sizeof(int));

	
	in.close();
}

void NetworkGrowthSimulator::read_remodeled_indicators(const char* filename)
{
	struct stat buf;
			
	if ( stat(filename, &buf) == 1 )
	{
		std::cerr << " File " << filename << " doesn't exist!\n" << std::endl;
		return;
	}
	
	std::ifstream in;

	in.open(filename, std::ios::in | std::ios::binary);

	int N, num_trial;
	
	in.read(reinterpret_cast<char *>(&N), sizeof(int));
	in.read(reinterpret_cast<char *>(&num_trial), sizeof(int));

	if (N != N_RA)
		std::cerr << "Number of HVC-RA neurons in file with remodeled indicators N = " << N 
				  << " is different from N_RA = " << N_RA << std::endl; 
			
	std::cout << "Trial number read from file with remodeled indicators: " << num_trial << std::endl;	  
	
	for (int i = 0; i < N; i++)
		in.read(reinterpret_cast<char *>(&remodeled_global[i]), sizeof(int));

	
	in.close();
}

void NetworkGrowthSimulator::read_activity_history(const char* filename)
{
	struct stat buf;
			
	if ( stat(filename, &buf) == 1 )
	{
		std::cerr << " File " << filename << " doesn't exist!\n" << std::endl;
		return;
	}
	
	std::ifstream in;

	in.open(filename, std::ios::in | std::ios::binary);

	int N, rate_window, num_trial;

	in.read(reinterpret_cast<char *>(&N), sizeof(int));
	in.read(reinterpret_cast<char *>(&rate_window), sizeof(int));
	in.read(reinterpret_cast<char *>(&num_trial), sizeof(int));
	
	if (N != N_RA)
		std::cerr << "Number of HVC-RA neurons in file with activity history N = " << N 
				  << " is different from N_RA = " << N_RA << std::endl; 
	
	if (rate_window != maturation_params.RATE_WINDOW_LONG)
		std::cerr << "Rate window in file with activity history rate_window = " << rate_window 
				  << " is different from RATE_WINDOW_LONG = " << maturation_params.RATE_WINDOW_LONG << std::endl; 
	
			
	std::cout << "Trial number read from file with activity history: " << num_trial << std::endl;	  
	
	
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < rate_window; j++)
		{
			//printf("weigths[%d][%d] = %1.10f\n", i, j, weights[i][j]);
			in.read(reinterpret_cast<char *>(&num_spikes_in_recent_trials_global[i][j]), sizeof(int));

		}
	}
	
	in.close();
}

void NetworkGrowthSimulator::read_super_synapses(const char* filename)
{
	struct stat buf;
	
	if ( stat(filename, &buf) == 1 )
	{
		std::cerr << " File " << filename << " doesn't exist!\n" << std::endl;
		return;
	}
	
	std::ifstream in;

	in.open(filename, std::ios::in | std::ios::binary);

	// write number of neurons to each file
	int N, num_trial;

	in.read(reinterpret_cast<char *>(&N), sizeof(int));
	in.read(reinterpret_cast<char *>(&num_trial), sizeof(int));

	if (N != N_RA)
		std::cerr << "Number of HVC-RA neurons in file with super synapses N = " << N 
				  << " is different from N_RA = " << N_RA << std::endl; 
	
	std::cout << "Trial number read from file with super synapses: " << num_trial << std::endl;	  
	
	for (int i = 0; i < N; i++)
	{
		// write number of connections to RA neurons
		int num_synapses;
		
		in.read(reinterpret_cast<char *>(&num_synapses), sizeof(int));
	
		supersynapses_global[i].resize(num_synapses);
		
		for (int j = 0; j < num_synapses; j++)
			in.read(reinterpret_cast<char *>(&supersynapses_global[i][j]), sizeof(int));
	}
	in.close();
}


void NetworkGrowthSimulator::read_active_synapses(const char* filename)
{
	struct stat buf;
			
	if ( stat(filename, &buf) == 1 )
	{
		std::cerr << " File " << filename << " doesn't exist!\n" << std::endl;
		return;
	}
	
	std::ifstream in;
	in.open(filename, std::ios::in | std::ios::binary);

	// write number of neurons to each file
	int N, num_trial;

	in.read(reinterpret_cast<char *>(&N), sizeof(int));
	in.read(reinterpret_cast<char *>(&num_trial), sizeof(int));


	if (N != N_RA)
		std::cerr << "Number of HVC-RA neurons in file with active synapses N = " << N 
				  << " is different from N_RA = " << N_RA << std::endl; 

	std::cout << "Trial number read from file with active synapses: " << num_trial << std::endl;	  
	
	for (int i = 0; i < N; i++)
	{
		// write number of connections to RA neurons
		int num_synapses;
		
		in.read(reinterpret_cast<char *>(&num_synapses), sizeof(int));
	
		active_synapses_global[i].resize(num_synapses);
		
		for (int j = 0; j < num_synapses; j++)
			in.read(reinterpret_cast<char *>(&active_synapses_global[i][j]), sizeof(int));
	}
	in.close();
}

void NetworkGrowthSimulator::read_axonal_delays(std::vector<std::vector<double>> &axonal_delays, const char* filename)
{
	struct stat buf;
	
	if ( stat(filename, &buf) == 1 )
	{
		std::cerr << " File " << filename << " doesn't exist!\n" << std::endl;
		return;
	}
	
	std::ifstream in;

	in.open(filename, std::ios::in | std::ios::binary);

	
	int num_neurons = static_cast<int>(axonal_delays.size());

	int N, num_trial;

	in.read(reinterpret_cast<char *>(&N), sizeof(int));
	in.read(reinterpret_cast<char *>(&num_trial), sizeof(int));
	
	
	
	if (N != num_neurons)
		std::cerr << "Number of neurons in file " << filename << " N = " << N 
					<< " is different from size of axonal delay array = " << num_neurons << std::endl;
	
	std::cout << "Trial number read from file with axonal delays: " << num_trial << std::endl;


	for (int i = 0; i < N; i++)
	{
		int num_targets;
		
		in.read(reinterpret_cast<char *>(&num_targets), sizeof(int));
		
		axonal_delays[i].resize(num_targets);
		
		for (int j = 0; j < num_targets; j++)
		{
			//printf("weigths[%d][%d] = %1.10f\n", i, j, weights[i][j]);
			in.read(reinterpret_cast<char *>(&axonal_delays[i][j]), sizeof(double));

		}
	}
	in.close();
}

void NetworkGrowthSimulator::read_weights(const char* filename)
{
	struct stat buf;
	
	if ( stat(filename, &buf) == 1 )
	{
		std::cerr << " File " << filename << " doesn't exist!\n" << std::endl;
		return;
	}
	
	std::ifstream in;

	in.open(filename, std::ios::in | std::ios::binary);

	int N, num_trial;

	in.read(reinterpret_cast<char *>(&N), sizeof(int));
	in.read(reinterpret_cast<char *>(&num_trial), sizeof(int));
	
	if (N != N_RA)
		std::cerr << "Number of HVC-RA neurons in file with weights N = " << N << " is different from N_RA = " << N_RA << std::endl;
	
	std::cout << "Trial number read from file with synaptic weights: " << num_trial << std::endl;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			in.read(reinterpret_cast<char *>(&weights_RA_RA_global[i][j]), sizeof(double));
	}
	
	in.close();
}

//~ void NetworkGrowthSimulator::read_weights(const char* filename)
//~ {
    //~ std::ifstream inp;
//~ 
    //~ inp.open(filename, std::ios::in | std::ios::binary);
//~ 
    //~ // read number of HVC(RA) neurons
//~ 
  //~ int N;
//~ 
    //~ inp.read(reinterpret_cast<char *>(&N), sizeof(N));
//~ 
    //~ if (N != N_RA)
        //~ std::cerr << "Number of HVC(RA) neurons read from file with supersynapses: N = " << N << "is different from N_RA = " << N_RA << std::endl;
//~ 
    //~ int num_trial;
    //~ 
    //~ inp.read(reinterpret_cast<char *>(&num_trial), sizeof(num_trial));
    //~ 
    //~ if (MPI_rank == 0)
        //~ std::cout << "Trial number read from file with synaptic weights: " << num_trial << std::endl;
//~ 
    //~ int counter = 0; // counter of neuron id in the process
//~ 
//~ 
//~ 
    //~ for (int i = 0; i < N_RA; i++)
    //~ {
        //~ bool source_neuron_found = false; // indicator that source neuron is found
        //~ 
        //~ if (counter < N_RA_local)
            //~ if (Id_RA_local[counter] == i)
            //~ {
//~ 
              //~ source_neuron_found = true;
            //~ }
       //~ 
        //~ for (int j = 0; j < N_RA; j++)
        //~ {
            //~ 
            //~ double synapse_weight; // weight of synapse
//~ 
            //~ inp.read(reinterpret_cast<char *>(&synapse_weight), sizeof(synapse_weight));
//~ 
            //~ if (source_neuron_found)
            //~ {
                //~ weights_local[counter][j] = synapse_weight;
            //~ }
//~ 
        //~ }
        //~ // increase counter
        //~ if (source_neuron_found)
            //~ counter++;
//~ 
    //~ }
    //~ inp.close();
//~ }
//~ 
//~ void NetworkGrowthSimulator::read_super_synapses(const char *filename)
//~ {
    //~ std::ifstream inp;
//~ 
    //~ inp.open(filename, std::ios::in | std::ios::binary);
//~ 
    //~ // read number of HVC(RA) neurons
//~ 
  //~ int N;
//~ 
    //~ inp.read(reinterpret_cast<char *>(&N), sizeof(N));
//~ 
    //~ if (N != N_RA)
        //~ std::cerr << "Number of HVC(RA) neurons read from file with supersynapses: N = " << N << "is different from N_RA = " << N_RA << std::endl;
//~ 
    //~ int counter = 0; // counter of neuron id in the process
//~ 
    //~ for (int i = 0; i < N_RA; i++)
    //~ {
        //~ // read neuron's ID number
        //~ int n_id;
        //~ bool source_neuron_found = false; // indicator that source neuron is found
        //~ 
        //~ inp.read(reinterpret_cast<char *>(&n_id), sizeof(n_id));
        //~ 
        //~ // read number of connections to RA neurons
//~ 
      //~ int num_super;
        //~ 
        //~ inp.read(reinterpret_cast<char *>(&num_super), sizeof(num_super));
        //~ 
        //~ if (counter < N_RA_local)
            //~ if (Id_RA_local[counter] == i)
            //~ {
                //~ supersynapses_local[counter].resize(num_super);
//~ 
              //~ source_neuron_found = true;
            //~ }
        //~ // read supersynapses
        //~ for (int j = 0; j < num_super; j++)
        //~ {
            //~ int target_id; // id of supersynaptic target
            //~ double synapse_weight; // weight of synapse
//~ 
            //~ inp.read(reinterpret_cast<char *>(&target_id), sizeof(target_id));
            //~ inp.read(reinterpret_cast<char *>(&synapse_weight), sizeof(synapse_weight));
//~ 
            //~ if (source_neuron_found)
            //~ {
                //~ supersynapses_local[counter][j] = target_id;
                //~ supersynapses_indicators_local[counter][target_id] = 1;
            //~ }
//~ 
        //~ }
        //~ // increase counter
        //~ if (source_neuron_found)
            //~ counter++;
//~ 
    //~ }
    //~ inp.close();
//~ }


void NetworkGrowthSimulator::send_growth_parameters()
{
	// create MPI datatypes for each structure
	
	//////////////////////////
	//// Synaptic parameters
	//////////////////////////
	
	const int countSynapticParameters = 13;
    
    MPI_Datatype MpiSynapticParameters;
    MPI_Datatype typeSynapticParameters[countSynapticParameters] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
												  MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
												  MPI_DOUBLE};
												  
    int blocklenSynapticParameters[countSynapticParameters] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    
    MPI_Aint dispSynapticParameters[countSynapticParameters] = { offsetof(SynapticParameters, Nss), offsetof(SynapticParameters, A_P),
											   offsetof(SynapticParameters, A_D), offsetof(SynapticParameters, T_0),
											   offsetof(SynapticParameters, T_P), offsetof(SynapticParameters, T_D),
											   offsetof(SynapticParameters, TAU_P), offsetof(SynapticParameters, TAU_D),
											   offsetof(SynapticParameters, BETA), offsetof(SynapticParameters, BETA_SUPERSYNAPSE),
											   offsetof(SynapticParameters, ACTIVATION_THRESHOLD), offsetof(SynapticParameters, SUPERSYNAPSE_THRESHOLD),
											   offsetof(SynapticParameters, WEIGHT_MAX)};

    MPI_Type_create_struct(countSynapticParameters, blocklenSynapticParameters, dispSynapticParameters, 
									typeSynapticParameters, &MpiSynapticParameters);
    MPI_Type_commit(&MpiSynapticParameters);
	
	//////////////////////////
	//// Maturation parameters
	//////////////////////////
	const int countMaturationParameters = 22;
    
    MPI_Datatype MpiMaturationParameters;
    MPI_Datatype typeMaturationParameters[countMaturationParameters] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
												  MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
												  MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
												  MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
												  
    int blocklenMaturationParameters[countMaturationParameters] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    
    MPI_Aint dispMaturationParameters[countMaturationParameters] = { offsetof(MaturationParameters, E_REST_MATURE), offsetof(MaturationParameters, E_REST_IMMATURE),
																	 offsetof(MaturationParameters, E_GABA_MATURE), offsetof(MaturationParameters, E_GABA_IMMATURE),
																	 offsetof(MaturationParameters, AD_MATURE),     offsetof(MaturationParameters, AD_IMMATURE),
																	 offsetof(MaturationParameters, GK_MATURE), 	offsetof(MaturationParameters, GK_IMMATURE),
																	 offsetof(MaturationParameters, GNA_MATURE), 	offsetof(MaturationParameters, GNA_IMMATURE),
																	 offsetof(MaturationParameters, GCA_MATURE), 	offsetof(MaturationParameters, GCA_IMMATURE),
																	 offsetof(MaturationParameters, GCAK_MATURE), 	offsetof(MaturationParameters, GCAK_IMMATURE),
																	 offsetof(MaturationParameters, GSL_MATURE), 	offsetof(MaturationParameters, GSL_IMMATURE),
																	 offsetof(MaturationParameters, GDL_MATURE), 	offsetof(MaturationParameters, GDL_IMMATURE),
																	 offsetof(MaturationParameters, RC_MATURE), 	offsetof(MaturationParameters, RC_IMMATURE),
																	 offsetof(MaturationParameters, DEATH_RATE_THRESHOLD), 	offsetof(MaturationParameters, RATE_WINDOW_LONG) };
											  

    MPI_Type_create_struct(countMaturationParameters, blocklenMaturationParameters, dispMaturationParameters, 
									typeMaturationParameters, &MpiMaturationParameters);
    MPI_Type_commit(&MpiMaturationParameters);
    
    
    //////////////////////////
	//// Noise parameters
	//////////////////////////
	const int countNoiseParameters = 4;
    
    MPI_Datatype MpiNoiseParameters;
    MPI_Datatype typeNoiseParameters[countNoiseParameters] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
												  
    int blocklenNoiseParameters[countNoiseParameters] = { 1, 1, 1, 1 };
    
    MPI_Aint dispNoiseParameters[countNoiseParameters] = { offsetof(NoiseParameters, white_noise_mean_soma), offsetof(NoiseParameters, white_noise_std_soma),
														   offsetof(NoiseParameters, white_noise_mean_dend), offsetof(NoiseParameters, white_noise_std_dend) };

    MPI_Type_create_struct(countNoiseParameters, blocklenNoiseParameters, dispNoiseParameters, 
									typeNoiseParameters, &MpiNoiseParameters);
    MPI_Type_commit(&MpiNoiseParameters);
    
    ///////////////////////////////////////
    //// Send parameters to every process
	///////////////////////////////////////
	
	MPI_Bcast(&synaptic_params,   1, MpiSynapticParameters,   0, MPI_COMM_WORLD);
	MPI_Bcast(&maturation_params, 1, MpiMaturationParameters, 0, MPI_COMM_WORLD);
	MPI_Bcast(&noise_params,      1, MpiNoiseParameters,      0, MPI_COMM_WORLD);

	//if (MPI_rank == 1)
		//this->print_simulation_parameters();
}

void NetworkGrowthSimulator::send_axonal_delays_RA2RA()
{
	int* sendcounts;
	int* displs;
	
	if (MPI_rank == 0)
	{
		sendcounts = new int[MPI_size];
		displs = new int[MPI_size];
		
		sendcounts[0] = N_RA_sizes[0];
		displs[0] = 0;
		
		for (int i = 1; i < MPI_size; i++)
		{
			sendcounts[i] = N_RA_sizes[i];
			displs[i] = displs[i-1] + sendcounts[i-1];
		}
    }

	// for rank 0 copy data from global arrays to local
	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_RA_local; i++)
			axonal_delays_RA_RA_local[i] = axonal_delays_RA_RA_global[i];
	}

	MPI_Status status;
	int start_row_id = N_RA_sizes[0];

	for (int i = 1; i < MPI_size; i++)
	{ 
		for (int j = 0; j < N_RA_sizes[i]; j++)
		{ 	
			if (MPI_rank == 0)
				MPI_Send(&axonal_delays_RA_RA_global[j + start_row_id][0], N_RA, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			
			else if (MPI_rank == i)
				MPI_Recv(&axonal_delays_RA_RA_local[j][0], N_RA, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		}
		
		start_row_id += N_RA_sizes[i];
	}

	if (MPI_rank == 0)
	{
        delete [] sendcounts;
        delete [] displs;
    }
}


void NetworkGrowthSimulator::send_connections_RA2RA()
{
	int* sendcounts;
	int* displs;
	
	if (MPI_rank == 0)
	{
		sendcounts = new int[MPI_size];
		displs = new int[MPI_size];
		
		sendcounts[0] = N_RA_sizes[0];
		displs[0] = 0;
		
		for (int i = 1; i < MPI_size; i++)
		{
			sendcounts[i] = N_RA_sizes[i];
			displs[i] = displs[i-1] + sendcounts[i-1];
		}
    }

	// for rank 0 copy data from global arrays to local
	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_RA_local; i++)
		{
			weights_RA_RA_local[i] = weights_RA_RA_global[i];
			axonal_delays_RA_RA_local[i] = axonal_delays_RA_RA_global[i];
		}
	}

	MPI_Status status;
	int start_row_id = N_RA_sizes[0];

	for (int i = 1; i < MPI_size; i++)
	{ 
		for (int j = 0; j < N_RA_sizes[i]; j++)
		{ 	
			if (MPI_rank == 0)
			{
				MPI_Send(&weights_RA_RA_global[j + start_row_id][0], N_RA, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				MPI_Send(&axonal_delays_RA_RA_global[j + start_row_id][0], N_RA, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}	
			else if (MPI_rank == i)
			{
				MPI_Recv(&weights_RA_RA_local[j][0], N_RA, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(&axonal_delays_RA_RA_local[j][0], N_RA, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
			}
		}
		
		start_row_id += N_RA_sizes[i];
	}

	if (MPI_rank == 0)
	{
        delete [] sendcounts;
        delete [] displs;
    }
}

void NetworkGrowthSimulator::send_connections_RAandI()
{
	MPI_Status status;
	
	int* sendcounts_syn_num_RA;
	int* sendcounts_syn_num_I;
	int* displs_syn_num_RA;
	int* displs_syn_num_I;
	int* syn_num_RA;
	int* syn_num_I;
	int* syn_num_RA_local = new int[N_RA_local];
	int* syn_num_I_local = new int[N_I_local];

	if (MPI_rank == 0)
	{
		sendcounts_syn_num_RA = new int[MPI_size];
		sendcounts_syn_num_I = new int[MPI_size];
		syn_num_RA = new int[N_RA];
		syn_num_I = new int[N_I];

		displs_syn_num_RA = new int[MPI_size];
		displs_syn_num_I = new int[MPI_size];

		// send connections to all processes

		sendcounts_syn_num_RA[0] = N_RA_sizes[0];
		sendcounts_syn_num_I[0] = N_I_sizes[0];
		displs_syn_num_RA[0] = 0;
		displs_syn_num_I[0] = 0;

		for (int i = 1; i < MPI_size; i++)
		{
			sendcounts_syn_num_RA[i] = N_RA_sizes[i];
			sendcounts_syn_num_I[i] = N_I_sizes[i];
			displs_syn_num_RA[i] = displs_syn_num_RA[i-1] + sendcounts_syn_num_RA[i-1];
			displs_syn_num_I[i] = displs_syn_num_I[i-1] + sendcounts_syn_num_I[i-1];
		}

		for (int i = 0; i < N_RA; i++)
		{
			syn_num_RA[i] = static_cast<int>(syn_ID_RA_I_global[i].size());
			//printf("Master process; syn_num_RA[%d] = %d\n", i, syn_num_RA[i]);
		}

		for (int i = 0; i < N_I; i++)
		{
			syn_num_I[i] = static_cast<int>(syn_ID_I_RA_global[i].size());

			//printf("Master process; syn_num_I[%d] = %d\n", i, syn_num_I[i]);
		}


    }

	// send number of connections for each neuron

	MPI_Scatterv(&syn_num_RA[0], sendcounts_syn_num_RA, displs_syn_num_RA, MPI_INT,
		&syn_num_RA_local[0], N_RA_local, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Scatterv(&syn_num_I[0], sendcounts_syn_num_I, displs_syn_num_I, MPI_INT,
		&syn_num_I_local[0], N_I_local, MPI_INT, 0, MPI_COMM_WORLD);

	//for (int i = 0; i < N_RA_local; i++)
	//	printf("My rank = %d; sun_num_RA[%d] = %d\n", MPI_rank, Id_RA_local[i], syn_num_RA_local[i]);
	
    //for (int i = 0; i < N_I_local; i++)
	//	printf("My rank = %d; sun_num_I[%d] = %d\n", MPI_rank, Id_I_local[i], syn_num_I_local[i]);

	for (int i = 0; i < N_RA_local; i++)
	{
		syn_ID_RA_I_local[i].resize(syn_num_RA_local[i]);
		weights_RA_I_local[i].resize(syn_num_RA_local[i]);
		axonal_delays_RA_I_local[i].resize(syn_num_RA_local[i]);
	}

	for (int i = 0; i < N_I_local; i++)
	{
		syn_ID_I_RA_local[i].resize(syn_num_I_local[i]);
		weights_I_RA_local[i].resize(syn_num_I_local[i]);
		axonal_delays_I_RA_local[i].resize(syn_num_I_local[i]);
	}

	// send connection ID and weights

	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_RA_local; i++)
		{
			weights_RA_I_local[i] = weights_RA_I_global[i];
			syn_ID_RA_I_local[i] = syn_ID_RA_I_global[i];
			axonal_delays_RA_I_local[i] = axonal_delays_RA_I_global[i];
		}

		for (int i = 0; i < N_I_local; i++)
		{
			weights_I_RA_local[i] = weights_I_RA_global[i];
			syn_ID_I_RA_local[i] = syn_ID_I_RA_global[i];
			axonal_delays_I_RA_local[i] = axonal_delays_I_RA_global[i];
		}

		int offset_RA = N_RA_local;
		int offset_I = N_I_local;

		for (int i = 1; i < MPI_size; i++)
		{
            //  send ID and weights of RA targets
			for (int j = 0; j < sendcounts_syn_num_RA[i]; j++)
			{
				//printf("Master. Here!\n");
				MPI_Send(&syn_ID_RA_I_global[offset_RA+j][0], syn_num_RA[offset_RA+j], MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&weights_RA_I_global[offset_RA+j][0], syn_num_RA[offset_RA+j], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				MPI_Send(&axonal_delays_RA_I_global[offset_RA+j][0], syn_num_RA[offset_RA+j], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				
				//for (int k = 0; k < syn_num_RA[offset_RA+j]; k++)
				//	printf("Master. RA neuron %d; RA2I %d; w = %f\n", offset_RA+j, syn_ID_RA_I_global[offset_RA+j][k],
                 //       weights_RA_I_global[offset_RA+j][k]);

			}
			offset_RA += sendcounts_syn_num_RA[i];

            // send ID and weights of I targets
			for (int j = 0; j < sendcounts_syn_num_I[i]; j++)
			{
				//printf("Master. Here!\n");
				MPI_Send(&syn_ID_I_RA_global[offset_I+j][0], syn_num_I[offset_I+j], MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&weights_I_RA_global[offset_I+j][0], syn_num_I[offset_I+j], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				MPI_Send(&axonal_delays_I_RA_global[offset_I+j][0], syn_num_I[offset_I+j], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

				//for (int k = 0; k < syn_num_I[offset_I+j]; k++)
				//	printf("Master. I neuron %d; I2RA %d; w = %f\n", offset_I+j, syn_ID_I_RA_global[offset_I+j][k],
                  //      weights_I_RA_global[offset_I+j][k]);

			}
			offset_I += sendcounts_syn_num_I[i];
			//printf("Master. Number of synapses;\n");

		}

	}
	else
	{
        // receive ID and weights of RA targets
		for (int i = 0; i < N_RA_local; i++)
		{
			//printf("Rank = %d. Here!\n", MPI_rank);
			MPI_Recv(&syn_ID_RA_I_local[i][0], syn_num_RA_local[i], MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&weights_RA_I_local[i][0], syn_num_RA_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&axonal_delays_RA_I_local[i][0], syn_num_RA_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

			//for (int j = 0; j < syn_num_RA_local[i]; j++)
			//	printf("My rank = %d; RA neuron %d; RA2I %d; w = %f\n", MPI_rank, Id_RA_local[i], syn_ID_RA_I_local[i][j],
              //      weights_RA_I_local[i][j]);
        }

        // receive ID and weights of I targets
        for (int i = 0; i < N_I_local; i++)
        {
            MPI_Recv(&syn_ID_I_RA_local[i][0], syn_num_I_local[i], MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&weights_I_RA_local[i][0], syn_num_I_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&axonal_delays_I_RA_local[i][0], syn_num_I_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

			//for (int j = 0; j < syn_num_I_local[i]; j++)
			//	printf("My rank = %d; I neuron %d; I2RA %d; w = %f\n", MPI_rank, Id_I_local[i], syn_ID_I_RA_local[i][j],
              //      weights_I_RA_local[i][j]);

		}
	}

	if (MPI_rank == 0)
	{
        delete [] sendcounts_syn_num_RA;
        delete [] sendcounts_syn_num_I;
        delete [] syn_num_RA;
        delete [] syn_num_I;
        delete [] displs_syn_num_RA;
        delete [] displs_syn_num_I;
    }
    
     delete [] syn_num_RA_local;
     delete [] syn_num_I_local;
}

void NetworkGrowthSimulator::send_maturation_properties()
{
	int* sendcounts;
	int* displs;
	
	if (MPI_rank == 0)
	{
		sendcounts = new int[MPI_size];
		displs = new int[MPI_size];
		
		sendcounts[0] = N_RA_sizes[0];
		displs[0] = 0;
		
		for (int i = 1; i < MPI_size; i++)
		{
			sendcounts[i] = N_RA_sizes[i];
			displs[i] = displs[i-1] + sendcounts[i-1];
		}
    }
	
	MPI_Scatterv(&maturation_scale_global[0], sendcounts, displs, MPI_INT,
		&maturation_scale_local[0], N_RA_local, MPI_INT, 0, MPI_COMM_WORLD);
	
	MPI_Scatterv(&mature_global[0], sendcounts, displs, MPI_INT,
		&mature_local[0], N_RA_local, MPI_INT, 0, MPI_COMM_WORLD);
	
	MPI_Scatterv(&GCa_global[0], sendcounts, displs, MPI_DOUBLE,
		&GCa_local[0], N_RA_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	MPI_Scatterv(&rest_potential_global[0], sendcounts, displs, MPI_DOUBLE,
		&rest_potential_local[0], N_RA_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if (MPI_rank == 0)
	{
		delete [] sendcounts;
		delete [] displs;
	}
}


void NetworkGrowthSimulator::send_mature_indicators()
{
	MPI_Bcast(&mature_global[0], N_RA, MPI_INT, 0, MPI_COMM_WORLD);
}

void NetworkGrowthSimulator::send_remodeled_indicators()
{
	int* sendcounts;
	int* displs;
	
	if (MPI_rank == 0)
	{
		sendcounts = new int[MPI_size];
		displs = new int[MPI_size];
		
		sendcounts[0] = N_RA_sizes[0];
		displs[0] = 0;
		
		for (int i = 1; i < MPI_size; i++)
		{
			sendcounts[i] = N_RA_sizes[i];
			displs[i] = displs[i-1] + sendcounts[i-1];
		}
    }
	
	MPI_Scatterv(&remodeled_global[0], sendcounts, displs, MPI_INT,
		&remodeled_local[0], N_RA_local, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (MPI_rank == 0)
	{
		delete [] sendcounts;
		delete [] displs;
	}
}

void NetworkGrowthSimulator::send_active_synapses()
{
	int* sendcounts;
	int* displs;
	
	int* syn_num;
	int* syn_num_local = new int[N_RA_local];

	if (MPI_rank == 0)
	{
		sendcounts = new int[MPI_size];
		syn_num = new int[N_RA];
		displs = new int[MPI_size];
		
		sendcounts[0] = N_RA_sizes[0];
		displs[0] = 0;
	
		for (int i = 1; i < MPI_size; i++)
		{
			sendcounts[i] = N_RA_sizes[i];
			displs[i] = displs[i-1] + sendcounts[i-1];
		}

		for (int i = 0; i < N_RA; i++)
			syn_num[i] = static_cast<int>(active_synapses_global[i].size());
    }

	// send number of connections for each neuron

	MPI_Scatterv(&syn_num[0], sendcounts, displs, MPI_INT,
		&syn_num_local[0], N_RA_local, MPI_INT, 0, MPI_COMM_WORLD);
	
	
	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_RA_local; i++)
			active_synapses_local[i] = active_synapses_global[i];
	
		int offset = N_RA_local;
		
		for (int i = 1; i < MPI_size; i++)
		{
			for (int j = 0; j < N_RA_sizes[i]; j++)
				MPI_Send(&active_synapses_global[offset+j][0], syn_num[offset + j], MPI_INT, i, j, MPI_COMM_WORLD);
			
			offset += N_RA_sizes[i];
		}
	}
	else
	{
		for (int i = 0; i < N_RA_local; i++)
			active_synapses_local[i].resize(syn_num_local[i]);
	
		MPI_Status status;
	
		for (int i = 0; i < N_RA_local; i++)
			MPI_Recv(&active_synapses_local[i][0], syn_num_local[i], MPI_INT, 0, i, MPI_COMM_WORLD, &status);
	}
	
	// set active synapse indicators
	
	
	//~ if (MPI_rank == 1)
	//~ {
		//~ for (int i = 0; i < N_RA_local; i++)
		//~ {
			//~ std::cout << "Neuron " << Id_RA_local[i] << "\n"; 
			//~ 
			//~ for (size_t j = 0; j < active_synapses_local[i].size(); j++)
				//~ std::cout << active_synapses_local[i][j] << " ";
			//~ std::cout << std::endl;
		//~ }
	//~ }
	
	
	
	if (MPI_rank == 0)
	{
		delete [] sendcounts;
		delete [] displs;
		delete [] syn_num;
	}
	
	delete [] syn_num_local;
}

void NetworkGrowthSimulator::send_super_synapses()
{
	int* sendcounts;
	int* displs;
	
	int* syn_num;
	int* syn_num_local = new int[N_RA_local];

	if (MPI_rank == 0)
	{
		sendcounts = new int[MPI_size];
		syn_num = new int[N_RA];
		displs = new int[MPI_size];
		
		sendcounts[0] = N_RA_sizes[0];
		displs[0] = 0;
	
		for (int i = 1; i < MPI_size; i++)
		{
			sendcounts[i] = N_RA_sizes[i];
			displs[i] = displs[i-1] + sendcounts[i-1];
		}

		for (int i = 0; i < N_RA; i++)
			syn_num[i] = static_cast<int>(supersynapses_global[i].size());
    }

	// send number of connections for each neuron

	MPI_Scatterv(&syn_num[0], sendcounts, displs, MPI_INT,
		&syn_num_local[0], N_RA_local, MPI_INT, 0, MPI_COMM_WORLD);
	
	
	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_RA_local; i++)
			supersynapses_local[i] = supersynapses_global[i];
	
		int offset = N_RA_local;
		
		for (int i = 1; i < MPI_size; i++)
		{
			for (int j = 0; j < N_RA_sizes[i]; j++)
				MPI_Send(&supersynapses_global[offset+j][0], syn_num[offset + j], MPI_INT, i, j, MPI_COMM_WORLD);
			
			offset += N_RA_sizes[i];
		}
	}
	else
	{
		for (int i = 0; i < N_RA_local; i++)
			supersynapses_local[i].resize(syn_num_local[i]);
	
		MPI_Status status;
	
		for (int i = 0; i < N_RA_local; i++)
			MPI_Recv(&supersynapses_local[i][0], syn_num_local[i], MPI_INT, 0, i, MPI_COMM_WORLD, &status);
	}
	
	if (MPI_rank == 0)
	{
		delete [] sendcounts;
		delete [] displs;
		delete [] syn_num;
	}
	
	delete [] syn_num_local;
}

void NetworkGrowthSimulator::send_activity_history()
{
	int* sendcounts;
	int* displs;
	
	if (MPI_rank == 0)
	{
		sendcounts = new int[MPI_size];
		displs = new int[MPI_size];
		
		sendcounts[0] = N_RA_sizes[0];
		displs[0] = 0;
		
		for (int i = 1; i < MPI_size; i++)
		{
			sendcounts[i] = N_RA_sizes[i];
			displs[i] = displs[i-1] + sendcounts[i-1];
		}
    }
	
	// for rank 0 copy data from global arrays to local
	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_RA_local; i++)
			std::copy(num_spikes_in_recent_trials_global[i].begin(), num_spikes_in_recent_trials_global[i].end(), 
							num_spikes_in_recent_trials_local[i].begin());	
	}

	MPI_Status status;
	int start_row_id = N_RA_sizes[0];

	for (int i = 1; i < MPI_size; i++)
	{ 
		for (int j = 0; j < N_RA_sizes[i]; j++)
		{ 	
			if (MPI_rank == 0)
				MPI_Send(&num_spikes_in_recent_trials_global[j + start_row_id][0], maturation_params.RATE_WINDOW_LONG, MPI_INT, i, 0, MPI_COMM_WORLD);
			
			else if (MPI_rank == i)
				MPI_Recv(&num_spikes_in_recent_trials_local[j][0], maturation_params.RATE_WINDOW_LONG, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
				
		}
		
		start_row_id += N_RA_sizes[i];
	}
	
	if (MPI_rank == 0)
	{
        delete [] sendcounts;
        delete [] displs;
    }
}


void NetworkGrowthSimulator::send_replacement_history()
{
	int* sendcounts;
	int* displs;
	
	if (MPI_rank == 0)
	{
		sendcounts = new int[MPI_size];
		displs = new int[MPI_size];
		
		sendcounts[0] = N_RA_sizes[0];
		displs[0] = 0;
		
		for (int i = 1; i < MPI_size; i++)
		{
			sendcounts[i] = N_RA_sizes[i];
			displs[i] = displs[i-1] + sendcounts[i-1];
		}
    }
    
	MPI_Scatterv(&num_trials_after_replacement_global[0], sendcounts, displs, MPI_INT,
		&num_trials_after_replacement_local[0], N_RA_local, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (MPI_rank == 0)
	{
        delete [] sendcounts;
        delete [] displs;
    }
}

void NetworkGrowthSimulator::disable_RA2I_immature_outputs()
{
	for (int i = 0; i < N_RA_local; i++)
	{
		if ( mature_global[Id_RA_local[i]] == 0)
			for (size_t j = 0; j < syn_ID_RA_I_local[i].size(); j++)
				weights_RA_I_local[i][j] = 0.0;
		
	}
}

//~ void NetworkGrowthSimulator::disable_RA2I_connections()
//~ {
	//~ for (int i = 0; i < N_RA_local; i++)
	//~ {
		//~ weights_RA_I_local[i].clear();
		//~ syn_ID_RA_I_local[i].clear();
	//~ }
//~ }
//~ 
//~ void NetworkGrowthSimulator::initialize_ideal_chain_connections(int num_layers)
//~ {
    //~ if (MPI_rank == 0)
    //~ {
        //~ // connections from HVC(RA) to HVC(I) and vice versa
        //~ 
        //~ for (int i = 0; i < num_layers; i++)
        //~ {
            //~ // from HVC(RA) to HVC(I)
            //~ for (int j = 0; j < Nss; j++)
            //~ {
                //~ double G = this->sample_Ge2i();
//~ 
                //~ weights_RA_I_global[i*Nss + j].push_back(G);
                //~ syn_ID_RA_I_global[i*Nss + j].push_back(i*Nss + j);
            //~ }
        //~ 
            //~ // from HVC(I) to HVC(RA)
            //~ for (int j = 0; j < Nss; j++)
            //~ {
                //~ double G = this->sample_Gi2e();
//~ 
                //~ weights_I_RA_global[i*Nss + j].push_back(G);
                //~ syn_ID_I_RA_global[i*Nss + j].push_back((i+1)*Nss + j);
//~ 
            //~ }
        //~ }
    //~ }
//~ 
    //~ this->send_connections();
//~ 
    //~ std::string fileRA2I = outputDirectory + "RA_I_connections.bin";
    //~ std::string fileI2RA = outputDirectory + "I_RA_connections.bin";
    //~ std::string filePajekFixed = outputDirectory + "fixed.net";
//~ 
    //~ this->write_invariable_synapses(fileRA2I.c_str(), fileI2RA.c_str());
    //~ this->write_pajek_fixed(filePajekFixed.c_str());
    //~ 
    //~ // initialize excitatory connections for all chain layers
    //~ for (int i = 0; i < num_layers; i++)
    //~ {
        //~ for (int j = 0; j < N_RA_local; j++)
        //~ {
            //~ if ( (Id_RA_local[j] >= Nss * i) && (Id_RA_local[j] < Nss * (i + 1)) )
            //~ {
                //~ int diff = Nss - (Id_RA_local[j] - Nss * i);
//~ 
                //~ for (int k = 0; k < Nss; k ++)
                    //~ weights_local[j][Id_RA_local[j] + diff + k] = WEIGHT_MAX;
            //~ }
        //~ }
    //~ }
    //~ 
//~ 
    //~ this->update_all_synapses();
//~ }
//~ 
//~ void NetworkGrowthSimulator::initialize_random_chain_connections(int num_layers)
//~ {
    //~ std::vector<std::vector<int>> chain; // neurons in the chain
//~ 
    //~ chain.resize(num_layers);
    //~ chain[0].resize(N_TR);
//~ 
    //~ for (int i = 0; i < N_TR; i++)
        //~ chain[0][i] = i;
//~ 
    //~ for (int i = 1; i < num_layers; i++)
        //~ chain[i].resize(Nss);
//~ 
    //~ if (MPI_rank == 0)
    //~ {
//~ 
        //~ // connections for HVC(RA) neurons
		//~ for (int i = 0; i < N_RA; i++)
     	//~ {
	 		//~ for (int j = 0; j < N_I; j++)
	        //~ {
	         	 //~ if (generator.random(1) < p_RA2I(i,j))
	             //~ {
		        	 //~ double G = this->sample_Ge2i();
//~ 
		             //~ weights_RA_I_global[i].push_back(G);
		             //~ syn_ID_RA_I_global[i].push_back(j);
		         //~ }
			 //~ }
//~ 
		 //~ }
		//~ // connections for HVC(I) neurons
//~ 
		//~ for (int i = 0; i < N_I; i++)
     	//~ {
	 		//~ for (int j = 0; j < N_RA; j++)
	        //~ {
	         	 //~ if (generator.random(1) < p_I2RA(i,j))
	             //~ {
		        	 //~ double G = this->sample_Gi2e();
//~ 
		             //~ weights_I_RA_global[i].push_back(G);
		             //~ syn_ID_I_RA_global[i].push_back(j);
		         //~ }
			 //~ }
		 //~ }
        //~ 
        //~ // chain connections
        //~ std::vector<int> neuronsInChain; // id of HVC(RA) neurons in chain
//~ 
        //~ for (int i = 0; i < N_TR; i++)
            //~ neuronsInChain.push_back(i);
//~ 
        //~ // recruit layers
        //~ for (int i = 1; i < num_layers; i++)
        //~ {
            //~ // recruit neurons in the layer
            //~ for (int j = 0; j < Nss; j++)
            //~ {
                //~ bool neuronAlreadyInChain = false; // indicator that selected neuron is already in the chain
//~ 
                //~ int target_id; // id of the recruited target
            //~ 
                //~ do
                //~ {
                    //~ // sample target id
                    //~ target_id = generator.sample_integer(N_TR, N_RA);
                    //~ 
                    //~ // check if target is already in the chain
                    //~ std::vector<int>::iterator pos = std::find(neuronsInChain.begin(), neuronsInChain.end(), target_id);
                    //~ 
                    //~ if (pos != neuronsInChain.end())
                        //~ neuronAlreadyInChain = true;
                    //~ else
                        //~ neuronAlreadyInChain = false;
//~ 
                //~ } while (neuronAlreadyInChain);
//~ 
                //~ chain[i][j] = target_id;
                //~ neuronsInChain.push_back(target_id);
            //~ }
        //~ }
//~ 
        //~ // print chain
        //~ std::cout << "Chain:\n" << std::endl;
        //~ for (int i = 0; i < num_layers; i++)
        //~ {
            //~ for (size_t j = 0; j < chain[i].size(); j++)
                //~ std::cout << chain[i][j] << "\t";
            //~ std::cout << std::endl;
        //~ }
	//~ }
//~ 
    //~ this->send_connections();
//~ 
    //~ std::string fileRA2I = outputDirectory + "RA_I_connections.bin";
    //~ std::string fileI2RA = outputDirectory + "I_RA_connections.bin";
    //~ std::string filePajekFixed = outputDirectory + "fixed.net";
//~ 
    //~ this->write_invariable_synapses(fileRA2I.c_str(), fileI2RA.c_str());
    //~ this->write_pajek_fixed(filePajekFixed.c_str());
//~ 
    //~ // send neurons in the chain
    //~ for (int i = 1; i < num_layers; i++)
        //~ MPI_Bcast(&chain[i][0], Nss, MPI_INT, 0, MPI_COMM_WORLD);
//~ 
    //~ // connect neurons according to the chain
//~ 
    //~ for (int i = 0; i < num_layers-1; i++)
    //~ {
        //~ // all neurons in the source layer
        //~ for (size_t j = 0; j < chain[i].size(); j++)
        //~ {
            //~ // connect to all neurons in the next layer
            //~ for (size_t k = 0; k < chain[i+1].size(); k++)
            //~ {
                //~ // determine the location of the source neuron
                //~ int rank, shift;
                //~ 
                //~ this->get_neuronRA_location(chain[i][j], &rank, &shift);
                //~ 
//~ 
                //~ if (MPI_rank == rank)
                //~ {
                    //~ weights_local[shift][chain[i+1][k]] = WEIGHT_MAX;
            //~ 
                    //~ std::cout << "Rank = " << MPI_rank << " shift = " << shift << " Id_RA_local = " << Id_RA_local[shift] << " source_id = " << chain[i][j] 
                              //~ << std::endl;
                //~ }
            //~ }
        //~ }
    //~ }
//~ 
    //~ this->update_all_synapses();
    //~ this->gather_data();
    //~ 
    //~ std::string fileActiveGraph = outputDirectory + "RA_RA_active_connections.bin";
    //~ std::string fileSuperGraph = outputDirectory + "RA_RA_super_connections.bin";
	    	//~ 
    //~ this->write_supersynapses(fileSuperGraph.c_str());
    //~ this->write_active_synapses(fileActiveGraph.c_str());
//~ }
//~ 
//~ void NetworkGrowthSimulator::initialize_test_connections(int num_RA_targets, int num_RA_target_groups)
//~ {
//~ 
    //~ if (MPI_rank == 0)
    //~ {
//~ 
        //~ if (num_RA_targets > N_RA - N_TR)
            //~ std::cerr << "Number of RA targets exceeds number of pool neurons!" << std::endl;
//~ 
        //~ // connect first training neuron to a single HVC(I) neuron
        //~ weights_RA_I_global[0].push_back(Gei_mean);
        //~ syn_ID_RA_I_global[0].push_back(0);
//~ 
		//~ // connect HVC(I) neuron to 
//~ 
      //~ double inhibition_strength = 0.6;  
            //~ 
        //~ for (int j = N_TR; j < num_RA_targets; j++)
        //~ {
//~ 
             //~ if ((j - N_TR) % num_RA_target_groups == 0)
                //~ inhibition_strength += Gie_mean;
                //~ 
             //~ weights_I_RA_global[0].push_back(inhibition_strength);
             //~ syn_ID_I_RA_global[0].push_back(j);
        //~ }
	//~ }
//~ 
    //~ this->send_connections();
//~ 
    //~ std::string fileRA2I = outputDirectory + "RA_I_connections.bin";
    //~ std::string fileI2RA = outputDirectory + "I_RA_connections.bin";
    //~ std::string filePajekFixed = outputDirectory + "fixed.net";
//~ 
    //~ this->write_invariable_synapses(fileRA2I.c_str(), fileI2RA.c_str());
    //~ this->write_pajek_fixed(filePajekFixed.c_str());
//~ 
//~ }
//~ 
//~ 
//~ 

//~ void NetworkGrowthSimulator::print_invariable_connections()
//~ {
    //~ if (MPI_rank == 0)
    //~ {
        //~ for (int i = 0; i < N_RA; i++)
            //~ std::cout << "RA neuron " << i << "has " << weights_RA_I_global[i].size() << " connections to I neurons\n" << std::endl;
//~ 
        //~ for (int i = 0; i < N_I; i++)
            //~ std::cout << "I neuron " << i << "has " << weights_I_RA_global[i].size() << " connections to RA neurons" << std::endl;
    //~ }
//~ }

//~ void NetworkGrowthSimulator::print_received_invariable_connections()
//~ {
    //~ /*
    //~ for (int i = 0; i < N_RA_local; i++)
    //~ {
        //~ std::cout << "RA neuron " << Id_RA_local[i] << "has " << weights_RA_I_local[i].size() << " connections to I neurons:\n";
        //~ 
        //~ for (size_t j = 0; j < weights_RA_I_local[i].size(); j++)
            //~ std::cout << syn_ID_RA_I_local[i][j] << "\t";
//~ 
        //~ std:: cout << std::endl;
    //~ }
    //~ */
    //~ 
    //~ for (int i = 0; i < N_I_local; i++)
    //~ {
        //~ std::cout << "I neuron " << Id_I_local[i] << "has " << weights_I_RA_local[i].size() << " connections to RA neurons:\n";
        //~ 
        //~ for (size_t j = 0; j < weights_I_RA_local[i].size(); j++)
            //~ std::cout << syn_ID_I_RA_local[i][j] << "\t";
//~ 
        //~ std:: cout << std::endl;
    //~ }
//~ }
//~ 
//~ void NetworkGrowthSimulator::set_all_mature()
//~ {
    //~ for (int i = 0; i < N_RA_local; i++)
        //~ mature_local[i] = 1;
//~ 
//~ }

int NetworkGrowthSimulator::check_bad_values()
{
	for (int i = 0; i < N_RA_local; i++)
	{
		if ( HVCRA_local[i].check_bad_values() == -1 )
		{
			std::cerr << "Bad values for HVC-RA neuron " << Id_RA_local[i] << std::endl;
			return -1;
		}
	}
	
	return 0;
	
}

void NetworkGrowthSimulator::reset_after_chain_test()
{
	for (int i = 0; i < N_RA; i++)
    {
		if (MPI_rank == 0)
		{
			spikes_in_trial_soma_global[i].clear();
			spikes_in_trial_dend_global[i].clear();
        }
    }

	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_I; i++)
			spikes_in_trial_interneuron_global[i].clear();
	}
	
    for (int i = 0; i < N_RA_local; i++)
    {
		
		// update deliveries queues
		delivery_queue_RA_RA_soma[i].clear();
		delivery_queue_RA_I[i].clear();
	    HVCRA_local[i].set_to_rest();

        spikes_in_trial_soma_local[i].clear();
        spikes_in_trial_dend_local[i].clear();
    }

    for (int i = 0; i < N_I_local; i++)
    {
		// update delivery queues
		delivery_queue_I_RA[i].clear();
        HVCI_local[i].set_to_rest();
        spikes_in_trial_interneuron_local[i].clear();
    }	
}

void NetworkGrowthSimulator::setToRest_afterEpoch()
{
	for (int i = 0; i < N_RA; i++)
    {
		if (MPI_rank == 0)
		{
			spikes_in_trial_soma_global[i].clear();
			spikes_in_trial_dend_global[i].clear();
        }
    }

	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_I; i++)
			spikes_in_trial_interneuron_global[i].clear();
	}
	
    for (int i = 0; i < N_RA_local; i++)
    {
		HVCRA_local[i].set_to_rest();

        spikes_in_trial_soma_local[i].clear();
        spikes_in_trial_dend_local[i].clear();
    }

    for (int i = 0; i < N_I_local; i++)
    {
        HVCI_local[i].set_to_rest();
        spikes_in_trial_interneuron_local[i].clear();
    }	
}

void NetworkGrowthSimulator::setToRest_after_trial_soma_pre_dend_post_delays()
{
	for (int i = 0; i < N_RA; i++)
    {
		if (MPI_rank == 0)
		{
			spikes_in_trial_soma_global[i].clear();
			spikes_in_trial_dend_global[i].clear();
        }
        
		// delete bursts 
        previous_dendritic_spike_times_global[i].clear();
    }

	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_I; i++)
			spikes_in_trial_interneuron_global[i].clear();
	}
	
    for (int i = 0; i < N_RA_local; i++)
    {
		// delete delivered spikes that are too old
		for (int j = 0; j < N_RA; j++)
			delivered_spike_times[i][j].clear();
			
		delivery_queue_RA_RA_soma[i].clear();
		delivery_queue_RA_I[i].clear();
		
	    HVCRA_local[i].set_to_rest();

        spikes_in_trial_soma_local[i].clear();
        spikes_in_trial_dend_local[i].clear();
    }

    for (int i = 0; i < N_I_local; i++)
    {
		
		// clear delivery queues
		delivery_queue_I_RA[i].clear();
		
        HVCI_local[i].set_to_rest();
        spikes_in_trial_interneuron_local[i].clear();
    }	
}
	

void NetworkGrowthSimulator::reset_after_trial_soma_pre_dend_post_delays()
{
	for (int i = 0; i < N_RA; i++)
    {
		if (MPI_rank == 0)
		{
			spikes_in_trial_soma_global[i].clear();
			spikes_in_trial_dend_global[i].clear();
        }
        
       
		
		// delete bursts that are too old
        if ( !previous_dendritic_spike_times_global[i].empty() )
        {
			//~ if (MPI_rank == 0)
			//~ {
				//~ std::cout << "Somatic spikes for neuron " << i << " before erasing: \n";
				//~ std::for_each(previous_somatic_spike_times_global[i].begin(), previous_somatic_spike_times_global[i].end(),
													//~ [](double &x){std::cout << x << " ";});
				//~ std::cout << std::endl;
				//~ 
			//~ }
			
			auto it = std::lower_bound(previous_dendritic_spike_times_global[i].begin(), previous_dendritic_spike_times_global[i].end(),
											TRIAL_DURATION - STDP_WINDOW);
			
			
			previous_dendritic_spike_times_global[i].erase(previous_dendritic_spike_times_global[i].begin(), it);
			
			// update the spike times
			std::for_each(previous_dendritic_spike_times_global[i].begin(), previous_dendritic_spike_times_global[i].end(),
													[](double &x){x -= TRIAL_DURATION;});
			
			//~ if (MPI_rank == 0)
			//~ {
				//~ std::cout << "Somatic spikes for neuron " << i << " after erasing: \n";
				//~ std::for_each(previous_somatic_spike_times_global[i].begin(), previous_somatic_spike_times_global[i].end(),
													//~ [](double &x){std::cout << x << " ";});
				//~ std::cout << std::endl;
				//~ 
			//~ }
										
		}
    }

	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_I; i++)
			spikes_in_trial_interneuron_global[i].clear();
	}
	
    for (int i = 0; i < N_RA_local; i++)
    {
		
		
		// delete delivered spikes that are too old
		for (int j = 0; j < N_RA; j++)
		{
			if ( !delivered_spike_times[i][j].empty() )
			{
				//~ if (MPI_rank == 0)
				//~ {
					//~ std::cout << "Somatic spikes for neuron " << i << " before erasing: \n";
					//~ std::for_each(previous_somatic_spike_times_global[i].begin(), previous_somatic_spike_times_global[i].end(),
														//~ [](double &x){std::cout << x << " ";});
					//~ std::cout << std::endl;
					//~ 
				//~ }
				
				auto it = std::lower_bound(delivered_spike_times[i][j].begin(), delivered_spike_times[i][j].end(),
												TRIAL_DURATION - STDP_WINDOW);
				
				
				delivered_spike_times[i][j].erase(delivered_spike_times[i][j].begin(), it);
				
				// update the spike times
				std::for_each(delivered_spike_times[i][j].begin(), delivered_spike_times[i][j].end(),
														[](double &x){x -= TRIAL_DURATION;});
				
				//~ if (MPI_rank == 0)
				//~ {
					//~ std::cout << "Somatic spikes for neuron " << i << " after erasing: \n";
					//~ std::for_each(previous_somatic_spike_times_global[i].begin(), previous_somatic_spike_times_global[i].end(),
														//~ [](double &x){std::cout << x << " ";});
					//~ std::cout << std::endl;
					//~ 
				//~ }
											
			}
		}
		
		// update deliveries queues
		std::for_each(delivery_queue_RA_RA_soma[i].begin(), delivery_queue_RA_RA_soma[i].end(),
													[](std::pair<double,int> &p){p.first -= TRIAL_DURATION;});
		
		std::for_each(delivery_queue_RA_I[i].begin(), delivery_queue_RA_I[i].end(),
													[](std::pair<double,int> &p){p.first -= TRIAL_DURATION;});
		
	    HVCRA_local[i].reset_time();

        spikes_in_trial_soma_local[i].clear();
        spikes_in_trial_dend_local[i].clear();
    }

    for (int i = 0; i < N_I_local; i++)
    {
		// update delivery queues
		std::for_each(delivery_queue_I_RA[i].begin(), delivery_queue_I_RA[i].end(),
													[](std::pair<double,int> &p){p.first -= TRIAL_DURATION;});
		
        HVCI_local[i].reset_time();
        spikes_in_trial_interneuron_local[i].clear();
    }	
}
	
void NetworkGrowthSimulator::reset_after_trial_soma_pre_dend_post_no_delays()
{
	for (int i = 0; i < N_RA; i++)
    {
		if (MPI_rank == 0)
		{
			spikes_in_trial_soma_global[i].clear();
			spikes_in_trial_dend_global[i].clear();
        }
        
        //~ // delete spikes that are too old
        //~ if ( !previous_somatic_spike_times_global[i].empty() )
        //~ {
			//~ /*
			//~ if (MPI_rank == 0)
			//~ {
				//~ std::cout << "Somatic spikes for neuron " << i << " before erasing: \n";
				//~ std::for_each(previous_somatic_spike_times_global[i].begin(), previous_somatic_spike_times_global[i].end(),
													//~ [](double &x){std::cout << x << " ";});
				//~ std::cout << std::endl;
				//~ 
			//~ }*/
			//~ 
			//~ auto it = std::lower_bound(previous_somatic_spike_times_global[i].begin(), previous_somatic_spike_times_global[i].end(),
											//~ TRIAL_DURATION - STDP_WINDOW);
			//~ 
			//~ 
			//~ previous_somatic_spike_times_global[i].erase(previous_somatic_spike_times_global[i].begin(), it);
			//~ 
			//~ // update the spike times
			//~ std::for_each(previous_somatic_spike_times_global[i].begin(), previous_somatic_spike_times_global[i].end(),
													//~ [](double &x){x -= TRIAL_DURATION;});
			//~ 
			//~ /*if (MPI_rank == 0)
			//~ {
				//~ std::cout << "Somatic spikes for neuron " << i << " after erasing: \n";
				//~ std::for_each(previous_somatic_spike_times_global[i].begin(), previous_somatic_spike_times_global[i].end(),
													//~ [](double &x){std::cout << x << " ";});
				//~ std::cout << std::endl;
				//~ 
			//~ }*/
										//~ 
		//~ }
		
		
		
		
		// delete bursts that are too old
        if ( !previous_dendritic_spike_times_global[i].empty() )
        {
			//~ if (MPI_rank == 0)
			//~ {
				//~ std::cout << "Somatic spikes for neuron " << i << " before erasing: \n";
				//~ std::for_each(previous_somatic_spike_times_global[i].begin(), previous_somatic_spike_times_global[i].end(),
													//~ [](double &x){std::cout << x << " ";});
				//~ std::cout << std::endl;
				//~ 
			//~ }
			
			auto it = std::lower_bound(previous_dendritic_spike_times_global[i].begin(), previous_dendritic_spike_times_global[i].end(),
											TRIAL_DURATION - STDP_WINDOW);
			
			
			previous_dendritic_spike_times_global[i].erase(previous_dendritic_spike_times_global[i].begin(), it);
			
			// update the spike times
			std::for_each(previous_dendritic_spike_times_global[i].begin(), previous_dendritic_spike_times_global[i].end(),
													[](double &x){x -= TRIAL_DURATION;});
			
			//~ if (MPI_rank == 0)
			//~ {
				//~ std::cout << "Somatic spikes for neuron " << i << " after erasing: \n";
				//~ std::for_each(previous_somatic_spike_times_global[i].begin(), previous_somatic_spike_times_global[i].end(),
													//~ [](double &x){std::cout << x << " ";});
				//~ std::cout << std::endl;
				//~ 
			//~ }
										
		}
    }

	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_I; i++)
			spikes_in_trial_interneuron_global[i].clear();
	}
	
    for (int i = 0; i < N_RA_local; i++)
    {
		
		 // delete spikes that are too old
        if ( !previous_somatic_spike_times_local[i].empty() )
        {
			//~ if (MPI_rank == 0)
			//~ {
				//~ std::cout << "Somatic spikes for neuron " << i << " before erasing: \n";
				//~ std::for_each(previous_somatic_spike_times_global[i].begin(), previous_somatic_spike_times_global[i].end(),
													//~ [](double &x){std::cout << x << " ";});
				//~ std::cout << std::endl;
				//~ 
			//~ }
			
			auto it = std::lower_bound(previous_somatic_spike_times_local[i].begin(), previous_somatic_spike_times_local[i].end(),
											TRIAL_DURATION - STDP_WINDOW);
			
			
			previous_somatic_spike_times_local[i].erase(previous_somatic_spike_times_local[i].begin(), it);
			
			// update the spike times
			std::for_each(previous_somatic_spike_times_local[i].begin(), previous_somatic_spike_times_local[i].end(),
													[](double &x){x -= TRIAL_DURATION;});
			
			//~ if (MPI_rank == 0)
			//~ {
				//~ std::cout << "Somatic spikes for neuron " << i << " after erasing: \n";
				//~ std::for_each(previous_somatic_spike_times_global[i].begin(), previous_somatic_spike_times_global[i].end(),
													//~ [](double &x){std::cout << x << " ";});
				//~ std::cout << std::endl;
				//~ 
			//~ }
										
		}
		
	    HVCRA_local[i].reset_time();

        spikes_in_trial_soma_local[i].clear();
        spikes_in_trial_dend_local[i].clear();
    }

    for (int i = 0; i < N_I_local; i++)
    {
        HVCI_local[i].reset_time();
        spikes_in_trial_interneuron_local[i].clear();
    }
	
	
}

/*void NetworkGrowthSimulator::reset_after_trial()
{
	for (int i = 0; i < N_RA; i++)
    {
		if (MPI_rank == 0)
		{
			spikes_in_trial_soma_global[i].clear();
			spikes_in_trial_dend_global[i].clear();
        }
        
        // delete spikes that are too old
        if ( !previous_somatic_spike_times_global[i].empty() )
        {
			//~ if (MPI_rank == 0)
			//~ {
				//~ std::cout << "Somatic spikes for neuron " << i << " before erasing: \n";
				//~ std::for_each(previous_somatic_spike_times_global[i].begin(), previous_somatic_spike_times_global[i].end(),
													//~ [](double &x){std::cout << x << " ";});
				//~ std::cout << std::endl;
				//~ 
			//~ }
			
			auto it = std::lower_bound(previous_somatic_spike_times_global[i].begin(), previous_somatic_spike_times_global[i].end(),
											TRIAL_DURATION - STDP_WINDOW);
			
			
			previous_somatic_spike_times_global[i].erase(previous_somatic_spike_times_global[i].begin(), it);
			
			// update the spike times
			std::for_each(previous_somatic_spike_times_global[i].begin(), previous_somatic_spike_times_global[i].end(),
													[](double &x){x -= TRIAL_DURATION;});
			
			//~ if (MPI_rank == 0)
			//~ {
				//~ std::cout << "Somatic spikes for neuron " << i << " after erasing: \n";
				//~ std::for_each(previous_somatic_spike_times_global[i].begin(), previous_somatic_spike_times_global[i].end(),
													//~ [](double &x){std::cout << x << " ";});
				//~ std::cout << std::endl;
				//~ 
			//~ }
										
		}
    }

	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_I; i++)
			spikes_in_trial_interneuron_global[i].clear();
	}
	
    for (int i = 0; i < N_RA_local; i++)
    {
	    HVCRA_local[i].reset_time();

        spikes_in_trial_soma_local[i].clear();
        spikes_in_trial_dend_local[i].clear();
    }

    for (int i = 0; i < N_I_local; i++)
    {
        HVCI_local[i].reset_time();
        spikes_in_trial_interneuron_local[i].clear();
    }
	
	
}*/

//~ 
//~ void NetworkGrowthSimulator::randomize_after_trial()
//~ {
    //~ for (int i = 0; i < N_RA; i++)
    //~ {
        //~ spikes_in_trial_soma_global[i].clear();
        //~ spikes_in_trial_dend_global[i].clear();
//~ 
    //~ }
//~ 
    //~ for (int i = 0; i < N_I; i++)
        //~ spikes_in_trial_interneuron_global[i].clear();
//~ 
    //~ for (int i = 0; i < N_RA_local; i++)
    //~ {
	    //~ HVCRA_local[i].set_to_rest();
//~ 
        //~ spikes_in_trial_soma_local[i].clear();
        //~ spikes_in_trial_dend_local[i].clear();
    //~ }
//~ 
    //~ for (int i = 0; i < N_I_local; i++)
    //~ {
        //~ HVCI_local[i].set_to_rest();
        //~ spikes_in_trial_interneuron_local[i].clear();
    //~ }
//~ }
//~ 
//~ 

void NetworkGrowthSimulator::set_training_current(double t)
{
    std::function<double (double)> f = std::bind(&training_current, t, _1);

	for (int i = 0; i < N_TR; i++)
	{
		int rank;
		int shift;
		
		this->get_neuronRA_location(training_neurons[i], &rank, &shift);
		
		if (MPI_rank == rank)
			HVCRA_local[shift].set_dend_current(f);
	}
}

 //~ 
//~ 
//~ void NetworkGrowthSimulator::set_training_current()
//~ {
	//~ double current_injection_time; // time of current injection to training neurons
	//~ 
    //~ if (MPI_rank == 0)
        //~ current_injection_time = trial_number * trial_duration + WAITING_TIME + generator.random(trial_duration-2*WAITING_TIME);
        //~ 
    //~ // send current injection time to all processes
    //~ MPI_Bcast(&current_injection_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//~ 
	//~ std::function<double (double)> f = std::bind(&training_current, current_injection_time, _1);
//~ 
	//~ for (int i = 0; i < N_TR; i++)
	//~ {
		//~ int rank;
		//~ int shift;
		//~ 
		//~ this->get_neuronRA_location(training_neurons[i], &rank, &shift);
		//~ 
		//~ if (MPI_rank == rank)
			//~ HVCRA_local[shift].set_dend_current(f);
	//~ }
//~ }

//~ void NetworkGrowthSimulator::test_ideal_chain(int num_layers, int num_trials)
//~ {
    //~ // initialize coordintates and ideal chain connections
    //~ this->initialize_coordinates();
    //~ this->write_all_coordinates();
    //~ this->initialize_ideal_chain_connections(num_layers);
   //~ 
    //~ this->gather_data();
//~ 
    //~ //this->print_received_invariable_connections();
//~ 
    //~ 
//~ 
    //~ // write active and super synapses to files
    //~ std::string fileActiveGraph = outputDirectory + "RA_RA_active_connections.bin";
    //~ std::string fileSuperGraph = outputDirectory + "RA_RA_super_connections.bin";
	    	//~ 
    //~ this->write_supersynapses(fileSuperGraph.c_str());
    //~ this->write_active_synapses(fileActiveGraph.c_str());
//~ 
    //~ // make all neurons mature
    //~ for (int i = 0; i < N_RA_local; i++)
        //~ gaba_potential_local[i] = E_GABA_MATURE;
//~ 
    //~ // unused array with average dendritic spike times. Needed for gathering mature data
    //~ std::vector<std::vector<double>> average_dendritic_spike_time; // array with average dendritic spike time in every trial
//~ 
	//~ average_dendritic_spike_time.resize(N_RA);
    //~ 
    //~ for (int i = 0; i < num_trials; i++)
    //~ {
        //~ this->mature_trial();
        //~ this->gather_mature_data(average_dendritic_spike_time);
//~ 
        //~ for (int n = 1; n <= num_layers; n++)
            //~ for (int j = 0; j < Nss; j++)
                //~ this->write_RA((outputDirectory + "RA/RA" + std::to_string(n*Nss + j) + "_trial" + std::to_string(i+1) + ".bin").c_str(), n*Nss + j);
        //~ 
        //~ for (int n = 0; n <= num_layers; n++)
            //~ for (int j = 0; j < Nss; j++)
                //~ this->write_I((outputDirectory + "I/I" + std::to_string(n*Nss + j) + "_trial" + std::to_string(i+1) + ".bin").c_str(), n*Nss + j);
//~ 
        //~ // write all occured spikes to files
        //~ this->write_soma_spike_times((outputDirectory + "soma_spikes_trial" + std::to_string(i+1) + ".bin").c_str());
        //~ this->write_dend_spike_times((outputDirectory + "dend_spikes_trial" + std::to_string(i+1) + ".bin").c_str());
        //~ 
        //~ this->randomize_after_trial();
    //~ }   
//~ }
//~ 
//~ void NetworkGrowthSimulator::continue_growth(std::string dataDir, int starting_trial, int save_freq_short, int save_freq_long)
//~ {
	//~ this->read_network_state(dataDir, starting_trial); // read data from file
    //~ 
    //~ outputDirectory = dataDir;
    //~ 
    //~ std::cout << "Trial number : " << trial_number << std::endl;
    //~ 
    //~ bool training = true;
    //~ 
    //~ this->chain_growth(training, save_freq_short, save_freq_long);
    //~ 
//~ }

//~ void NetworkGrowthSimulator::test_grown_chain(int num_trials, std::string dataDir, int starting_trial, std::string outputDir)
//~ {
    //~ this->read_network_state(dataDir, starting_trial); // read data from file
    //~ trial_number = 0;
    //~ outputDirectory = outputDir;
    //~ this->test_mature_chain(num_trials); // test networkc
//~ }

//~ void NetworkGrowthSimulator::test_random_chain(int num_layers, int num_trials)
//~ {
    //~ this->initialize_coordinates();
    //~ this->write_all_coordinates();
    //~ this->initialize_random_chain_connections(num_layers);
//~ 
    //~ // set all neurons to be mature
    //~ for (int i = 0; i < N_RA_local; i++)
    //~ {
        //~ mature_local[i] = 1;
        //~ gaba_potential_local[i] = E_GABA_MATURE;
    //~ }
//~ 
    //~ this->test_mature_chain(num_trials);
//~ }

//~ void NetworkGrowthSimulator::test_mature_chain(int num_trials)
//~ {
    //~ std::string file_soma_spikes = outputDirectory + "soma_spikes_in_trial.bin"; // file with somatic spikes in trial
    //~ std::string file_dend_spikes = outputDirectory + "dend_spikes_in_trial.bin"; // file with dendritic spikes in trial
    //~ std::string file_chain_test = outputDirectory + "mature_chain_test.bin"; // trial with mature chain test info
//~ 
	//~ std::vector<std::vector<double>> average_dendritic_spike_time; // array with average dendritic spike time in every trial
	//~ std::vector<std::vector<int>> num_dendritic_spikes_in_trials; // number of dendritic spikes produced in all trials
	//~ std::vector<std::vector<int>> num_somatic_spikes_in_trials; // number of somatic spikes produced in all trials
//~ 
	//~ average_dendritic_spike_time.resize(N_RA);
    //~ num_dendritic_spikes_in_trials.resize(N_RA);
    //~ num_somatic_spikes_in_trials.resize(N_RA);
//~ 
    //~ if (MPI_rank == 0)
    //~ {
        //~ for (int j = 0; j < N_RA; j++)
        //~ {
            //~ num_dendritic_spikes_in_trials[j].resize(num_trials);
            //~ num_somatic_spikes_in_trials[j].resize(num_trials);
        //~ }
    //~ }
    //~ // neurons for gabaMaturation 300117	
    //~ /*std::vector<int> RAtoWrite = {71, 186, 187, 84, 44, 175, 219, 238, 224, 70, 288, 117, 99, 276, 23, 24, 165, 
                                  //~ 128, 184, 155, 114, 203, 257, 65, 273, 183, 294, 19, 35, 97, 142, 233, 6, 192, 
                                  //~ 248, 295, 38, 69, 207, 268, 49, 263, 132, 101, 33, 206, 90, 252, 77, 43, 293, 36, 
                                  //~ 5, 180, 282, 65, 34, 267, 208, 66, 146, 179};
    //~ */
    //~ /*
    //~ // neurons for gabaMaturation 310117	
    //~ std::vector<int> RAtoWrite = {201, 209, 124, 275, 40, 87, 66, 282, 222, 285, 115, 58, 183, 123, 244, 96, 226,
                                  //~ 110, 15, 20, 178, 215, 192, 128, 280, 38, 7, 235, 273, 258, 227, 132, 169, 172, 
                                  //~ 243, 100, 188};
    //~ */
    //~ // neurons for gabaMaturation 010217	
    //~ /*std::vector<int> RAtoWrite = {179, 129, 128, 268, 130, 142, 15, 115, 273, 19, 23, 282, 29, 261, 290, 163, 292, 
                                  //~ 37, 167, 168, 169, 199, 172, 51, 182, 60, 68, 69, 256, 201, 207, 208, 209, 82, 85, 
                                  //~ 87, 90, 92, 122, 144, 226, 227, 131, 101, 81, 259, 231, 110, 114, 243, 117, 120, 250, 123, 124, 213};
    //~ */
	//~ /*
    //~ // neurons for gabaMaturation270317 huxley
    //~ std::vector<int> RAtoWrite = {179, 66, 11, 123, 173, 129, 148, 287, 199, 174, 285, 298, 144, 20, 161, 165, 205, 89, 17}; 
    //~ */
    //~ // neurons for gabaMaturation270317 huxley
    //~ //std::vector<int> RAtoWrite = {179, 66, 11, 123, 173, 129, 148, 287, 199, 174, 285, 298, 144, 20, 161, 165, 205, 89, 17}; 
    //~ //std::vector<int> ItoWrite;
//~ 
	//~ /*
    //~ // neurons for gabaMaturation010417 hodgkin
    //~ std::vector<int> RAtoWrite = {281, 156, 84, 52, 16, 92, 238, 75, 47, 10, 283, 171, 115, 194, 225, 78, 268, 221, 289, 104,
                                  //~ 185, 285, 287, 21, 58, 55, 229, 222, 145, 239, 123, 173, 295, 179, 240, 134, 280, 42, 228, 178, 
                                  //~ 208, 244, 294, 130, 45, 4, 217, 143, 87, 226, 148, 233, 190, 223, 255, 138, 29, 192, 290, 12, 
                                  //~ 142, 129, 150, 48, 69, 271, 174, 17, 167, 168, 273, 68, 35, 95, 163, 207, 128, 172, 231, 258, 
                                  //~ 99, 30, 100}; 
    //~ */
    //~ /*
    //~ // neurons for gabaMaturation280317 huxley
    //~ std::vector<int> RAtoWrite = {111, 253, 62, 265, 260, 8, 291, 160, 143, 64, 271, 128, 134, 84, 38, 72, 267, 34, 137, 77, 
                                  //~ 20, 188, 200, 136, 173, 13, 206, 5, 118};
    //~ */
    //~ /*
    //~ // neurons for gabaMaturation040417 huxley
    //~ std::vector<int> RAtoWrite = {85, 197, 201, 44, 262, 247, 228, 249, 185, 46, 199, 212, 64, 140, 174, 210, 236, 77, 129, 
								  //~ 15, 39, 298, 168, 216, 142, 295, 204, 13, 23, 34, 280, 186, 299, 121, 54, 269, 292, 105, 9, 
								  //~ 35, 57, 251, 100, 69, 260, 182, 136, 237, 134, 26, 66, 157, 286, 135, 193, 45, 219, 80, 20, 
								  //~ 126, 196, 211, 6, 190, 257, 81, 104, 36, 253, 25, 90, 115, 30, 183, 63, 109, 266, 202, 94, 113, 
								  //~ 222, 187, 246, 86, 206, 232, 160, 125, 240, 117, 282, 152, 19, 259, 198, 128};
    //~ */
    //~ /*
    //~ // neurons for gabaMaturation130417 huxley
    //~ std::vector<int> RAtoWrite = {51, 48, 146, 172, 132, 277, 203, 175, 275, 28, 31, 37, 140, 235, 67, 245, 21, 50, 138, 93, 76,
									//~ 228, 46, 225, 187, 231, 156, 210, 246, 148, 7, 49, 195, 74, 124, 255, 169, 152, 269, 206, 260, 
									//~ 94, 83, 259, 57, 171, 114, 23, 222, 248, 113, 165, 20, 104, 116, 59, 257, 25, 26, 89, 252, 151, 
									//~ 229, 253, 106, 176, 115, 183, 283, 30, 112, 226, 267, 139, 238, 158, 167, 95, 84, 268, 162, 111, 164, 163};
	//~ */
	//~ 
	//~ // neurons for gabaMaturation280317 hodgkin
    //~ //std::vector<int> RAtoWrite = {102, 18, 13, 141, 88, 269, 40, 286, 256, 63, 119, 262, 41, 109, 195, 128, 276, 153, 288, 271,
	//~ //							  51, 204, 231, 120, 173, 215, 124, 280, 273, 163, 259, 127, 146, 261, 21, 66, 65, 267, 98, 34, 22, 221, 
	//~ //							  44, 290, 125, 235, 249, 46, 90, 138, 137}; 
	//~ 
    //~ // neurons for gabaMaturation300317 hodgkin
    //~ //std::vector<int> RAtoWrite = {297, 28, 278, 286, 225, 194, 292, 78, 15, 284, 14, 299, 240, 122, 59, 228, 145, 80, 239, 254,
		//~ //						  79, 98, 298, 65, 197, 248, 91, 211, 133, 215, 192, 223, 128, 136, 68, 200, 234, 120, 188, 100, 
			//~ //					  172, 89}; 
	//~ 
	//~ /*
	//~ // neurons for gabaMaturation090417 hodgkin
    //~ std::vector<int> RAtoWrite = {187, 209, 32, 155, 136, 271, 74, 267, 135, 181, 260, 122, 117, 279, 227, 11, 177, 169, 166, 148, 
		//~ 251, 126, 116, 139, 86, 275, 151, 132, 248, 61, 119, 52, 90, 120, 174, 16, 214, 33, 273, 101, 189, 76, 173, 96, 261, 21, 112,
		 //~ 150, 178, 179, 84, 17, 68, 192, 31, 92, 224, 254, 56, 257, 144, 171, 288, 223, 199, 159, 237, 118, 121, 168, 228, 263, 161, 234, 
		 //~ 128, 71, 75, 200, 218, 99, 167, 217, 95, 48};
	//~ */
	//~ 
     //~ // neurons for gabaMaturation170417 huxley
     //~ /*
    //~ std::vector<int> RAtoWrite = {62, 215, 293, 184, 270, 41, 61, 146, 277, 122, 132, 288, 133, 254, 275, 249, 245, 140, 67, 33, 125,
								   //~ 21, 32, 77, 228, 192 , 263, 46, 194, 210, 287, 230, 131, 145, 99, 71, 56, 64, 266, 204, 147, 12, 163, 
								   //~ 280, 252, 162, 59, 268, 104, 190, 25, 183, 253, 114, 165, 181, 109, 83, 286, 4, 240, 128, 152, 241,
								    //~ 269, 170, 23, 271, 103, 19, 34, 276, 151, 191, 68, 87, 13, 11, 142, 278, 166, 70, 294, 5, 189, 35, 58,
								     //~ 299, 217, 42, 130, 9};
    //~ 
    //~ // 090617_lionx_2
    //~ std::vector<int> RAtoWrite = {150, 21, 286, 60, 68, 142, 287, 125, 117, 115, 191, 94, 47, 162, 137, 250, 90, 32, 164, 267, 65, 69, 95,
								  //~ 135, 241, 228, 124, 55, 13, 237, 88, 171, 173, 214, 119, 160, 143, 81, 53, 148, 260, 175, 299,
								   //~ 71, 227, 96, 283, 273, 93, 235, 157, 98, 291, 39, 244, 200, 72, 242, 106, 245, 85, 158, 40, 
								   //~ 104, 139, 102, 251, 15, 243, 207, 266, 166, 209, 276, 66, 159, 263, 31, 133, 87, 295, 231, 177,
								    //~ 10, 138, 79, 24, 136, 284, 195, 223, 34, 75, 28, 49, 247, 259, 107, 256, 219, 199, 203, 33, 
								    //~ 37, 169, 280, 184, 25, 109, 196, 212, 226, 202, 132, 91, 11, 113, 180, 272, 59, 151, 126, 26, 
								    //~ 298, 74, 269, 257, 217, 57, 185, 86, 220, 154, 80};
	//~ // 090617_lionx_3
    //~ std::vector<int> RAtoWrite = {130, 274, 286, 162, 32, 164, 47, 205, 241, 88, 99, 221, 51, 65, 198, 58, 119, 214, 265, 5, 126, 292, 67, 
								  //~ 98, 288, 129, 18, 43, 96, 291, 22, 189, 180, 131, 108, 25, 280, 290, 271, 66, 169, 182, 163, 177, 133, 49, 
	//~							      202, 31, 207, 224, 136, 231, 41, 79, 270, 223, 24, 192, 89, 284, 30, 259, 210, 239, 195, 63, 61, 247, 122};							  
    //~ 
    //~ // 160617_lionx_2
    //~ std::vector<int> RAtoWrite = {21, 150, 287, 179, 130, 155, 117, 250, 7, 137, 90, 135, 47, 164, 241, 162, 69, 198, 81, 99, 65, 218, 171, 173, 
								  //~ 262, 246, 214, 83, 201, 53, 118, 228, 227, 204, 134, 288, 283, 285, 76, 149, 71, 106, 39, 233, 213, 282, 120, 
								  //~ 178, 14, 30, 176, 296, 114, 122, 192, 184, 62, 61, 289, 253, 248, 224, 243, 37, 33, 209, 169, 109, 182, 290, 207, 
								  //~ 258, 196, 280, 104, 251, 28, 87, 226, 66, 133, 75, 31, 295, 136, 231, 70, 263, 177, 138, 195, 79, 223, 24, 284, 247, 259};
    //~ 
    
    //~ // 160617_lionx_3
    //~ std::vector<int> RAtoWrite = {117, 274, 130, 167, 254, 164, 137, 241, 51, 23, 161, 221, 65, 67, 1, 55, 292, 282, 168, 232, 44, 154, 156, 14, 204, 269, 185, 
								  //~ 184, 57, 106, 225, 26, 280, 169, 131, 182, 271};
    
    //~ // 120617_lionx_1
    //~ std::vector<int> RAtoWrite = {174, 113, 150, 26, 277, 9, 269, 151, 117, 39, 167, 155, 217, 274, 240, 162, 241, 252, 58, 57, 180, 137, 33, 69, 175, 228, 183, 124, 129, 213, 214, 131, 
								  //~ 143, 56, 23, 166, 118, 227, 244, 108, 258, 266, 36, 133, 152, 83, 149, 87, 209, 89, 177, 54, 295, 176, 253, 41, 122, 239};
    //~ 
    //~ // 160617_lionx_5
    //~ std::vector<int> RAtoWrite = {246, 136, 59, 150, 284, 185, 195, 121, 202, 118, 68, 75, 117, 26, 217, 247, 36, 167, 33, 180, 279, 131, 169, 80, 14, 265, 241, 
								//~ 111, 5, 18, 292, 205, 152, 67, 35, 144, 207};
    //~ 
    //~ 
    //~ // 120617_huxley
    //~ std::vector<int> RAtoWrite = {277, 284, 285, 74, 274, 226, 259, 79, 220, 256, 178, 269, 162, 287, 263, 196, 63, 137, 169, 210, 154, 290, 114, 61, 292, 248, 89, 18, 41, 67, 224, 279};
    //~ */
    //~ 
    //~ // 120617_lionx_2
    //~ //std::vector<int> RAtoWrite = {195, 274, 259, 79, 298, 130, 90, 118, 256, 162, 177, 247, 254, 12, 176, 67, 295, 292, 122, 144};
    //~ 
    //~ std::vector<int> RAtoWrite = {};
    //~ 
    //~ std::vector<int> ItoWrite;
//~ 
    //~ //for (int i = 0; i < N_I; i++)
    //~ //    ItoWrite.push_back(i);
//~ 
    //~ for (int i = 0; i < num_trials; i++)
	//~ {
        //~ if (MPI_rank == 0)
            //~ std::cout << "Trial " << i << std::endl;
//~ 
		//~ this->mature_trial();
		//~ this->gather_data();
		//~ //this->gather_mature_data(average_dendritic_spike_time);
		//~ 
		//~ if (MPI_rank == 0)
		//~ {
			//~ for (int i = 0; i < N_RA; i++)
			//~ {
				//~ if (spikes_in_trial_dend_global[i].size() > 0)
				//~ {
					//~ double average_spike_time = std::accumulate(spikes_in_trial_dend_global[i].begin(), spikes_in_trial_dend_global[i].end(), 0.0) / static_cast<double>(spikes_in_trial_dend_global[i].size());
//~ 
					//~ //printf("Average dendritic spike time = %f\n", average_spike_time);
					//~ average_dendritic_spike_time[i].push_back(average_spike_time);
				//~ }
			//~ }
		//~ }
//~ 
        //~ if (MPI_rank == 0)
        //~ {
            //~ for (int j = 0; j < N_RA; j++)
            //~ {
                //~ num_dendritic_spikes_in_trials[j][i] = static_cast<int>(spikes_in_trial_dend_global[j].size());
                //~ num_somatic_spikes_in_trials[j][i] = static_cast<int>(spikes_in_trial_soma_global[j].size());
            //~ }
//~ 
        //~ }
//~ 
        //~ for (size_t j = 0; j < RAtoWrite.size(); j++)
            //~ this->write_RA((outputDirectory + "RA/RA" + std::to_string(RAtoWrite[j]) + "_trial" + std::to_string(i+1) + ".bin").c_str(), RAtoWrite[j]);
        //~ 
        //~ for (size_t j = 0; j < ItoWrite.size(); j++)
            //~ this->write_I((outputDirectory + "I/I" + std::to_string(ItoWrite[j]) + "_trial" + std::to_string(i+1) + ".bin").c_str(), ItoWrite[j]);
        //~ 
       //~ 
		//~ this->write_soma_spike_times((outputDirectory + "soma_spikes_in_trial_" + std::to_string(trial_number) + "_.bin").c_str());
		//~ this->write_dend_spike_times((outputDirectory + "dend_spikes_in_trial_" + std::to_string(trial_number) + "_.bin").c_str());
        //~ 
		//~ this->randomize_after_trial();
		//~ this->set_time_for_neurons(0.0);
		//~ trial_number++;
	//~ }
	//~ 
	//~ // process dendritic spikes
//~ 
	//~ std::vector<double> mean_burst_time; // average of dendritic spike time
	//~ std::vector<double> std_burst_time; // standard deviation of dendritic spike time
    //~ std::vector<double> average_num_dend_spikes_in_trial; // average number of dendritic spikes in trial
    //~ std::vector<double> average_num_soma_spikes_in_trial; // average number of somatic spikes in trials
//~ 
    //~ std::vector<int> num_trials_with_dend_spikes; // number of trials in which neuron produced dendritic spikes
//~ 
	//~ mean_burst_time.resize(N_RA);
	//~ std_burst_time.resize(N_RA);
    //~ average_num_dend_spikes_in_trial.resize(N_RA);
    //~ average_num_soma_spikes_in_trial.resize(N_RA);
    //~ num_trials_with_dend_spikes.resize(N_RA);
//~ 
	//~ if (MPI_rank == 0)
	//~ {
		//~ for (int i = 0; i < N_RA; i++)
		//~ {
            //~ average_num_dend_spikes_in_trial[i] = std::accumulate(num_dendritic_spikes_in_trials[i].begin(), num_dendritic_spikes_in_trials[i].end(), 0.0) 
                                                //~ / static_cast<double>(num_trials);
//~ 
            //~ average_num_soma_spikes_in_trial[i] = std::accumulate(num_somatic_spikes_in_trials[i].begin(), num_somatic_spikes_in_trials[i].end(), 0.0) 
                                                //~ / static_cast<double>(num_trials);
            //~ 
            //~ num_trials_with_dend_spikes[i] = static_cast<int>(average_dendritic_spike_time[i].size());
			//~ 
            //~ /*
            //~ if (i == 288)
            //~ {
                //~ std::cout << "Number of somatic spikes in all trials: " << std::endl;
                //~ for (int j = 0; j < num_trials; j++)
                    //~ std::cout << num_somatic_spikes_in_trials[i][j] << '\t';
//~ 
                //~ std::cout << '\n' << std::endl;
                //~ 
                //~ std::cout << "Number of dendritic spikes in all trials: " << std::endl;
                //~ for (int j = 0; j < num_trials; j++)
                    //~ std::cout << num_dendritic_spikes_in_trials[i][j] << "\t";
//~ 
                //~ std::cout << std::endl << std::endl;
                //~ 
                //~ std::cout << "Number of trials in which neuron produced dendritic bursts: " << num_trials_with_dend_spikes[i] << "\n" << std::endl;
                //~ std::cout << "All dendritic burst times in these trials: " << std::endl;
                //~ 
                //~ for (size_t j = 0; j < average_dendritic_spike_time[i].size(); j++)
                    //~ std::cout << average_dendritic_spike_time[i][j] << "\t";
//~ 
                //~ std::cout << std::endl;
//~ 
            //~ }
            //~ */
//~ 
            //~ if (num_trials_with_dend_spikes[i] > 0)
            //~ {
//~ 
                //~ //for (int j = 0; j < (int) average_dendritic_spike_time[i].size(); j++)
                  //~ //  printf("average_dendritic_spike_time[%d][%d] = %f\n", i, j, average_dendritic_spike_time[i][j]);
				//~ mean_burst_time[i] = std::accumulate(average_dendritic_spike_time[i].begin(), average_dendritic_spike_time[i].end(), 0.0) / 
                                        //~ static_cast<double>(num_trials_with_dend_spikes[i]);
//~ 
            //~ }
//~ 
			//~ else
				//~ mean_burst_time[i] = -1;
			//~ // calculate standard deviation of burst times
//~ 
			//~ double accum = 0;
			//~ double mean = mean_burst_time[i];
//~ 
			//~ std::for_each(average_dendritic_spike_time[i].begin(), average_dendritic_spike_time[i].end(), [&accum, mean](const double t)
			//~ {
				//~ accum += (t - mean) * (t - mean);
			//~ });
//~ 
			//~ if (static_cast<int>(average_dendritic_spike_time[i].size() > 1))
				//~ std_burst_time[i] = sqrt(accum / (static_cast<double>(average_dendritic_spike_time[i].size()) - 1));
			//~ else
				//~ std_burst_time[i] = -1;
//~ 
//~ 
		//~ }
	//~ }
//~ 
	//~ this->write_chain_test(num_trials, num_trials_with_dend_spikes, average_num_dend_spikes_in_trial, average_num_soma_spikes_in_trial, 
                           //~ mean_burst_time, std_burst_time, file_chain_test.c_str());
//~ }
//~ 
//~ void NetworkGrowthSimulator::mature_trial()
//~ {
	//~ int some_RA_neuron_fired_soma_local;
//~ 
	//~ int some_I_neuron_fired_local;
	//~ int some_RA_neuron_fired_soma_global;
//~ 
	//~ int some_I_neuron_fired_global;
    //~ 
    //~ internal_time = 0;
    //~ network_time = internal_time + network_update_frequency;
//~ 
	//~ // set training current
	//~ double current_injection_time = 100; //ms
    //~ this->set_training_current(current_injection_time);
//~ 
	//~ // initialize update arrays and fired indicators
	//~ std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	//~ std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	//~ std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	//~ std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	//~ std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	//~ std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
//~ 
    //~ some_RA_neuron_fired_soma_local = 0;
	//~ some_RA_neuron_fired_soma_global = 0;
//~ 
    //~ some_I_neuron_fired_local = 0;
    //~ some_I_neuron_fired_global = 0;
	//~ 
    //~ // evolve dynamics
    //~ for (int t = 1; t < size; t++)
	//~ {
		//~ internal_time += timeStep;
		//~ 
		//~ for (int i = 0; i < N_RA_local; i++)
		//~ {
            //~ // set GABA potential
            //~ HVCRA_local[i].set_Ei(gaba_potential_local[i]);
            //~ 
            //~ // Debraband step
            //~ HVCRA_local[i].Debraband_step_no_target_update();
            //~ 
            //~ // if some neuron produced somatic spike
            //~ if (HVCRA_local[i].get_fired_soma())
            //~ {
                //~ spikes_in_trial_soma_local[i].push_back(internal_time);
                //~ 
                //~ // update conductances of targets
                //~ 
				//~ some_RA_neuron_fired_soma_local = 1;
				//~ // loop over all inhibitory targets of fired neurons
				//~ size_t num_I_targets = syn_ID_RA_I_local[i].size();
				//~ for (size_t j = 0; j < num_I_targets; j++)
				//~ {
					//~ int syn_ID = syn_ID_RA_I_local[i][j];
					//~ update_Ge_I_local[syn_ID] += weights_RA_I_local[i][j];
//~ 
				//~ }
				//~ 
				//~ // loop over all excitatory targets
				//~ size_t num_RA_targets = active_synapses_local[i].size();
				//~ 
				//~ //std::cout << "Neuron fired: " << Id_RA_local[i] << " num_RA_targets: " << num_RA_targets << std::endl;
//~ 
				//~ for (size_t j = 0; j < num_RA_targets; j++)
				//~ {
					//~ int syn_ID = active_synapses_local[i][j];
					//~ update_Ge_RA_local[syn_ID] += weights_local[i][syn_ID];
					//~ //std::cout << "Neuron fired: " << Id_RA_local[i] << " target: " << syn_ID << " synaptic weight: " << weights_local[i][syn_ID] << std::endl;
				//~ }
                //~ 
            //~ } 
//~ 
            //~ if (HVCRA_local[i].get_fired_dend())
            //~ {
                //~ spikes_in_trial_dend_local[i].push_back(internal_time);
            //~ }
		//~ }
//~ 
		//~ for (int i = 0; i < N_I_local; i++)
		//~ {
            //~ HVCI_local[i].DP8_step_no_target_update();
            //~ 
            //~ //  if some I neuron spikes, change conductance update array
            //~ if (HVCI_local[i].get_fired())
            //~ {
                //~ //printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
                //~ some_I_neuron_fired_local = 1;
                //~ spikes_in_trial_interneuron_local[i].push_back(internal_time);
//~ 
                //~ size_t num_RA_targets = syn_ID_I_RA_local[i].size();
                //~ // loop over all targets of fired neurons
                //~ for (size_t j = 0; j < num_RA_targets; j++)
                //~ {
                    //~ int syn_ID = syn_ID_I_RA_local[i][j];
                    //~ update_Gi_RA_local[syn_ID] += weights_I_RA_local[i][j];
                    //~ //printf("Rank = %d; i = %d; update_Gi_RA_local[%d] = %f; weights_I_RA_local[%d][%d] = %f\n", MPI_rank, i, syn_ID,
                     //~ //   update_Gi_RA_local[syn_ID], weights_I_RA_local[fired_ID][j], fired_ID, j);
                //~ }
            //~ }
		//~ }
//~ 
        //~ // if we need to update network state
        //~ // get if any neurons fired in some process
//~ 
        //~ if (internal_time > network_time)
        //~ {
            //~ MPI_Allreduce(&some_RA_neuron_fired_soma_local, &some_RA_neuron_fired_soma_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            //~ MPI_Allreduce(&some_I_neuron_fired_local, &some_I_neuron_fired_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        //~ 
            //~ if (some_I_neuron_fired_global > 0)
            //~ {
            //~ // sum update array and send to all processes
//~ 
                //~ MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//~ 
                //~ for (int i = 0; i < N_RA_local; i++)
                //~ {
                    //~ HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
                //~ }
//~ 
				//~ // reset fired indicators and arrays
				//~ some_I_neuron_fired_global = 0;
				//~ some_I_neuron_fired_local = 0;
//~ 
//~ 
				//~ std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				//~ std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            //~ }
//~ 
            //~ //if (some_RA_neuron_fired_global == 1)
            //~ //    printf("Rank %d; some_RA_neuron_fired_global: %d\n", MPI_rank, some_RA_neuron_fired_global);
//~ 
            //~ // if somatic compartment of any neuron in the pool fired, update synaptic conductances
             //~ if (some_RA_neuron_fired_soma_global > 0)
             //~ {
            //~ // sum all update arrays and send to all processes
//~ 
                //~ MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                //~ MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//~ 
                //~ // now update excitatory conductances of all neurons
                //~ for (int i = 0; i < N_RA_local; i++)
                //~ {
                    //~ HVCRA_local[i].raiseE(update_Ge_RA_global[Id_RA_local[i]]); // update conductance
				//~ }
//~ 
                //~ for (int i = 0; i < N_I_local; i++)
                //~ {
                    //~ HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);
                //~ }
//~ 
				//~ // reset conductance arrays and fired indicators
				//~ some_RA_neuron_fired_soma_global = 0;
				//~ some_RA_neuron_fired_soma_local = 0;
//~ 
				//~ std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				//~ std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
				//~ std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				//~ std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
            //~ }
//~ 
            //~ network_time += network_update_frequency;
        //~ } // end network update
//~ 
//~ 
        //~ //MPI_Barrier(MPI_COMM_WORLD);
    //~ }
 //~ // end evolve dynamics
//~ }


//~ 
//~ void NetworkGrowthSimulator::run_trials_no_save(int num_trials)
//~ {
    //~ bool training = true;
//~ 
	//~ while (trial_number < num_trials)
    //~ {
        //~ this->trial(training);
//~ 
        //~ trial_number++;
        //~ this->randomize_after_trial();
    //~ }
//~ }
//~ 
//~ 
//~ void NetworkGrowthSimulator::run_trials_with_save(int num_trials)
//~ {
	//~ std::string fileNumSynapses = outputDirectory + "num_synapses.bin";
	//~ 
    //~ std::string fileActiveGraph = outputDirectory + "RA_RA_connections.bin";
    //~ std::string fileSuperGraph = outputDirectory + "RA_RA_super_connections.bin";
    //~ std::string fileTimeSoma = outputDirectory + "spike_times_soma.bin";
    //~ std::string fileTimeDend = outputDirectory + "spike_times_dend.bin";
    //~ std::string fileTimeInterneuron = outputDirectory + "spike_times_interneuron.bin";
    //~ std::string fileWeightsGraph = outputDirectory + "weights.bin";
    //~ std::string fileMaturationGraph = outputDirectory + "mature.bin";
//~ 
    //~ bool training = true;
//~ 
	//~ while (trial_number < num_trials)
    //~ {
        //~ this->trial(training);
//~ 
        //~ this->gather_data();
//~ 
        //~ this->write_num_synapses(fileNumSynapses.c_str());
        //~ this->write_soma_spike_times(fileTimeSoma.c_str());
        //~ this->write_dend_spike_times(fileTimeDend.c_str());
        //~ this->write_interneuron_spike_times(fileTimeInterneuron.c_str());
       //~ 
        //~ this->write_supersynapses(fileSuperGraph.c_str());
        //~ this->write_maturation_info(fileMaturationGraph.c_str());
    //~ 
        //~ this->write_weights(fileWeightsGraph.c_str());
        //~ this->write_active_synapses(fileActiveGraph.c_str());
        //~ this->write_maturation_info((outputDirectory + "mature" + std::to_string(trial_number) + ".bin").c_str());
        //~ 
        //~ trial_number++;
        //~ this->randomize_after_trial();
    //~ }
//~ 
//~ }

void NetworkGrowthSimulator::chain_growth_with_inhibition_tracking(bool training, int save_freq_long, int num_trials, 
																double time_resolution_conductance, std::string outputDirectory)
{	
	std::vector<double> spread_times(N_RA);
	
	if (MPI_rank == 0){
		for (int i = 0; i < N_TR; i++){
			spread_times[training_neurons[i]] = training_spread_times[i];
			
			std::cout << "Neuron " << training_neurons[i] << " spread time = " << spread_times[training_neurons[i]] << std::endl;
		}
	}
	
	MPI_Bcast(&spread_times[0], N_RA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&training_spread_times[0], N_TR, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (trial_number = 0; trial_number < num_trials; trial_number++)
    {
        if (MPI_rank == 0)
            std::cout << "Trial " << trial_number << std::endl;
        
	    this->trial_1stSoma_pre_1stSoma_post_delays_fixedSpread_with_inhibition_tracking(training, spread_times, time_resolution_conductance);
	    
	    int bad_values_indicator_local = this->check_bad_values();
	    int bad_values_indicator_global;
	    
	    MPI_Allreduce(&bad_values_indicator_local, &bad_values_indicator_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	   
		if (bad_values_indicator_global < 0)
		{
			std::cout << "Bad values in neuron dynamics!" << std::endl;
			return;
        }  	

		

		if (trial_number % save_freq_long == 0)
		{
			this->gather_graph_state_data();
			this->gather_full_state_data();
			
			if (MPI_rank == 0)
				this->write_full_network_state("_" + std::to_string(trial_number), outputDirectory);
		}
		
		this->gather_inhibition();
		this->gather_graph_state_data();

		if (MPI_rank == 0)
			this->write_inhibition_tracking_state(trial_number, time_resolution_conductance, outputDirectory);
			
		
		// advance internal time of each neuron
		for (int i = 0; i < N_RA_local; i++)
			num_trials_after_replacement_local[i] += 1;
        
		std::vector<int> neurons_to_replace;
		 
		this->check_neuron_activity(neurons_to_replace);
		this->update_neuron_properties_sameDendrite_diffMaturationRate();
		 
		if ( !neurons_to_replace.empty() )
			this->replace_neurons(neurons_to_replace);
			
		this->setToRest_afterEpoch();
    }
}


void NetworkGrowthSimulator::chain_growth(bool training, int save_freq_short, int save_freq_long, std::string outputDirectory)
{
    std::vector<int> RAtoWrite{692, 40, 570, 4, 103, 482, 67, 454};
    std::vector<int> ItoWrite{0, 1, 2, 3};
    

    std::vector<int> source{0, 1, 2, 3};
    std::vector<int> target{};

    for (int i = 0; i < N_RA; i++)
        target.push_back(i);

    // make all neurons able to make output connections
    //this->set_all_mature();

	// set recording
	//this->set_recording(RAtoWrite, ItoWrite, outputDirectory);

	
	std::vector<double> spread_times(N_RA);
	
	if (MPI_rank == 0){
		for (int i = 0; i < N_TR; i++){
			spread_times[training_neurons[i]] = training_spread_times[i];
			
			std::cout << "Neuron " << training_neurons[i] << " spread time = " << spread_times[training_neurons[i]] << std::endl;
		}
	}
	
	MPI_Bcast(&spread_times[0], N_RA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&training_spread_times[0], N_TR, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    while (true)
    {
       /* if (trial_number == trial_to_add_neurons)
        {
            if (MPI_rank == 0)
                std::cout << "Adding new neurons!" << std::endl;
            
            this->add_new_neurons(N);
        }
        */
		//break;
		
        if (MPI_rank == 0)
            std::cout << "Trial " << trial_number << std::endl;
        
		//this->trial_burst_stdp(training);
		//this->trial_somatic_stdp_no_delays(training);
		
		//this->set_neuron_properties();
		//this->trial_soma_pre_dend_post_stdp_no_delays(training);
	    //this->trial_soma_pre_dend_post_stdp_delays(training);
	    //this->trial_event_pre_dend_post_delays_sudden_maturation(training);
	    //this->trial_burst_pre_dend_post_delays_sudden_maturation_noImmatureOut(training);
	    
	    //this->trial_noImmatureOut_fixedSpread(training, spread_times);
	    
	    //this->trial_burst_pre_dend_event_post_delays_sudden_maturation_noImmatureOut_fixedSpread(training, spread_times);
	    //this->trial_soma_pre_dend_post_stdp_delays_sudden_maturation(training);
	    
	    this->trial_1stSoma_pre_1stSoma_post_delays_fixedSpread(training, spread_times);
	    
	    int bad_values_indicator_local = this->check_bad_values();
	    int bad_values_indicator_global;
	    
	    MPI_Allreduce(&bad_values_indicator_local, &bad_values_indicator_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	   
		if (bad_values_indicator_global < 0)
		{
			std::cout << "Bad values in neuron dynamics!" << std::endl;
			return;
        }  
	    
        if (trial_number % save_freq_short == 0)
        {
			this->gather_graph_state_data();

			if (MPI_rank == 0)
				this->write_graph_network_state(outputDirectory);
           	
           	
           	//this->write_maturation_info(fileMaturationGraph.c_str());
		
            //this->write_weights_time_sequence_from_source_to_target(source, target, fileWeightsSourceToTarget.c_str());
            //this->write_maturation_time_sequence(target, fileMaturationTimeSequence.c_str());
        }
	
			

		if (trial_number % save_freq_long == 0)
		{
			
			this->gather_graph_state_data();
			this->gather_full_state_data();
			
			if (MPI_rank == 0)
				this->write_full_network_state("_" + std::to_string(trial_number), outputDirectory);
		}
		
		//break;    
		//for (int i = 0; i < N_I; i++)
	    //{
	      //  fileAllIneurons = Idir + "I" + std::to_string(i) + ".bin";
		//	pool.write_I(fileAllIneurons.c_str(), i);
	    //}
	    
		//this->reset_after_trial();
		
		// advance internal time of each neuron
		for (int i = 0; i < N_RA_local; i++)
			num_trials_after_replacement_local[i] += 1;
        
		
		 
		std::vector<int> neurons_to_replace;
		 
		this->check_neuron_activity(neurons_to_replace);
		
		
		this->update_neuron_properties_sameDendrite_diffMaturationRate();
		 
		if ( !neurons_to_replace.empty() )
		{
			// write network state before replacement
			//this->gather_graph_state_data();
			//this->gather_full_state_data();
			
			//if (MPI_rank == 0)
				//this->write_full_network_state("_" + std::to_string(trial_number) + "beforeReplacement", outputDirectory);
	
			this->replace_neurons(neurons_to_replace);
			
			// write network state after replacement
			//this->gather_graph_state_data();
			//this->gather_full_state_data();
			
			//if (MPI_rank == 0)
			//{
				//this->write_full_network_state("_" + std::to_string(trial_number) + "afterReplacement", outputDirectory);
				//this->write_replaced_neurons(neurons_to_replace, (outputDirectory + "replaced_neurons.bin").c_str());
			//}
		}
		
		this->setToRest_afterEpoch();
		
		//this->update_neuron_properties();
		
		
		//this->reset_after_trial_soma_pre_dend_post_no_delays();
		//this->reset_after_trial_soma_pre_dend_post_delays();
		
		//this->setToRest_after_trial_soma_pre_dend_post_delays();
		
		
		// replace neurons if needed
		
		//~ // gather all neurons that are to be replaced
		//~ this->gather_neurons_2replace();
    //~ 
		//~ // if some neurons are to be replaced, replace them
		//~ if (replace_real_id_global.size() > 0)
			//~ this->replace_neurons();
			
        trial_number++;
        
        //this->randomize_after_trial();
		//break;
    }
}


/*void NetworkGrowthSimulator::trial_somatic_stdp_no_delays(bool training)
{
	// spike indicators
	int some_RA_neuron_spiked_local;
	int some_I_neuron_fired_local;
	
	int some_RA_neuron_spiked_global;
	int some_I_neuron_fired_global;

    std::vector<double> spike_times_fired_dend_local;
    
    std::vector<int> RA_neurons_spiked_global;
    std::vector<int> RA_neurons_spiked_local;
    
    // arrays with updates for conductances
    std::vector<double> update_Ge_RA_local(N_RA);
	std::vector<double> update_Ge_RA_global(N_RA);

	std::vector<double> update_Gi_RA_local(N_RA);
	std::vector<double> update_Gi_RA_global(N_RA);

	std::vector<double> update_Ge_I_local(N_I);
	std::vector<double> update_Ge_I_global(N_I);

	// sample training innervation time and send to all processes
    double training_kick_time;
    
	if (MPI_rank == 0)
	{
		training_kick_time = WAITING_TIME + noise_generator.random(TRIAL_DURATION - 2*WAITING_TIME);
		std::cout << "training_kick_time = " << training_kick_time << std::endl;
    }
    MPI_Bcast(&training_kick_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    bool training_excited = false; // indicator that training neurons were already excited
    
    double internal_time = 0;
    double network_time = NETWORK_UPDATE_FREQUENCY;

   

    //printf("network time = %f\n", network_time);

    
	// initialize update arrays and fired indicators
	std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);

	some_RA_neuron_spiked_local = 0;
    some_RA_neuron_spiked_global = 0;

    some_I_neuron_fired_local = 0;
    some_I_neuron_fired_global = 0;
	
	this->set_neuron_properties();
	
	
    // evolve dynamics
    for (int t = 1; t < static_cast<int>(TRIAL_DURATION / TIMESTEP); t++)
	{
		internal_time += TIMESTEP;
	
		if ( ( training ) && ( !training_excited ) && ( internal_time >= training_kick_time ) )
		{
			for (int i = 0; i < N_TR; i++)
			{
				int rank;
				int shift;
				
				this->get_neuronRA_location(training_neurons[i], &rank, &shift);
				
				if (MPI_rank == rank)
					HVCRA_local[shift].raiseE(G_TRAINING_KICK);
			}
			
			training_excited = true;
		}
		
		for (int i = 0; i < N_RA_local; i++)
		{
            // Debraband step
            HVCRA_local[i].Debraband_step_no_target_update();
            
            // if some neuron produced somatic spike, do LTD for all previous dendritic spikes
            if ( HVCRA_local[i].get_fired_soma() )
            {
				
				// show spikes of pool neurons only
                auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
                
                if (p.first == p.second)
					std::cout << "HVC-RA " << Id_RA_local[i] << " soma spike at " << internal_time << std::endl;
                	
				//std::cout << "HVC-RA " << Id_RA_local[i] << " soma spike at " << internal_time << std::endl;
                	
				
                spikes_in_trial_soma_local[i].push_back(internal_time);
				//previous_soma_spike_times_local[i].push_back(internal_time);
				
				RA_neurons_spiked_local.push_back(Id_RA_local[i]);
                
				some_RA_neuron_spiked_local = 1;
				
				// loop over all inhibitory targets of fired neurons
				size_t num_I_targets = syn_ID_RA_I_local[i].size();
				for (size_t j = 0; j < num_I_targets; j++)
				{
					int syn_ID = syn_ID_RA_I_local[i][j];
					update_Ge_I_local[syn_ID] += weights_RA_I_local[i][j];

				}
				
				
				// loop over all excitatory targets
				size_t num_RA_targets = active_synapses_local[i].size();
				for (size_t j = 0; j < num_RA_targets; j++)
				{
					int syn_ID = active_synapses_local[i][j];
					update_Ge_RA_local[syn_ID] += weights_RA_RA_local[i][syn_ID];
				}

				 // if neuron is saturated apply LTD only to supersynapses
                if (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss)
                {
                    for (size_t k = 0; k < supersynapses_local[i].size(); k++)
                    {
                        int supersynapse_id = supersynapses_local[i][k];
                        
                        if ( Id_RA_local[i] != supersynapse_id )
                        {
							// find previous somatic spikes that are not too old
							if ( !previous_somatic_spike_times_global[supersynapse_id].empty() )
							{
								std::vector<double>::iterator it_relevant_spikes = std::lower_bound(previous_somatic_spike_times_global[supersynapse_id].begin(),
																									previous_somatic_spike_times_global[supersynapse_id].end(),
																									internal_time - STDP_WINDOW);
																					
								for (auto it = it_relevant_spikes; it != previous_somatic_spike_times_global[supersynapse_id].end(); it++)
								{
									double dt = *it - internal_time;

									//std::cout << "dt in saturated LTD = " << dt << std::endl;

									
									double w_before = weights_RA_RA_local[i][supersynapse_id];     
						   
									LTD(weights_RA_RA_local[i][supersynapse_id], dt);
									
									
									//~ std::cout << "LTD from saturated " << Id_RA_local[i] << " -> " << supersynapse_id
											  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][supersynapse_id] - w_before
											  //~ << std::endl;
											  
									//printf("LTD from saturated neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n",
									//            Id_RA_local[i], supersynapse_id, spikes_in_trial_soma_local[i].back(), spike_times_dend_global[supersynapse_id],
									//            dt, weights_local[i][supersynapse_id] - w);
									
									update_synapse(i, supersynapse_id);	
									
								}
							}
                           
                        }
                    }

					// if some supersynapse desaturated, update all synapses
                	if (static_cast<int>(supersynapses_local[i].size()) < synaptic_params.Nss)
						for (int j = 0; j < N_RA; j++)
							this->update_synapse(i, j);
                }
                // if not saturated apply LTD rule with the last dendritic spike of all neurons and add glutamate to all neuron except the fired one
                else
                {
                    for (int j = 0; j < N_RA; j++)
                    {
                        if ( Id_RA_local[i] != j )
                        {
							// find previous somatic spikes that are not too old
							if ( !previous_somatic_spike_times_global[j].empty() )
							{
								std::vector<double>::iterator it_relevant_spikes = std::lower_bound(previous_somatic_spike_times_global[j].begin(),
																									previous_somatic_spike_times_global[j].end(),
																									internal_time - STDP_WINDOW);
																					
								for (auto it = it_relevant_spikes; it != previous_somatic_spike_times_global[j].end(); it++)
								{
									double dt = *it - internal_time;
                            
									//std::cout << "dt in LTD = " << dt << std::endl;
									double w_before = weights_RA_RA_local[i][j];

									LTD(weights_RA_RA_local[i][j], dt);
									
									//~ // show LTD of pool neurons only
									//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
									//~ 
									//~ if (p.first == p.second)
									//~ {
										//~ std::cout << "LTD from  " << Id_RA_local[i] << " -> " << j
											  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][j] - w_before
											  //~ << std::endl;
									//~ }
									
									
									//~ std::cout << "LTD from  " << Id_RA_local[i] << " -> " << j
											  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][j] - w_before
											  //~ << std::endl;
									//printf("LTD from neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n", Id_RA_local[i], j,
									//	spikes_in_trial_soma_local[i].back(), spike_times_dend_global[j], dt, weights_local[i][j] - w);
									
									update_synapse(i, j);
											
									
								}
							}
                        }
                    }
                }

                //for (int j = 0; j < last_soma_spikes_local[i].size(); j++)
                //{
                  //  printf("My rank = %d; Soma spikes of neuron %d:  spike_time = %f\n", MPI_rank, Id_RA_local[i],
                    //    last_soma_spikes_local[i][j]);

                //}

				//std::cout << "HVC-RA " << Id_RA_local[i] << " soma spike at " << internal_time << std::endl;
                
                
              
                

            } // end if get fired soma
            
            // if some neuron produced dendritic spike, store this neuron in array
            if ( HVCRA_local[i].get_fired_dend() )
            {
				// show spikes of pool neurons only
                auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
                
                if (p.first == p.second)
					std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
                
				
                spikes_in_trial_dend_local[i].push_back(internal_time); // dend spike time relative to trial onset
                
                //std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
                
                
                //printf("My rank = %d; RA neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], internal_time);
               
                //printf("My rank = %d; Dend neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], spike_times_dend_local[i]);

               
            } // end if get_fired_dend

		} // end for i -> N_RA_local

		for (int i = 0; i < N_I_local; i++)
		{
            HVCI_local[i].DP8_step_no_target_update();
            
            //  if some I neuron spikes, change conductance update array
            if (HVCI_local[i].get_fired())
            {
                //printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
                some_I_neuron_fired_local = 1;
                spikes_in_trial_interneuron_local[i].push_back(internal_time);

                size_t num_RA_targets = syn_ID_I_RA_local[i].size();
                // loop over all targets of fired neurons
                for (size_t j = 0; j < num_RA_targets; j++)
                {
                    int syn_ID = syn_ID_I_RA_local[i][j];
                    update_Gi_RA_local[syn_ID] += weights_I_RA_local[i][j];
                    //printf("Rank = %d; i = %d; update_Gi_RA_local[%d] = %f; weights_I_RA_local[%d][%d] = %f\n", MPI_rank, i, syn_ID,
                     //   update_Gi_RA_local[syn_ID], weights_I_RA_local[fired_ID][j], fired_ID, j);
                }
                
                //std::cout << "HVC-I " << Id_I_local[i] << " spike at " << internal_time << std::endl;
                
            }
		} // end if i -> N_I_local

        // if we need to update network state or if we reached the end of the trial
        // get if any neurons fired in some process

        if ( ( internal_time > network_time ) || (t == size-1) )
        {
            //if (MPI_rank == 0)
            //{
            //    std::cout << "internal time = " << internal_time << " > network_time = " << network_time << std::endl; 
            //}

            MPI_Allreduce(&some_RA_neuron_spiked_local, &some_RA_neuron_spiked_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_I_neuron_fired_local, &some_I_neuron_fired_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            if (some_I_neuron_fired_global > 0)
            {
            // sum update array and send to all processes

                MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                for (int i = 0; i < N_RA_local; i++)
                    HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
            	
				// update conductance arrays and fired indicators
				some_I_neuron_fired_local = 0;
            	some_I_neuron_fired_global = 0;
				
				std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            }


            //if (some_RA_neuron_fired_global == 1)
            //    printf("Rank %d; some_RA_neuron_fired_global: %d\n", MPI_rank, some_RA_neuron_fired_global);

            // if somatic compartment of any neuron in the pool fired, update synaptic conductances
             if ( some_RA_neuron_spiked_global > 0 )
             {
            // sum all update arrays and send to all processes

                MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                // now update excitatory conductances of all neurons
                for (int i = 0; i < N_RA_local; i++)
                    HVCRA_local[i].raiseExcWeight(update_Ge_RA_global[Id_RA_local[i]]); // update conductance
				
                for (int i = 0; i < N_I_local; i++)
                    HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);
                
				// update conductance arrays and fired indicators
            	some_RA_neuron_spiked_local = 0;
	        	some_RA_neuron_spiked_global = 0;
				
				std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
				
				std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
            
				// gather all spiked HVC-RA neurons
				this->gather_spiked_or_bursted_neurons(RA_neurons_spiked_local, RA_neurons_spiked_global);
				
				//~ if (MPI_rank == 0)
				//~ {
					//~ std::cout << "Spiked neurons:\n";
					//~ 
					//~ for (size_t j = 0; j < RA_neurons_spiked_global.size(); j++)
						//~ std::cout << RA_neurons_spiked_global[j] << " ";
					//~ std::cout << std::endl;
				//~ }
				
				// add somatic spikes for spiked neurons
				for (size_t j = 0; j < RA_neurons_spiked_global.size(); j++)
					previous_somatic_spike_times_global[RA_neurons_spiked_global[j]].push_back(internal_time);
                  
              	/* 
                if (MPI_rank == 0)
                {
					for (size_t i = 0; i < RA_neurons_fired_dend_global.size(); i++)
                        printf("neuron %d bursted; spike_time_dend = %f\n", RA_neurons_fired_dend_global[i], 
								spikes_in_trial_dend_global[RA_neurons_fired_dend_global[i]].back());
                }
                */

			/*
                // apply LTP rule to the previous somatic spikes of pool neurons and dendritic spike of fired neurons

                for (int i = 0; i < N_RA_local; i++)
                {
                    // if neuron is saturated apply LTP only if spiked neurons are among supersynapse targets
                    if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
                    {
                        for (size_t j = 0; j < RA_neurons_spiked_global.size(); j++)
                        {
                            int fired_ID = RA_neurons_spiked_global[j];
                            
                            // do not allow self-to-self synapses to emerge
							if ( Id_RA_local[i] != fired_ID )
							{
								std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
											supersynapses_local[i].end(), fired_ID);

								if ( pos!=supersynapses_local[i].end() )
								{
									// find previous somatic spikes that are not too old
									if ( !previous_somatic_spike_times_global[Id_RA_local[i]].empty() )
									{
										std::vector<double>::iterator it_relevant_spikes = std::lower_bound(previous_somatic_spike_times_global[Id_RA_local[i]].begin(),
																											previous_somatic_spike_times_global[Id_RA_local[i]].end(),
																											internal_time - STDP_WINDOW);
																							
										for (auto it = it_relevant_spikes; it != previous_somatic_spike_times_global[Id_RA_local[i]].end(); it++)
										{
											double dt = internal_time - *it;
									
										
											//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
											
											if (dt <= synaptic_params.T_0)
											{
												double w_before = weights_RA_RA_local[i][fired_ID];
												
												LTD(weights_RA_RA_local[i][fired_ID], dt);
												
												//~ std::cout   << "LTD from saturated " << Id_RA_local[i] << " -> " << fired_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][fired_ID] - w_before
															//~ << std::endl;
											}
											else
											{
												double w_before = weights_RA_RA_local[i][fired_ID];
												
												
												LTP(weights_RA_RA_local[i][fired_ID], dt);
												
												//~ std::cout   << "LTP from saturated " << Id_RA_local[i] << " -> " << fired_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][fired_ID] - w_before
															//~ << std::endl;
												
											}	
											//double w = weights_local[i][fired_ID];
											
											//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
											 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
											 //           dt, weights_local[i][fired_ID] - w);
											
											update_synapse(i, fired_ID);
											
										}
									}
								}
							}
                        }

                    }
                    // if not saturated apply LTP for all dendritic spikes
                    else
                    {
                        for (size_t j = 0; j < RA_neurons_spiked_global.size(); j++)
                        {
                            int fired_ID = RA_neurons_spiked_global[j];
                            // don't allow self-to-self connections
                            if ( fired_ID != Id_RA_local[i] )
                            {
                                // find previous somatic spikes that are not too old
								if ( !previous_somatic_spike_times_global[Id_RA_local[i]].empty() )
								{
									std::vector<double>::iterator it_relevant_spikes = std::lower_bound(previous_somatic_spike_times_global[Id_RA_local[i]].begin(),
																										previous_somatic_spike_times_global[Id_RA_local[i]].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it = it_relevant_spikes; it != previous_somatic_spike_times_global[Id_RA_local[i]].end(); it++)
									{
										double dt = internal_time - *it;
								

										//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
										
									    if (dt <= synaptic_params.T_0)
									    {
											double w_before = weights_RA_RA_local[i][fired_ID];
												
											LTD(weights_RA_RA_local[i][fired_ID], dt);
											
											//~ // show LTD on pool neurons only
											//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), fired_ID);
											//~ 
											//~ if (p.first == p.second)
											//~ {
												//~ std::cout << "LTD from  " << Id_RA_local[i] << " -> " << fired_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][fired_ID] - w_before
													  //~ << std::endl;
											//~ }
											
											//~ std::cout   << "LTD from " << Id_RA_local[i] << " -> " << fired_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][fired_ID] - w_before
															//~ << std::endl;
											
										}
										else
										{
											double w_before = weights_RA_RA_local[i][fired_ID];
											
											LTP(weights_RA_RA_local[i][fired_ID], dt);
											
											//~ // show LTP on pool neurons only
											//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), fired_ID);
											//~ 
											//~ if (p.first == p.second)
											//~ {
												//~ std::cout << "LTP from  " << Id_RA_local[i] << " -> " << fired_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][fired_ID] - w_before
													  //~ << std::endl;
											//~ }
											
											
											//~ std::cout   << "LTP from " << Id_RA_local[i] << " -> " << fired_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][fired_ID] - w_before
															//~ << std::endl;
										}	
										//double w = weights_local[i][fired_ID];
										
										//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
										 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
										 //           dt, weights_local[i][fired_ID] - w);
										
										update_synapse(i, fired_ID);
									
									}
									
									
								}
                                
                            }
                        }
                   }

                } // end for i -> N_RA_local

				       

               // check if we need axon remodeling

                for (int i = 0; i < N_RA_local; i++)
                {
                    if ( (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss) && (remodeled_local[i] == 0) )
                    {
					    this->axon_remodeling(i);

				    }
                }
            
				// clear arrays
				RA_neurons_spiked_local.clear();
				RA_neurons_spiked_global.clear();
            
            } // end if some HVC-RA neuron spiked

           
            network_time += NETWORK_UPDATE_FREQUENCY;
            
            //printf("network_time = %f\n", network_time);
        }


        //MPI_Barrier(MPI_COMM_WORLD);
    }
    this->potentiation_decay();
    //printf("After potentiation decay")
    this->update_all_synapses();
	
	// calculate new gaba reverse potential
	// update rate
	for (int i = 0; i < N_RA_local; i++)
	{
		num_spikes_in_recent_trials_local[i].push_front(static_cast<int>(spikes_in_trial_soma_local[i].size()));
        
        //firing_rate_short_local[i] = std::accumulate(num_spikes_in_recent_trials[i].begin(), num_spikes_in_recent_trials[i].begin() + RATE_WINDOW_SHORT, 0.0)
          //                          / static_cast<double>(RATE_WINDOW_SHORT);

		// calculate firing rate in large window:
        firing_rate_long_local[i] = std::accumulate(num_spikes_in_recent_trials_local[i].begin(), num_spikes_in_recent_trials_local[i].end(), 0.0)
                                                                                            / static_cast<double>(maturation_params.RATE_WINDOW_LONG);
                                                                                            
		//if (std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].end(), 0.0) >= 1)
		//{
		//	std::cout << "Recent bursts of neuron " << Id_RA_local[i] << " :" << std::endl;
			
		//	for (size_t j = 0; j < num_bursts_in_recent_trials[i].size(); j++)
		//		std::cout << num_bursts_in_recent_trials[i][j] << "\t";
			
		//	std::cout << std::endl;
		//}
	}
    // update GABA potential based on firing rates
	//this->update_Ei();
    
    // gather all neurons that are to be replaced
    //this->gather_neurons_2replace();
    
    // if some neurons are to be replaced, replace them
    //if (replace_real_id_global.size() > 0)
    //    this->replace_neurons();
    
    
     // advance internal time of each neuron
    for (int i = 0; i < N_RA_local; i++)
        num_trials_after_replacement_local[i] += 1;
        
    this->update_neuron_properties();
    
}*/

void NetworkGrowthSimulator::trial_burst_pre_dend_event_post_delays_sudden_maturation_noImmatureOut(bool training)
{
	double event_window = 15.0; // window in which all somatic spikes are considered as one event
	
	// indicators for conductance updates and bursting
	int some_RA_inh_conductance_was_updated_global = 0;
	int some_RA_inh_conductance_was_updated_local = 0;
	
	int some_RA_exc_conductance_was_updated_global = 0;
	int some_RA_exc_conductance_was_updated_local = 0;
	
	int some_I_exc_conductance_was_updated_global = 0;
	int some_I_exc_conductance_was_updated_local = 0;
	
	int some_RA_neuron_bursted_global = 0;
	int some_RA_neuron_bursted_local = 0;
	
	 // arrays with updates for conductances
    std::vector<double> update_Ge_RA_local(N_RA);
	std::vector<double> update_Ge_RA_global(N_RA);

	std::vector<double> update_Gi_RA_local(N_RA);
	std::vector<double> update_Gi_RA_global(N_RA);

	std::vector<double> update_Ge_I_local(N_I);
	std::vector<double> update_Ge_I_global(N_I);
	
	std::vector<int> RA_neurons_bursted_local; // local array of HVC-RA neurons that bursted
	std::vector<int> RA_neurons_bursted_global; // global array of HVC-RA neurons that bursted
	
	// sample training innervation time and send to all processes
    double training_kick_time;
    
	if (MPI_rank == 0)
	{
		training_kick_time = WAITING_TIME + noise_generator.random(TRIAL_DURATION - 2*WAITING_TIME);
		std::cout << "training_kick_time = " << training_kick_time << std::endl;
    }
    MPI_Bcast(&training_kick_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    bool training_excited = false; // indicator that training neurons were already excited
    
    double internal_time = 0;
    double network_time = NETWORK_UPDATE_FREQUENCY;

   

    //printf("network time = %f\n", network_time);

    
	// initialize update arrays and fired indicators
	std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
	
	
    // evolve dynamics
    int num_steps = static_cast<int>(TRIAL_DURATION / TIMESTEP);
    
    for (int t = 1; t < num_steps; t++)
	{
		internal_time += TIMESTEP;
	
		if ( ( training ) && ( !training_excited ) && ( internal_time >= training_kick_time ) )
		{
			for (int i = 0; i < N_TR; i++)
			{
				int rank;
				int shift;
				
				this->get_neuronRA_location(training_neurons[i], &rank, &shift);
				
				if (MPI_rank == rank)
					HVCRA_local[shift].raiseE(G_TRAINING_KICK);
			}
			
			training_excited = true;
		}
		
		//////////////////////
		// HVC-RA neurons
		//////////////////////
		
		for (int i = 0; i < N_RA_local; i++)
		{
            // Debraband step
            HVCRA_local[i].Debraband_step_no_target_update();
            
           
			///////////////////////////////
			// HVC-RA somatic spike
			///////////////////////////////
			// if some neuron produced somatic spike
            if (HVCRA_local[i].get_fired_soma())
            {
                spikes_in_trial_soma_local[i].push_back(internal_time);
                
				// update delivery queues if neuron is mature (no output for immature neuron and no wiring mechanism)
				if ( mature_global[Id_RA_local[i]] == 1 ){
				
					// for inhibitory neurons
					// loop over all inhibitory targets of fired neurons
					size_t num_I_targets = syn_ID_RA_I_local[i].size();
					
					for (size_t j = 0; j < num_I_targets; j++)
					{
						double delivery_time = internal_time + axonal_delays_RA_I_local[i][j];
						
						// if queue is empty, just add item to the queue
						if ( delivery_queue_RA_I[i].empty() )
							delivery_queue_RA_I[i].push_back(std::pair<double,int>(delivery_time, j));	
						// otherwise add item so that queue is sorted
						else
						{
							auto it = std::upper_bound(delivery_queue_RA_I[i].begin(), delivery_queue_RA_I[i].end(), std::pair<double,int>(delivery_time, j));
							delivery_queue_RA_I[i].insert(it, std::pair<double,int>(delivery_time, j));
						}
					}
						//std::cout << "neuron " << Id_RA_local[i] << " spike time = " << internal_time << " axonal delay = " << axonal_delays_RA_I[Id_RA_local[i]][syn_ID] << " delivery time RA to I: " << delivery_queue_RA_I[i].back().first << " delivery target id: " << syn_ID_RA_I_local[i][delivery_queue_RA_I[i].back().second] << std::endl;
					
					
					// if neuron is saturated, deliver spikes only to supersynaptic targets
					if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
					{
						for (size_t j = 0; j < supersynapses_local[i].size(); j++)
						{
							int target_id = supersynapses_local[i][j];
							
							double delivery_time = internal_time + axonal_delays_RA_RA_local[i][target_id];
							
							// if queue is empty, just add item to the queue
							if ( delivery_queue_RA_RA_soma[i].empty() )
								delivery_queue_RA_RA_soma[i].push_back(std::pair<double,int>(delivery_time, target_id));	
							// otherwise add item so that queue is sorted
							else
							{
								auto it = std::upper_bound(delivery_queue_RA_RA_soma[i].begin(), delivery_queue_RA_RA_soma[i].end(), std::pair<double,int>(delivery_time, target_id));
								delivery_queue_RA_RA_soma[i].insert(it, std::pair<double,int>(delivery_time, target_id));
							}
						}
					}
					else // deliver spikes to everyone except itself
					{	
						for (int j = 0; j < N_RA; j++)
						{
							if ( j != Id_RA_local[i])
							{
							
								double delivery_time = internal_time + axonal_delays_RA_RA_local[i][j];
							
								// if queue is empty, just add item to the queue
								if ( delivery_queue_RA_RA_soma[i].empty() )
									delivery_queue_RA_RA_soma[i].push_back(std::pair<double,int>(delivery_time, j));	
								// otherwise add item so that queue is sorted
								else
								{
									auto it = std::upper_bound(delivery_queue_RA_RA_soma[i].begin(), delivery_queue_RA_RA_soma[i].end(), std::pair<double,int>(delivery_time, j));
									delivery_queue_RA_RA_soma[i].insert(it, std::pair<double,int>(delivery_time, j));
								}
							}
						}
					}
				}
                
            } // end if local HVC-RA neuron spiked  
            
             //////////////////////////////////
			// HVC-RA -> HVC-I delivery queue
			//////////////////////////////////
            // check if somatic spike was delivered to some interneuron
            // loop through the delivery queue to check if current time exceeds the spike delivery time
            std::vector<std::pair<double,int>>::iterator it = delivery_queue_RA_I[i].begin();
			
            for (; it != delivery_queue_RA_I[i].end(); it++)
            {
				if (internal_time >= it->first)
				{
					
					some_I_exc_conductance_was_updated_local = 1;
					
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_RA_I_local[i][pos_in_local_target_array];
					
					update_Ge_I_local[target_id] += weights_RA_I_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(RA) neuron " << Id_RA_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to HVC(I) neuron " << syn_ID_RA_I_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_I[i].erase(delivery_queue_RA_I[i].begin(), it);
			
			//////////////////////////////////////////////////
			// HVC-RA -> HVC-RA somatic spike delivery queue
			//////////////////////////////////////////////////
			// check if somatic spike was delivered to some HVC(RA) neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
            it = delivery_queue_RA_RA_soma[i].begin();
            
            for (; it != delivery_queue_RA_RA_soma[i].end(); it++)
            {
				double delivery_time = it->first;
				
				if ( internal_time >= delivery_time )
				{
					int target_id = it->second;
					
					
					
					// update conductance of target if synapse is active
					if ( active_indicators_local[i][target_id] == 1 )
					{
						some_RA_exc_conductance_was_updated_local = 1;
						update_Ge_RA_local[target_id] += weights_RA_RA_local[i][target_id];
					}
					
					
					
					//////////////
					//// LTD
					//////////////
					// do LTD on dendritic spike time of target neuron if it is not mature neuron
					//~ if ( mature_global[target_id] != 1 )
					//~ {
						// update delivered spikes only if current spike occured later than previous delivered
						// spike + event_window
						bool new_event_occured = false;
					
						if ( !delivered_spike_times[i][target_id].empty() ){
							if ( delivery_time > delivered_spike_times[i][target_id].back() + event_window ){
								delivered_spike_times[i][target_id].push_back(delivery_time);
								new_event_occured = true;
							}
						}
						else{
							delivered_spike_times[i][target_id].push_back(delivery_time);
							new_event_occured = true;
						}
						
						if ( new_event_occured ){
							 // if neuron is saturated apply LTD only if target is among super synapses
							if (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss)
							{
								if ( supersynapses_indicators_local[i][target_id] == 1 )
								{	
									// find previous dendritic spikes that are not too old
									if ( !previous_dendritic_spike_times_global[target_id].empty() )
									{
										std::vector<double>::iterator it_relevant_dend_spikes = std::lower_bound(previous_dendritic_spike_times_global[target_id].begin(),
																											previous_dendritic_spike_times_global[target_id].end(),
																											internal_time - STDP_WINDOW);
																							
										for (auto it_dend_spike = it_relevant_dend_spikes; it_dend_spike != previous_dendritic_spike_times_global[target_id].end(); it_dend_spike++)
										{
											double dt = *it_dend_spike - delivery_time;

											//std::cout << "dt in saturated LTD = " << dt << std::endl;

											
											double w_before = weights_RA_RA_local[i][target_id];     
								   
											LTD(weights_RA_RA_local[i][target_id], dt);
											
											
											//~ std::cout << "LTD from saturated " << Id_RA_local[i] << " -> " << supersynapse_id
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][supersynapse_id] - w_before
													  //~ << std::endl;
													  
											//printf("LTD from saturated neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n",
											//            Id_RA_local[i], supersynapse_id, spikes_in_trial_soma_local[i].back(), spike_times_dend_global[supersynapse_id],
											//            dt, weights_local[i][supersynapse_id] - w);
											
											update_synapse(i, target_id);	
											
										}   
									}
								}

								// if some supersynapse desaturated, update all synapses
								if (static_cast<int>(supersynapses_local[i].size()) < synaptic_params.Nss)
									for (int j = 0; j < N_RA; j++)
										this->update_synapse(i, j);
							}
							// if not saturated apply LTD rule 
							else
							{
								// find previous dendritic spikes that are not too old
								if ( !previous_dendritic_spike_times_global[target_id].empty() )
								{
									std::vector<double>::iterator it_relevant_dend_spikes = std::lower_bound(previous_dendritic_spike_times_global[target_id].begin(),
																										previous_dendritic_spike_times_global[target_id].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it_dend_spike = it_relevant_dend_spikes; it_dend_spike != previous_dendritic_spike_times_global[target_id].end(); it_dend_spike++)
									{
										double dt = *it_dend_spike - delivery_time;

										//std::cout << "dt in saturated LTD = " << dt << std::endl;

										
										double w_before = weights_RA_RA_local[i][target_id];     
							   
										LTD(weights_RA_RA_local[i][target_id], dt);
										
										//~ // show LTD to pool neurons only
										auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), target_id);
										
										if (p.first == p.second)
										{
											std::cout << "LTD from  " << Id_RA_local[i] << " -> " << target_id
												  << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][target_id] - w_before
												  << " w_after = " << weights_RA_RA_local[i][target_id] << std::endl;
										}
									//~ 
										update_synapse(i, target_id);	
										
									}   
								}
							}
						}
					//~ } // end if target neuron is not mature
				} // end if internal_time >= delivery spike time 
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_RA_soma[i].erase(delivery_queue_RA_RA_soma[i].begin(), it);
			
			
			///////////////////////////////
			// HVC-RA dendritic spike
			///////////////////////////////
			//if some neuron produced dendritic spike, store this neuron in array
			
			if ( HVCRA_local[i].get_fired_dend() )
			{
				
				// show spikes of pool neurons only
				auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
				
				if (p.first == p.second)
					std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
				
				
				spikes_in_trial_dend_local[i].push_back(internal_time); // dend spike time relative to trial onset
				
				//if ( mature_global[Id_RA_local[i]] != 1 )
				//{

				// update dendritic spikes only if current spike occured later than previous 
				// dendritic spike + event_window
				
				if ( previous_dendritic_spike_times_global[Id_RA_local[i]].empty() )
				{
					RA_neurons_bursted_local.push_back(Id_RA_local[i]);
					some_RA_neuron_bursted_local = 1;
				}
				else if ( internal_time > previous_dendritic_spike_times_global[Id_RA_local[i]].back() + event_window)
				{
					RA_neurons_bursted_local.push_back(Id_RA_local[i]);
					some_RA_neuron_bursted_local = 1;
				}
				//}
				//std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
			} // end if get_fired_dend
	
            
		} // end for i = 0 -> N_RA_local (loop through all local HVC-RA neurons)
		
		//////////////////////
		// HVC-I neurons
		//////////////////////
		for (int i = 0; i < N_I_local; i++)
		{
			HVCI_local[i].DP8_step_no_target_update();
			
			
			///////////////////////////////
			// HVC-I spike
			///////////////////////////////
			//  if some I neuron spikes, update delivery queue
			if (HVCI_local[i].get_fired())
			{
				//printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
				spikes_in_trial_interneuron_local[i].push_back(internal_time);

				size_t num_RA_targets = syn_ID_I_RA_local[i].size();
				// loop over all targets of fired neurons
				for (size_t j = 0; j < num_RA_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_I_RA_local[i][j];
					
					//std::cout << "HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << syn_ID_I_RA_local[i][j] << " spike time " << internal_time << " delivery time = " << delivery_time << std::endl; 
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_I_RA[i].empty() )
						delivery_queue_I_RA[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_I_RA[i].begin(), delivery_queue_I_RA[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_I_RA[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
			}
			
			//////////////////////////////////
			// HVC-I -> HVC-RA delivery queue
			//////////////////////////////////
			// check if interneuron spike was delivered to some HVC-RA neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
			std::vector<std::pair<double,int>>::iterator it = delivery_queue_I_RA[i].begin();
			
			for (; it != delivery_queue_I_RA[i].end(); it++)
			{
				if (internal_time >= it->first)
				{
					some_RA_inh_conductance_was_updated_local = 1;
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_I_RA_local[i][pos_in_local_target_array];
					
					//std::cout << "Delivered spike HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << target_id << " delivered_time " << internal_time << " delivery time in queue = " << it->first << std::endl; 
					
					
					update_Gi_RA_local[target_id] += weights_I_RA_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(I) neuron " << Id_I_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to neuron " << syn_ID_I_RA_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_I_RA[i].erase(delivery_queue_I_RA[i].begin(), it);
		
		}  // end for i = 0 -> N_I_local (loop through all local HVC-I neurons)
		
		/////////////////////////////
		// Network Synchronization //
		/////////////////////////////
		if ( ( internal_time > network_time ) || ( t == num_steps-1 ) )
        {
            //if (MPI_rank == 0)
            //{
            //    std::cout << "internal time = " << internal_time << " > network_time = " << network_time << std::endl; 
            //}

            MPI_Allreduce(&some_RA_exc_conductance_was_updated_local, &some_RA_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_RA_inh_conductance_was_updated_local, &some_RA_inh_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_I_exc_conductance_was_updated_local, &some_I_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_RA_neuron_bursted_local, &some_RA_neuron_bursted_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            
        
			
            ////////////////////////////////////////////
            //// Process HVC-I -> HVC-RA interactions
            ////////////////////////////////////////////
            if ( some_RA_inh_conductance_was_updated_global > 0 )
            {
            // sum update array and send to all processes

                MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                for (int i = 0; i < N_RA_local; i++)
                {
					//std::cout << "HVC-RA neuron " << Id_RA_local[i] << " raised inhibitory conductance by " << update_Gi_RA_global[Id_RA_local[i]] << std::endl;
                    HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
            	}
				// update conductance arrays and delivered indicators
				some_RA_inh_conductance_was_updated_global = 0;
				some_RA_inh_conductance_was_updated_local = 0;
				
				std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            }
            
            ////////////////////////////////////////////
            //// Process HVC-RA -> HVC-I interactions
            ////////////////////////////////////////////
			if ( some_I_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons

				for (int i = 0; i < N_I_local; i++)
					HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);

				// update conductance arrays and fired indicators
				some_I_exc_conductance_was_updated_local = 0;
				some_I_exc_conductance_was_updated_global = 0;

				std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
			}
			
			////////////////////////////////////////////
            //// Process HVC-RA -> HVC-RA interactions
            ////////////////////////////////////////////
			if ( some_RA_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons
				for (int i = 0; i < N_RA_local; i++)
					HVCRA_local[i].raiseExcWeight(update_Ge_RA_global[Id_RA_local[i]]); // update conductance

				// update conductance arrays and fired indicators
				some_RA_exc_conductance_was_updated_global = 0;
				some_RA_exc_conductance_was_updated_local = 0;
				
				std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
			}	

			////////////////////////////////////////////
            //// Process HVC-RA dendritic spikes
            ////////////////////////////////////////////
			if ( some_RA_neuron_bursted_global > 0 )
            {
				// update conductance arrays and fired indicators
            	some_RA_neuron_bursted_local = 0;
	        	some_RA_neuron_bursted_global = 0;
	        	
	        	// gather all bursted HVC-RA neurons
				this->gather_spiked_or_bursted_neurons(RA_neurons_bursted_local, RA_neurons_bursted_global);
				
				// add dendritic spikes for bursted neurons
				for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
					previous_dendritic_spike_times_global[RA_neurons_bursted_global[j]].push_back(internal_time);
				
				////////////////////////////////////////////////////////////////////////////////////
				// LTP or LTD of previously delivered somatic spike times on bursted HVC-RA neuron
				////////////////////////////////////////////////////////////////////////////////////
				
				for (int i = 0; i < N_RA_local; i++)
                {
					int presyn_ID = Id_RA_local[i]; // real id of presynaptic neuron
					
                    // if neuron is saturated apply LTP only if spiked neurons are among supersynapse targets
                    if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
                    {
                        for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_bursted_global[j]; // id of postsynaptic neuron
                            
                            // do not allow self-to-self synapses to emerge
							if ( presyn_ID != postsyn_ID )
							{
								std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
											supersynapses_local[i].end(), postsyn_ID );

								if ( pos!=supersynapses_local[i].end() )
								{
									// find previous somatic spikes that are not too old
									if ( !delivered_spike_times[i][postsyn_ID].empty() )
									{
										std::vector<double>::iterator it_relevant_spikes = std::lower_bound(delivered_spike_times[i][postsyn_ID].begin(),
																											delivered_spike_times[i][postsyn_ID].end(),
																											internal_time - STDP_WINDOW);
																							
										for (auto it = it_relevant_spikes; it != delivered_spike_times[i][postsyn_ID].end(); it++)
										{
											double dt = internal_time - *it;
									
										
											//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
											
											if (dt <= synaptic_params.T_0)
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												LTD(weights_RA_RA_local[i][postsyn_ID], dt);
												
												//~ std::cout   << "LTD from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
											}
											else
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												if ( rescaled_indicators_local[i][postsyn_ID] == 0 )
													LTP(weights_RA_RA_local[i][postsyn_ID], dt);
												else
													LTP_toRescaled(weights_RA_RA_local[i][postsyn_ID], dt);
												//~ std::cout   << "LTP from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
												
											}	
											//double w = weights_local[i][fired_ID];
											
											//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
											 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
											 //           dt, weights_local[i][fired_ID] - w);
											
											update_synapse(i, postsyn_ID);
											
										}
									}
								}
							}
                        }

                    }
                    // if not saturated apply LTP for all dendritic spikes
                    else
                    {
                        for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_bursted_global[j];
                            // don't allow self-to-self connections
                            if ( postsyn_ID != presyn_ID )
                            {
                                // find previous somatic spikes that are not too old
								if ( !delivered_spike_times[i][postsyn_ID].empty() )
								{
									std::vector<double>::iterator it_relevant_spikes = std::lower_bound(delivered_spike_times[i][postsyn_ID].begin(),
																										delivered_spike_times[i][postsyn_ID].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it = it_relevant_spikes; it != delivered_spike_times[i][postsyn_ID].end(); it++)
									{
										double dt = internal_time - *it;
								

										//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
										
									    if (dt <= synaptic_params.T_0)
									    {
											double w_before = weights_RA_RA_local[i][postsyn_ID];
												
											LTD(weights_RA_RA_local[i][postsyn_ID], dt);
											
											// show LTD on pool neurons only
											auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											
											if (p.first == p.second)
											{
												std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											}
											
											//~ std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
											
										}
										else
										{
											double w_before = weights_RA_RA_local[i][postsyn_ID];
											
											if ( rescaled_indicators_local[i][postsyn_ID] == 0 )
												LTP(weights_RA_RA_local[i][postsyn_ID], dt);
											else
												LTP_toRescaled(weights_RA_RA_local[i][postsyn_ID], dt);
										
											
											// show LTP on pool neurons only
											auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											
											if (p.first == p.second)
											{
												std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											}
											
											
											//~ std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
										}	
										//double w = weights_local[i][fired_ID];
										
										//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
										 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
										 //           dt, weights_local[i][fired_ID] - w);
										
										update_synapse(i, postsyn_ID);
									
									}
									
									
								}
                                
                            }
                        }
                   }

                } // end for i -> N_RA_local

				       

               // check if we need axon remodeling

                for (int i = 0; i < N_RA_local; i++)
                {
                    if ( (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss) && (remodeled_local[i] == 0) )
                    {
					    this->axon_remodeling(i);

				    }
                }
            
				// clear bursts
				RA_neurons_bursted_local.clear();
				RA_neurons_bursted_global.clear();
				
            } // end if some HVC-RA bursted
            
            
            network_time += NETWORK_UPDATE_FREQUENCY;
        }
    }
    
    this->potentiation_decay_sudden_maturation();
    //printf("After potentiation decay")
    //this->update_all_synapses_sudden_maturation();
    this->update_all_synapses();
	
	// update maturation info
	std::vector<int> RA_matured_local;
	std::vector<int> RA_matured_global;
	
	for (int i = 0; i < N_RA_local; i++)
	{
		num_spikes_in_recent_trials_local[i].push_front(static_cast<int>(spikes_in_trial_soma_local[i].size()));
        
        //firing_rate_short_local[i] = std::accumulate(num_spikes_in_recent_trials[i].begin(), num_spikes_in_recent_trials[i].begin() + RATE_WINDOW_SHORT, 0.0)
          //                          / static_cast<double>(RATE_WINDOW_SHORT);

		// calculate firing rate in large window:
		int maturation_window = 300;
		double maturation_threshold = 0.7;
		
        firing_rate_long_local[i] = std::accumulate(num_spikes_in_recent_trials_local[i].begin(), num_spikes_in_recent_trials_local[i].end(), 0.0)
                                                                                            / static_cast<double>(maturation_params.RATE_WINDOW_LONG);
        if ( mature_global[Id_RA_local[i]] != 1)
        {
			
			std::vector<double> recent_firing_robustness(maturation_window);
			
			for (int j = 0; j < maturation_window; j++ )
				if ( num_spikes_in_recent_trials_local[i][j] > 0 )
					recent_firing_robustness[j] = 1;
			
			double firing_robustness = std::accumulate(recent_firing_robustness.begin(), recent_firing_robustness.end(), 0.0) / static_cast<double>(maturation_window);
			
			if ( firing_robustness > maturation_threshold )
			{
				RA_matured_local.push_back(Id_RA_local[i]);
					
				this->set_neuron_mature(i);
			}                                                                                                                          
			//if (std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].end(), 0.0) >= 1)
			//{
			//	std::cout << "Recent bursts of neuron " << Id_RA_local[i] << " :" << std::endl;
				
			//	for (size_t j = 0; j < num_bursts_in_recent_trials[i].size(); j++)
			//		std::cout << num_bursts_in_recent_trials[i][j] << "\t";
				
			//	std::cout << std::endl;
			//}
		}
	}
	
	this->gather_spiked_or_bursted_neurons(RA_matured_local, RA_matured_global);
	
	for (size_t i = 0; i < RA_matured_global.size(); i++)
	{
		int neuron_id = RA_matured_global[i];
		
		mature_global[neuron_id] = 1;
		this->rescale_synapses_to_mature(neuron_id);
	}
}


void NetworkGrowthSimulator::trial_burst_pre_dend_post_delays_sudden_maturation_noImmatureOut(bool training)
{
	double event_window = 10.0; // window in which all somatic spikes are considered as one event
	
	// indicators for conductance updates and bursting
	int some_RA_inh_conductance_was_updated_global = 0;
	int some_RA_inh_conductance_was_updated_local = 0;
	
	int some_RA_exc_conductance_was_updated_global = 0;
	int some_RA_exc_conductance_was_updated_local = 0;
	
	int some_I_exc_conductance_was_updated_global = 0;
	int some_I_exc_conductance_was_updated_local = 0;
	
	int some_RA_neuron_bursted_global = 0;
	int some_RA_neuron_bursted_local = 0;
	
	 // arrays with updates for conductances
    std::vector<double> update_Ge_RA_local(N_RA);
	std::vector<double> update_Ge_RA_global(N_RA);

	std::vector<double> update_Gi_RA_local(N_RA);
	std::vector<double> update_Gi_RA_global(N_RA);

	std::vector<double> update_Ge_I_local(N_I);
	std::vector<double> update_Ge_I_global(N_I);
	
	std::vector<int> RA_neurons_bursted_local; // local array of HVC-RA neurons that bursted
	std::vector<int> RA_neurons_bursted_global; // global array of HVC-RA neurons that bursted
	
	// sample training innervation time and send to all processes
    double training_kick_time;
    
	if (MPI_rank == 0)
	{
		training_kick_time = WAITING_TIME + noise_generator.random(TRIAL_DURATION - 2*WAITING_TIME);
		std::cout << "training_kick_time = " << training_kick_time << std::endl;
    }
    MPI_Bcast(&training_kick_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    bool training_excited = false; // indicator that training neurons were already excited
    
    double internal_time = 0;
    double network_time = NETWORK_UPDATE_FREQUENCY;

   

    //printf("network time = %f\n", network_time);

    
	// initialize update arrays and fired indicators
	std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
	
	
    // evolve dynamics
    int num_steps = static_cast<int>(TRIAL_DURATION / TIMESTEP);
    
    for (int t = 1; t < num_steps; t++)
	{
		internal_time += TIMESTEP;
	
		if ( ( training ) && ( !training_excited ) && ( internal_time >= training_kick_time ) )
		{
			for (int i = 0; i < N_TR; i++)
			{
				int rank;
				int shift;
				
				this->get_neuronRA_location(training_neurons[i], &rank, &shift);
				
				if (MPI_rank == rank)
					HVCRA_local[shift].raiseE(G_TRAINING_KICK);
			}
			
			training_excited = true;
		}
		
		//////////////////////
		// HVC-RA neurons
		//////////////////////
		
		for (int i = 0; i < N_RA_local; i++)
		{
            // Debraband step
            HVCRA_local[i].Debraband_step_no_target_update();
            
           
			///////////////////////////////
			// HVC-RA somatic spike
			///////////////////////////////
			// if some neuron produced somatic spike
            if (HVCRA_local[i].get_fired_soma())
            {
                spikes_in_trial_soma_local[i].push_back(internal_time);
                
				// update delivery queues if neuron is mature (no output for immature neuron and no wiring mechanism)
				if ( mature_global[Id_RA_local[i]] == 1 ){
				
					// for inhibitory neurons
					// loop over all inhibitory targets of fired neurons
					size_t num_I_targets = syn_ID_RA_I_local[i].size();
					
					for (size_t j = 0; j < num_I_targets; j++)
					{
						double delivery_time = internal_time + axonal_delays_RA_I_local[i][j];
						
						// if queue is empty, just add item to the queue
						if ( delivery_queue_RA_I[i].empty() )
							delivery_queue_RA_I[i].push_back(std::pair<double,int>(delivery_time, j));	
						// otherwise add item so that queue is sorted
						else
						{
							auto it = std::upper_bound(delivery_queue_RA_I[i].begin(), delivery_queue_RA_I[i].end(), std::pair<double,int>(delivery_time, j));
							delivery_queue_RA_I[i].insert(it, std::pair<double,int>(delivery_time, j));
						}
					}
						//std::cout << "neuron " << Id_RA_local[i] << " spike time = " << internal_time << " axonal delay = " << axonal_delays_RA_I[Id_RA_local[i]][syn_ID] << " delivery time RA to I: " << delivery_queue_RA_I[i].back().first << " delivery target id: " << syn_ID_RA_I_local[i][delivery_queue_RA_I[i].back().second] << std::endl;
					
					
					// if neuron is saturated, deliver spikes only to supersynaptic targets
					if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
					{
						for (size_t j = 0; j < supersynapses_local[i].size(); j++)
						{
							int target_id = supersynapses_local[i][j];
							
							double delivery_time = internal_time + axonal_delays_RA_RA_local[i][target_id];
							
							// if queue is empty, just add item to the queue
							if ( delivery_queue_RA_RA_soma[i].empty() )
								delivery_queue_RA_RA_soma[i].push_back(std::pair<double,int>(delivery_time, target_id));	
							// otherwise add item so that queue is sorted
							else
							{
								auto it = std::upper_bound(delivery_queue_RA_RA_soma[i].begin(), delivery_queue_RA_RA_soma[i].end(), std::pair<double,int>(delivery_time, target_id));
								delivery_queue_RA_RA_soma[i].insert(it, std::pair<double,int>(delivery_time, target_id));
							}
						}
					}
					else // deliver spikes to everyone except itself
					{	
						for (int j = 0; j < N_RA; j++)
						{
							if ( j != Id_RA_local[i])
							{
							
								double delivery_time = internal_time + axonal_delays_RA_RA_local[i][j];
							
								// if queue is empty, just add item to the queue
								if ( delivery_queue_RA_RA_soma[i].empty() )
									delivery_queue_RA_RA_soma[i].push_back(std::pair<double,int>(delivery_time, j));	
								// otherwise add item so that queue is sorted
								else
								{
									auto it = std::upper_bound(delivery_queue_RA_RA_soma[i].begin(), delivery_queue_RA_RA_soma[i].end(), std::pair<double,int>(delivery_time, j));
									delivery_queue_RA_RA_soma[i].insert(it, std::pair<double,int>(delivery_time, j));
								}
							}
						}
					}
				}
                
            } // end if local HVC-RA neuron spiked  
            
             //////////////////////////////////
			// HVC-RA -> HVC-I delivery queue
			//////////////////////////////////
            // check if somatic spike was delivered to some interneuron
            // loop through the delivery queue to check if current time exceeds the spike delivery time
            std::vector<std::pair<double,int>>::iterator it = delivery_queue_RA_I[i].begin();
			
            for (; it != delivery_queue_RA_I[i].end(); it++)
            {
				if (internal_time >= it->first)
				{
					
					some_I_exc_conductance_was_updated_local = 1;
					
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_RA_I_local[i][pos_in_local_target_array];
					
					update_Ge_I_local[target_id] += weights_RA_I_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(RA) neuron " << Id_RA_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to HVC(I) neuron " << syn_ID_RA_I_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_I[i].erase(delivery_queue_RA_I[i].begin(), it);
			
			//////////////////////////////////////////////////
			// HVC-RA -> HVC-RA somatic spike delivery queue
			//////////////////////////////////////////////////
			// check if somatic spike was delivered to some HVC(RA) neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
            it = delivery_queue_RA_RA_soma[i].begin();
            
            for (; it != delivery_queue_RA_RA_soma[i].end(); it++)
            {
				double delivery_time = it->first;
				
				if ( internal_time >= delivery_time )
				{
					int target_id = it->second;
					
					
					
					// update conductance of target if synapse is active
					if ( active_indicators_local[i][target_id] == 1 )
					{
						some_RA_exc_conductance_was_updated_local = 1;
						update_Ge_RA_local[target_id] += weights_RA_RA_local[i][target_id];
					}
					
					
					
					//////////////
					//// LTD
					//////////////
					// do LTD on dendritic spike time of target neuron if it is not mature neuron
					//~ if ( mature_global[target_id] != 1 )
					//~ {
						// update delivered spikes only if current spike occured later than previous delivered
						// spike + event_window
						bool new_event_occured = false;
					
						if ( !delivered_spike_times[i][target_id].empty() ){
							if ( delivery_time > delivered_spike_times[i][target_id].back() + event_window ){
								delivered_spike_times[i][target_id].push_back(delivery_time);
								new_event_occured = true;
							}
						}
						else{
							delivered_spike_times[i][target_id].push_back(delivery_time);
							new_event_occured = true;
						}
						
						if ( new_event_occured ){
							 // if neuron is saturated apply LTD only if target is among super synapses
							if (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss)
							{
								if ( supersynapses_indicators_local[i][target_id] == 1 )
								{	
									// find previous dendritic spikes that are not too old
									if ( !previous_dendritic_spike_times_global[target_id].empty() )
									{
										std::vector<double>::iterator it_relevant_dend_spikes = std::lower_bound(previous_dendritic_spike_times_global[target_id].begin(),
																											previous_dendritic_spike_times_global[target_id].end(),
																											internal_time - STDP_WINDOW);
																							
										for (auto it_dend_spike = it_relevant_dend_spikes; it_dend_spike != previous_dendritic_spike_times_global[target_id].end(); it_dend_spike++)
										{
											double dt = *it_dend_spike - delivery_time;

											//std::cout << "dt in saturated LTD = " << dt << std::endl;

											
											double w_before = weights_RA_RA_local[i][target_id];     
								   
											LTD(weights_RA_RA_local[i][target_id], dt);
											
											
											//~ std::cout << "LTD from saturated " << Id_RA_local[i] << " -> " << supersynapse_id
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][supersynapse_id] - w_before
													  //~ << std::endl;
													  
											//printf("LTD from saturated neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n",
											//            Id_RA_local[i], supersynapse_id, spikes_in_trial_soma_local[i].back(), spike_times_dend_global[supersynapse_id],
											//            dt, weights_local[i][supersynapse_id] - w);
											
											update_synapse(i, target_id);	
											
										}   
									}
								}

								// if some supersynapse desaturated, update all synapses
								if (static_cast<int>(supersynapses_local[i].size()) < synaptic_params.Nss)
									for (int j = 0; j < N_RA; j++)
										this->update_synapse(i, j);
							}
							// if not saturated apply LTD rule 
							else
							{
								// find previous dendritic spikes that are not too old
								if ( !previous_dendritic_spike_times_global[target_id].empty() )
								{
									std::vector<double>::iterator it_relevant_dend_spikes = std::lower_bound(previous_dendritic_spike_times_global[target_id].begin(),
																										previous_dendritic_spike_times_global[target_id].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it_dend_spike = it_relevant_dend_spikes; it_dend_spike != previous_dendritic_spike_times_global[target_id].end(); it_dend_spike++)
									{
										double dt = *it_dend_spike - delivery_time;

										//std::cout << "dt in saturated LTD = " << dt << std::endl;

										
										double w_before = weights_RA_RA_local[i][target_id];     
							   
										LTD(weights_RA_RA_local[i][target_id], dt);
										
										//~ // show LTD to pool neurons only
										auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), target_id);
										
										if (p.first == p.second)
										{
											std::cout << "LTD from  " << Id_RA_local[i] << " -> " << target_id
												  << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][target_id] - w_before
												  << " w_after = " << weights_RA_RA_local[i][target_id] << std::endl;
										}
									//~ 
										update_synapse(i, target_id);	
										
									}   
								}
							}
						}
					//~ } // end if target neuron is not mature
				} // end if internal_time >= delivery spike time 
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_RA_soma[i].erase(delivery_queue_RA_RA_soma[i].begin(), it);
			
			
			///////////////////////////////
			// HVC-RA dendritic spike
			///////////////////////////////
			//if some neuron produced dendritic spike, store this neuron in array
			
			if ( HVCRA_local[i].get_fired_dend() )
			{
				// show spikes of pool neurons only
				auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
				
				if (p.first == p.second)
					std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
				
				
				spikes_in_trial_dend_local[i].push_back(internal_time); // dend spike time relative to trial onset
				
				//if ( mature_global[Id_RA_local[i]] != 1 )
				//{
					RA_neurons_bursted_local.push_back(Id_RA_local[i]);
				
					some_RA_neuron_bursted_local = 1;
				//}
				//std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
			} // end if get_fired_dend
	
            
		} // end for i = 0 -> N_RA_local (loop through all local HVC-RA neurons)
		
		//////////////////////
		// HVC-I neurons
		//////////////////////
		for (int i = 0; i < N_I_local; i++)
		{
			HVCI_local[i].DP8_step_no_target_update();
			
			
			///////////////////////////////
			// HVC-I spike
			///////////////////////////////
			//  if some I neuron spikes, update delivery queue
			if (HVCI_local[i].get_fired())
			{
				//printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
				spikes_in_trial_interneuron_local[i].push_back(internal_time);

				size_t num_RA_targets = syn_ID_I_RA_local[i].size();
				// loop over all targets of fired neurons
				for (size_t j = 0; j < num_RA_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_I_RA_local[i][j];
					
					//std::cout << "HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << syn_ID_I_RA_local[i][j] << " spike time " << internal_time << " delivery time = " << delivery_time << std::endl; 
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_I_RA[i].empty() )
						delivery_queue_I_RA[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_I_RA[i].begin(), delivery_queue_I_RA[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_I_RA[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
			}
			
			//////////////////////////////////
			// HVC-I -> HVC-RA delivery queue
			//////////////////////////////////
			// check if interneuron spike was delivered to some HVC-RA neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
			std::vector<std::pair<double,int>>::iterator it = delivery_queue_I_RA[i].begin();
			
			for (; it != delivery_queue_I_RA[i].end(); it++)
			{
				if (internal_time >= it->first)
				{
					some_RA_inh_conductance_was_updated_local = 1;
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_I_RA_local[i][pos_in_local_target_array];
					
					//std::cout << "Delivered spike HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << target_id << " delivered_time " << internal_time << " delivery time in queue = " << it->first << std::endl; 
					
					
					update_Gi_RA_local[target_id] += weights_I_RA_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(I) neuron " << Id_I_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to neuron " << syn_ID_I_RA_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_I_RA[i].erase(delivery_queue_I_RA[i].begin(), it);
		
		}  // end for i = 0 -> N_I_local (loop through all local HVC-I neurons)
		
		/////////////////////////////
		// Network Synchronization //
		/////////////////////////////
		if ( ( internal_time > network_time ) || ( t == num_steps-1 ) )
        {
            //if (MPI_rank == 0)
            //{
            //    std::cout << "internal time = " << internal_time << " > network_time = " << network_time << std::endl; 
            //}

            MPI_Allreduce(&some_RA_exc_conductance_was_updated_local, &some_RA_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_RA_inh_conductance_was_updated_local, &some_RA_inh_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_I_exc_conductance_was_updated_local, &some_I_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_RA_neuron_bursted_local, &some_RA_neuron_bursted_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            
        
			
            ////////////////////////////////////////////
            //// Process HVC-I -> HVC-RA interactions
            ////////////////////////////////////////////
            if ( some_RA_inh_conductance_was_updated_global > 0 )
            {
            // sum update array and send to all processes

                MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                for (int i = 0; i < N_RA_local; i++)
                {
					//std::cout << "HVC-RA neuron " << Id_RA_local[i] << " raised inhibitory conductance by " << update_Gi_RA_global[Id_RA_local[i]] << std::endl;
                    HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
            	}
				// update conductance arrays and delivered indicators
				some_RA_inh_conductance_was_updated_global = 0;
				some_RA_inh_conductance_was_updated_local = 0;
				
				std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            }
            
            ////////////////////////////////////////////
            //// Process HVC-RA -> HVC-I interactions
            ////////////////////////////////////////////
			if ( some_I_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons

				for (int i = 0; i < N_I_local; i++)
					HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);

				// update conductance arrays and fired indicators
				some_I_exc_conductance_was_updated_local = 0;
				some_I_exc_conductance_was_updated_global = 0;

				std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
			}
			
			////////////////////////////////////////////
            //// Process HVC-RA -> HVC-RA interactions
            ////////////////////////////////////////////
			if ( some_RA_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons
				for (int i = 0; i < N_RA_local; i++)
					HVCRA_local[i].raiseExcWeight(update_Ge_RA_global[Id_RA_local[i]]); // update conductance

				// update conductance arrays and fired indicators
				some_RA_exc_conductance_was_updated_global = 0;
				some_RA_exc_conductance_was_updated_local = 0;
				
				std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
			}	

			////////////////////////////////////////////
            //// Process HVC-RA dendritic spikes
            ////////////////////////////////////////////
			if ( some_RA_neuron_bursted_global > 0 )
            {
				// update conductance arrays and fired indicators
            	some_RA_neuron_bursted_local = 0;
	        	some_RA_neuron_bursted_global = 0;
	        	
	        	// gather all bursted HVC-RA neurons
				this->gather_spiked_or_bursted_neurons(RA_neurons_bursted_local, RA_neurons_bursted_global);
				
				// add dendritic spikes for bursted neurons
				for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
					previous_dendritic_spike_times_global[RA_neurons_bursted_global[j]].push_back(internal_time);
				
				////////////////////////////////////////////////////////////////////////////////////
				// LTP or LTD of previously delivered somatic spike times on bursted HVC-RA neuron
				////////////////////////////////////////////////////////////////////////////////////
				
				for (int i = 0; i < N_RA_local; i++)
                {
					int presyn_ID = Id_RA_local[i]; // real id of presynaptic neuron
					
                    // if neuron is saturated apply LTP only if spiked neurons are among supersynapse targets
                    if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
                    {
                        for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_bursted_global[j]; // id of postsynaptic neuron
                            
                            // do not allow self-to-self synapses to emerge
							if ( presyn_ID != postsyn_ID )
							{
								std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
											supersynapses_local[i].end(), postsyn_ID );

								if ( pos!=supersynapses_local[i].end() )
								{
									// find previous somatic spikes that are not too old
									if ( !delivered_spike_times[i][postsyn_ID].empty() )
									{
										std::vector<double>::iterator it_relevant_spikes = std::lower_bound(delivered_spike_times[i][postsyn_ID].begin(),
																											delivered_spike_times[i][postsyn_ID].end(),
																											internal_time - STDP_WINDOW);
																							
										for (auto it = it_relevant_spikes; it != delivered_spike_times[i][postsyn_ID].end(); it++)
										{
											double dt = internal_time - *it;
									
										
											//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
											
											if (dt <= synaptic_params.T_0)
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												LTD(weights_RA_RA_local[i][postsyn_ID], dt);
												
												//~ std::cout   << "LTD from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
											}
											else
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												if ( rescaled_indicators_local[i][postsyn_ID] == 0 )
													LTP(weights_RA_RA_local[i][postsyn_ID], dt);
												else
													LTP_toRescaled(weights_RA_RA_local[i][postsyn_ID], dt);
												//~ std::cout   << "LTP from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
												
											}	
											//double w = weights_local[i][fired_ID];
											
											//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
											 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
											 //           dt, weights_local[i][fired_ID] - w);
											
											update_synapse(i, postsyn_ID);
											
										}
									}
								}
							}
                        }

                    }
                    // if not saturated apply LTP for all dendritic spikes
                    else
                    {
                        for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_bursted_global[j];
                            // don't allow self-to-self connections
                            if ( postsyn_ID != presyn_ID )
                            {
                                // find previous somatic spikes that are not too old
								if ( !delivered_spike_times[i][postsyn_ID].empty() )
								{
									std::vector<double>::iterator it_relevant_spikes = std::lower_bound(delivered_spike_times[i][postsyn_ID].begin(),
																										delivered_spike_times[i][postsyn_ID].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it = it_relevant_spikes; it != delivered_spike_times[i][postsyn_ID].end(); it++)
									{
										double dt = internal_time - *it;
								

										//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
										
									    if (dt <= synaptic_params.T_0)
									    {
											double w_before = weights_RA_RA_local[i][postsyn_ID];
												
											LTD(weights_RA_RA_local[i][postsyn_ID], dt);
											
											// show LTD on pool neurons only
											auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											
											if (p.first == p.second)
											{
												std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											}
											
											//~ std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
											
										}
										else
										{
											double w_before = weights_RA_RA_local[i][postsyn_ID];
											
											if ( rescaled_indicators_local[i][postsyn_ID] == 0 )
												LTP(weights_RA_RA_local[i][postsyn_ID], dt);
											else
												LTP_toRescaled(weights_RA_RA_local[i][postsyn_ID], dt);
										
											
											// show LTP on pool neurons only
											auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											
											if (p.first == p.second)
											{
												std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											}
											
											
											//~ std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
										}	
										//double w = weights_local[i][fired_ID];
										
										//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
										 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
										 //           dt, weights_local[i][fired_ID] - w);
										
										update_synapse(i, postsyn_ID);
									
									}
									
									
								}
                                
                            }
                        }
                   }

                } // end for i -> N_RA_local

				       

               // check if we need axon remodeling

                for (int i = 0; i < N_RA_local; i++)
                {
                    if ( (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss) && (remodeled_local[i] == 0) )
                    {
					    this->axon_remodeling(i);

				    }
                }
            
				// clear bursts
				RA_neurons_bursted_local.clear();
				RA_neurons_bursted_global.clear();
				
            } // end if some HVC-RA bursted
            
            
            network_time += NETWORK_UPDATE_FREQUENCY;
        }
    }
    
    this->potentiation_decay_sudden_maturation();
    //printf("After potentiation decay")
    //this->update_all_synapses_sudden_maturation();
    this->update_all_synapses();
	
	// update maturation info
	std::vector<int> RA_matured_local;
	std::vector<int> RA_matured_global;
	
	for (int i = 0; i < N_RA_local; i++)
	{
		num_spikes_in_recent_trials_local[i].push_front(static_cast<int>(spikes_in_trial_soma_local[i].size()));
        
        //firing_rate_short_local[i] = std::accumulate(num_spikes_in_recent_trials[i].begin(), num_spikes_in_recent_trials[i].begin() + RATE_WINDOW_SHORT, 0.0)
          //                          / static_cast<double>(RATE_WINDOW_SHORT);

		// calculate firing rate in large window:
		int maturation_window = 300;
		double maturation_threshold = 0.7;
		
        firing_rate_long_local[i] = std::accumulate(num_spikes_in_recent_trials_local[i].begin(), num_spikes_in_recent_trials_local[i].end(), 0.0)
                                                                                            / static_cast<double>(maturation_params.RATE_WINDOW_LONG);
        if ( mature_global[Id_RA_local[i]] != 1)
        {
			
			std::vector<double> recent_firing_robustness(maturation_window);
			
			for (int j = 0; j < maturation_window; j++ )
				if ( num_spikes_in_recent_trials_local[i][j] > 0 )
					recent_firing_robustness[j] = 1;
			
			double firing_robustness = std::accumulate(recent_firing_robustness.begin(), recent_firing_robustness.end(), 0.0) / static_cast<double>(maturation_window);
			
			if ( firing_robustness > maturation_threshold )
			{
				RA_matured_local.push_back(Id_RA_local[i]);
					
				this->set_neuron_mature(i);
			}                                                                                                                          
			//if (std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].end(), 0.0) >= 1)
			//{
			//	std::cout << "Recent bursts of neuron " << Id_RA_local[i] << " :" << std::endl;
				
			//	for (size_t j = 0; j < num_bursts_in_recent_trials[i].size(); j++)
			//		std::cout << num_bursts_in_recent_trials[i][j] << "\t";
				
			//	std::cout << std::endl;
			//}
		}
	}
	
	this->gather_spiked_or_bursted_neurons(RA_matured_local, RA_matured_global);
	
	for (size_t i = 0; i < RA_matured_global.size(); i++)
	{
		int neuron_id = RA_matured_global[i];
		
		mature_global[neuron_id] = 1;
		this->rescale_synapses_to_mature(neuron_id);
	}
}

void NetworkGrowthSimulator::trial_1stSoma_pre_1stSoma_post_delays_fixedSpread_with_inhibition_tracking(bool training, std::vector<double>& spread_times, 
																					double time_resolution_conductance)
{
	double event_window = 30.0; // window in which all somatic spikes are considered as one event
	
	// indicators for conductance updates and bursting
	int some_RA_inh_conductance_was_updated_global = 0;
	int some_RA_inh_conductance_was_updated_local = 0;
	
	int some_RA_exc_conductance_was_updated_global = 0;
	int some_RA_exc_conductance_was_updated_local = 0;
	
	int some_I_exc_conductance_was_updated_global = 0;
	int some_I_exc_conductance_was_updated_local = 0;
	
	int some_RA_neuron_new_postEvent_global = 0;
	int some_RA_neuron_new_postEvent_local = 0;
	
	 // arrays with updates for conductances
    std::vector<double> update_Ge_RA_local(N_RA);
	std::vector<double> update_Ge_RA_global(N_RA);

	std::vector<double> update_Gi_RA_local(N_RA);
	std::vector<double> update_Gi_RA_global(N_RA);

	std::vector<double> update_Ge_I_local(N_I);
	std::vector<double> update_Ge_I_global(N_I);
	
	std::vector<int> RA_neurons_new_postEvent_local; // local array of HVC-RA neurons that produced new postsynaptic event
	std::vector<int> RA_neurons_new_postEvent_global; // global array of HVC-RA neurons that produced new postsynaptic event
	
	std::vector<std::vector<std::vector<double>>> delivered_preEvent_times_local(N_RA_local); // local array with the times when presynaptic events 
															 // source HVC-RA neuron were delivered to other HVC-RA neurons
		
	for (int i = 0; i < N_RA_local; i++)
		delivered_preEvent_times_local[i].resize(N_RA);
		
		
	std::vector<std::vector<std::pair<double,int>>> delivery_queue_RA_RA_soma_local(N_RA_local); // local queue with spike delivery times for HVC-RA -> HVC-RA interactions
	std::vector<std::vector<std::pair<double,int>>> delivery_queue_RA_I_local(N_RA_local); // local queue with spike delivery times for HVC-RA -> HVC-I interactions
	std::vector<std::vector<std::pair<double,int>>> delivery_queue_I_RA_local(N_I_local); // local queue with spike delivery times for HVC-I -> HVC-RA interactions
		
	std::vector<std::vector<double>> previous_postEvent_times_global(N_RA); // array with the previous somatic spike times of neurons that are not too old
		
	// sample training innervation time and send to all processes
  // training neurons innervation
	std::vector<bool> indicators_training(N_RA);
	std::fill(indicators_training.begin(), indicators_training.end(), false);
	
	for (int i = 0; i < N_TR; i++)
		indicators_training[training_neurons[i]] = true;
		
	std::vector<bool> indicators_current_injected(N_RA);
	std::fill(indicators_current_injected.begin(), indicators_current_injected.end(), false);
	
	// sample mean training innervation time and send to all processes
    double mean_injection_time;
    
	if (MPI_rank == 0)
	{
		mean_injection_time = WAITING_TIME + noise_generator.random(TRIAL_DURATION - 2*WAITING_TIME);
		std::cout << "mean_injection_time = " << mean_injection_time << std::endl;
    }
    MPI_Bcast(&mean_injection_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
    double internal_time = 0;
    double network_time = NETWORK_UPDATE_FREQUENCY;

   

    //printf("network time = %f\n", network_time);

    
	// initialize update arrays and fired indicators
	std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
	
	
    // evolve dynamics
    int num_steps = static_cast<int>(TRIAL_DURATION / TIMESTEP);
    
    for (int i = 0; i < N_RA_local; i++)
		Ginh_local[i][0] = HVCRA_local[i].get_Ginh_d();
    
    int position_in_local_inhibition_array = 1;
    
    for (int t = 1; t < num_steps; t++)
	{
		internal_time += TIMESTEP;
		
		//////////////////////
		// HVC-RA neurons
		//////////////////////
		
		for (int i = 0; i < N_RA_local; i++)
		{
			if ( (training) && ( indicators_training[Id_RA_local[i]] ) && ( !indicators_current_injected[Id_RA_local[i]] ) )
			{
				if ( internal_time > spread_times[Id_RA_local[i]] + mean_injection_time )
				{
					HVCRA_local[i].raiseE(G_TRAINING_KICK);
					indicators_current_injected[Id_RA_local[i]] = true;
				}
			}
			
            // Debraband step
            HVCRA_local[i].Debraband_step_no_target_update();
            
			///////////////////////////////
			// HVC-RA somatic spike
			///////////////////////////////
			// if some neuron produced somatic spike
            if (HVCRA_local[i].get_fired_soma())
            {
                spikes_in_trial_soma_local[i].push_back(internal_time);

				// update somatic spikes only if current spike occured later than previous 
				// somatic spike + event_window
				if ( previous_postEvent_times_global[Id_RA_local[i]].empty() )
				{
					RA_neurons_new_postEvent_local.push_back(Id_RA_local[i]);
					some_RA_neuron_new_postEvent_local = 1;
					//std::cout << "HVC-RA " << Id_RA_local[i] << " postsynaptic event at " << internal_time << std::endl;
				}
				else if ( internal_time > previous_postEvent_times_global[Id_RA_local[i]].back() + event_window)
				{
					RA_neurons_new_postEvent_local.push_back(Id_RA_local[i]);
					some_RA_neuron_new_postEvent_local = 1;
					//std::cout << "HVC-RA " << Id_RA_local[i] << " postsynaptic event at " << internal_time << std::endl;
				}
							
				// for inhibitory neurons
				// loop over all inhibitory targets of fired neurons
				size_t num_I_targets = syn_ID_RA_I_local[i].size();
				
				for (size_t j = 0; j < num_I_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_RA_I_local[i][j];
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_RA_I_local[i].empty() )
						delivery_queue_RA_I_local[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_RA_I_local[i].begin(), delivery_queue_RA_I_local[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_RA_I_local[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
					//std::cout << "neuron " << Id_RA_local[i] << " spike time = " << internal_time << " axonal delay = " << axonal_delays_RA_I[Id_RA_local[i]][syn_ID] << " delivery time RA to I: " << delivery_queue_RA_I[i].back().first << " delivery target id: " << syn_ID_RA_I_local[i][delivery_queue_RA_I[i].back().second] << std::endl;
				
				
				// if neuron is saturated, deliver spikes only to supersynaptic targets
				if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
				{
					for (size_t j = 0; j < supersynapses_local[i].size(); j++)
					{
						int target_id = supersynapses_local[i][j];
						
						double delivery_time = internal_time + axonal_delays_RA_RA_local[i][target_id];
						
						// if queue is empty, just add item to the queue
						if ( delivery_queue_RA_RA_soma_local[i].empty() )
							delivery_queue_RA_RA_soma_local[i].push_back(std::pair<double,int>(delivery_time, target_id));	
						// otherwise add item so that queue is sorted
						else
						{
							auto it = std::upper_bound(delivery_queue_RA_RA_soma_local[i].begin(), delivery_queue_RA_RA_soma_local[i].end(), std::pair<double,int>(delivery_time, target_id));
							delivery_queue_RA_RA_soma_local[i].insert(it, std::pair<double,int>(delivery_time, target_id));
						}
					}
				}
				else // deliver spikes to everyone except itself
				{	
					for (int j = 0; j < N_RA; j++)
					{
						if ( j != Id_RA_local[i])
						{
						
							double delivery_time = internal_time + axonal_delays_RA_RA_local[i][j];
						
							// if queue is empty, just add item to the queue
							if ( delivery_queue_RA_RA_soma_local[i].empty() )
								delivery_queue_RA_RA_soma_local[i].push_back(std::pair<double,int>(delivery_time, j));	
							// otherwise add item so that queue is sorted
							else
							{
								auto it = std::upper_bound(delivery_queue_RA_RA_soma_local[i].begin(), delivery_queue_RA_RA_soma_local[i].end(), std::pair<double,int>(delivery_time, j));
								delivery_queue_RA_RA_soma_local[i].insert(it, std::pair<double,int>(delivery_time, j));
							}
						}
					}
				}
			
            } // end if local HVC-RA neuron spiked  
            
             //////////////////////////////////
			// HVC-RA -> HVC-I delivery queue
			//////////////////////////////////
            // check if somatic spike was delivered to some interneuron
            // loop through the delivery queue to check if current time exceeds the spike delivery time
            std::vector<std::pair<double,int>>::iterator it = delivery_queue_RA_I_local[i].begin();
			
            for (; it != delivery_queue_RA_I_local[i].end(); it++)
            {
				if (internal_time >= it->first)
				{
					
					some_I_exc_conductance_was_updated_local = 1;
					
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_RA_I_local[i][pos_in_local_target_array];
					
					update_Ge_I_local[target_id] += weights_RA_I_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(RA) neuron " << Id_RA_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to HVC(I) neuron " << syn_ID_RA_I_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_I_local[i].erase(delivery_queue_RA_I_local[i].begin(), it);
			
			//////////////////////////////////////////////////
			// HVC-RA -> HVC-RA somatic spike delivery queue
			//////////////////////////////////////////////////
			// check if somatic spike was delivered to some HVC(RA) neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
            it = delivery_queue_RA_RA_soma_local[i].begin();
            
            for (; it != delivery_queue_RA_RA_soma_local[i].end(); it++)
            {
				double delivery_time = it->first;
				
				if ( internal_time >= delivery_time )
				{
					int target_id = it->second;
					
					// update conductance of target if synapse is active
					if ( active_indicators_local[i][target_id] == 1 )
					{
						some_RA_exc_conductance_was_updated_local = 1;
						update_Ge_RA_local[target_id] += weights_RA_RA_local[i][target_id];
					}
					
					//////////////
					//// LTD
					//////////////
					// update delivered spikes only if current spike occured later than previous delivered
					// spike + event_window
					bool new_preEvent_occured = false;
				
					if ( !delivered_preEvent_times_local[i][target_id].empty() ){
						if ( delivery_time > delivered_preEvent_times_local[i][target_id].back() + event_window ){
							delivered_preEvent_times_local[i][target_id].push_back(delivery_time);
							new_preEvent_occured = true;
						}
					}
					else{
						delivered_preEvent_times_local[i][target_id].push_back(delivery_time);
						new_preEvent_occured = true;
					}
					
					if ( new_preEvent_occured ){
						 // if neuron is saturated apply LTD only if target is among super synapses
						if (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss)
						{
							if ( supersynapses_indicators_local[i][target_id] == 1 )
							{	
								// find previous dendritic spikes that are not too old
								if ( !previous_postEvent_times_global[target_id].empty() )
								{
									std::vector<double>::iterator it_relevant_postEvent = std::lower_bound(previous_postEvent_times_global[target_id].begin(),
																										previous_postEvent_times_global[target_id].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it_postEvent = it_relevant_postEvent; it_postEvent != previous_postEvent_times_global[target_id].end(); it_postEvent++)
									{
										double dt = *it_postEvent - delivery_time;

										//std::cout << "dt in saturated LTD = " << dt << std::endl;

										
										double w_before = weights_RA_RA_local[i][target_id];     
							   
										LTD(weights_RA_RA_local[i][target_id], dt);
										
										
										//~ std::cout << "LTD from saturated " << Id_RA_local[i] << " -> " << supersynapse_id
												  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][supersynapse_id] - w_before
												  //~ << std::endl;
												  
										//printf("LTD from saturated neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n",
										//            Id_RA_local[i], supersynapse_id, spikes_in_trial_soma_local[i].back(), spike_times_dend_global[supersynapse_id],
										//            dt, weights_local[i][supersynapse_id] - w);
										
										update_synapse(i, target_id);	
										
									}   
								}
							}

							// if some supersynapse desaturated, update all synapses
							if (static_cast<int>(supersynapses_local[i].size()) < synaptic_params.Nss)
								for (int j = 0; j < N_RA; j++)
									this->update_synapse(i, j);
						}
						// if not saturated apply LTD rule 
						else
						{
							// find previous dendritic spikes that are not too old
							if ( !previous_postEvent_times_global[target_id].empty() )
							{
								std::vector<double>::iterator it_relevant_postEvent = std::lower_bound(previous_postEvent_times_global[target_id].begin(),
																									previous_postEvent_times_global[target_id].end(),
																									internal_time - STDP_WINDOW);
																					
								for (auto it_postEvent = it_relevant_postEvent; it_postEvent != previous_postEvent_times_global[target_id].end(); it_postEvent++)
								{
									double dt = *it_postEvent - delivery_time;

									//std::cout << "dt in saturated LTD = " << dt << std::endl;

									
									double w_before = weights_RA_RA_local[i][target_id];     
						   
									LTD(weights_RA_RA_local[i][target_id], dt);
									
									//~ // show LTD to pool neurons only
									//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), target_id);
									//~ 
									//~ if (p.first == p.second)
									//~ {
										//~ std::cout << "LTD from  " << Id_RA_local[i] << " -> " << target_id
											  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][target_id] - w_before
											  //~ << " w_after = " << weights_RA_RA_local[i][target_id] << std::endl;
									//~ }
								//~ 
									update_synapse(i, target_id);	
									
								}   
							}
						}
					}
				} // end if internal_time >= delivery spike time 
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_RA_soma_local[i].erase(delivery_queue_RA_RA_soma_local[i].begin(), it);
			
			
			///////////////////////////////
			// HVC-RA dendritic spike
			///////////////////////////////
			//if some neuron produced dendritic spike, store this neuron in array
			
			if ( HVCRA_local[i].get_fired_dend() )
			{
				//~ // show spikes of pool neurons only
				//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
				//~ 
				//~ if (p.first == p.second)
					//~ std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
				//~ 
				//std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
				
				spikes_in_trial_dend_local[i].push_back(internal_time); // dend spike time relative to trial onset
				
				
				//std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
			} // end if get_fired_dend
	
            
		} // end for i = 0 -> N_RA_local (loop through all local HVC-RA neurons)
		
		
		if ( internal_time >= static_cast<double>(position_in_local_inhibition_array) * time_resolution_conductance)
		{
			if (position_in_local_inhibition_array >= Ginh_local[0].size())
				std::cout << "internal time = " << internal_time
						  << " position_in_local_inhibition_array = " << position_in_local_inhibition_array
						  << " exceeds size of local conductance array = " << Ginh_local[0].size() << "\n" << std::endl;
			else
				for (int i = 0; i < N_RA_local; i++)
					Ginh_local[i][position_in_local_inhibition_array] = HVCRA_local[i].get_Ginh_d();
			
			position_in_local_inhibition_array += 1;
		}
		
		//////////////////////
		// HVC-I neurons
		//////////////////////
		for (int i = 0; i < N_I_local; i++)
		{
			HVCI_local[i].DP8_step_no_target_update();
			
			
			///////////////////////////////
			// HVC-I spike
			///////////////////////////////
			//  if some I neuron spikes, update delivery queue
			if (HVCI_local[i].get_fired())
			{
				//printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
				spikes_in_trial_interneuron_local[i].push_back(internal_time);

				size_t num_RA_targets = syn_ID_I_RA_local[i].size();
				// loop over all targets of fired neurons
				for (size_t j = 0; j < num_RA_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_I_RA_local[i][j];
					
					//std::cout << "HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << syn_ID_I_RA_local[i][j] << " spike time " << internal_time << " delivery time = " << delivery_time << std::endl; 
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_I_RA_local[i].empty() )
						delivery_queue_I_RA_local[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_I_RA_local[i].begin(), delivery_queue_I_RA_local[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_I_RA_local[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
			}
			
			//////////////////////////////////
			// HVC-I -> HVC-RA delivery queue
			//////////////////////////////////
			// check if interneuron spike was delivered to some HVC-RA neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
			std::vector<std::pair<double,int>>::iterator it = delivery_queue_I_RA_local[i].begin();
			
			for (; it != delivery_queue_I_RA_local[i].end(); it++)
			{
				if (internal_time >= it->first)
				{
					some_RA_inh_conductance_was_updated_local = 1;
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_I_RA_local[i][pos_in_local_target_array];
					
					//std::cout << "Delivered spike HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << target_id << " delivered_time " << internal_time << " delivery time in queue = " << it->first << std::endl; 
					
					
					update_Gi_RA_local[target_id] += weights_I_RA_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(I) neuron " << Id_I_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to neuron " << syn_ID_I_RA_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_I_RA_local[i].erase(delivery_queue_I_RA_local[i].begin(), it);
		
		}  // end for i = 0 -> N_I_local (loop through all local HVC-I neurons)
		
		/////////////////////////////
		// Network Synchronization //
		/////////////////////////////
		if ( ( internal_time > network_time ) || ( t == num_steps-1 ) )
        {
            //if (MPI_rank == 0)
            //{
            //    std::cout << "internal time = " << internal_time << " > network_time = " << network_time << std::endl; 
            //}

            MPI_Allreduce(&some_RA_exc_conductance_was_updated_local, &some_RA_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_RA_inh_conductance_was_updated_local, &some_RA_inh_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_I_exc_conductance_was_updated_local, &some_I_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_RA_neuron_new_postEvent_local, &some_RA_neuron_new_postEvent_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            
        
			
            ////////////////////////////////////////////
            //// Process HVC-I -> HVC-RA interactions
            ////////////////////////////////////////////
            if ( some_RA_inh_conductance_was_updated_global > 0 )
            {
            // sum update array and send to all processes

                MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                for (int i = 0; i < N_RA_local; i++)
                {
					//std::cout << "HVC-RA neuron " << Id_RA_local[i] << " raised inhibitory conductance by " << update_Gi_RA_global[Id_RA_local[i]] << std::endl;
                    HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
            	}
				// update conductance arrays and delivered indicators
				some_RA_inh_conductance_was_updated_global = 0;
				some_RA_inh_conductance_was_updated_local = 0;
				
				std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            }
            
            ////////////////////////////////////////////
            //// Process HVC-RA -> HVC-I interactions
            ////////////////////////////////////////////
			if ( some_I_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons

				for (int i = 0; i < N_I_local; i++)
					HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);

				// update conductance arrays and fired indicators
				some_I_exc_conductance_was_updated_local = 0;
				some_I_exc_conductance_was_updated_global = 0;

				std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
			}
			
			////////////////////////////////////////////
            //// Process HVC-RA -> HVC-RA interactions
            ////////////////////////////////////////////
			if ( some_RA_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons
				for (int i = 0; i < N_RA_local; i++)
					HVCRA_local[i].raiseExcWeight(update_Ge_RA_global[Id_RA_local[i]]); // update conductance

				// update conductance arrays and fired indicators
				some_RA_exc_conductance_was_updated_global = 0;
				some_RA_exc_conductance_was_updated_local = 0;
				
				std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
			}	

			////////////////////////////////////////////
            //// Process HVC-RA dendritic spikes
            ////////////////////////////////////////////
			if ( some_RA_neuron_new_postEvent_global > 0 )
            {
				// update conductance arrays and fired indicators
            	some_RA_neuron_new_postEvent_local = 0;
	        	some_RA_neuron_new_postEvent_global = 0;
	        	
	        	// gather all HVC-RA neurons that produced new postsynaptic event
				this->gather_spiked_or_bursted_neurons(RA_neurons_new_postEvent_local, RA_neurons_new_postEvent_global);
				
				// update postsynaptic event arrays
				for (size_t j = 0; j < RA_neurons_new_postEvent_global.size(); j++)
					previous_postEvent_times_global[RA_neurons_new_postEvent_global[j]].push_back(internal_time);
				
				////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// LTP or LTD of previously delivered somatic spike times on HVC-RA neuron with new postsynaptic events  ///
				////////////////////////////////////////////////////////////////////////////////////////////////////////////
				
				for (int i = 0; i < N_RA_local; i++)
                {
					int presyn_ID = Id_RA_local[i]; // real id of presynaptic neuron
					
                    // if neuron is saturated apply LTP only if spiked neurons are among supersynapse targets
                    if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
                    {
                        for (size_t j = 0; j < RA_neurons_new_postEvent_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_new_postEvent_global[j]; // id of postsynaptic neuron
                            
                            // do not allow self-to-self synapses to emerge
							if ( presyn_ID != postsyn_ID )
							{
								std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
											supersynapses_local[i].end(), postsyn_ID );

								if ( pos!=supersynapses_local[i].end() )
								{
									// find previous somatic spikes that are not too old
									if ( !delivered_preEvent_times_local[i][postsyn_ID].empty() )
									{
										std::vector<double>::iterator it_relevant_preEvents = std::lower_bound(delivered_preEvent_times_local[i][postsyn_ID].begin(),
																											delivered_preEvent_times_local[i][postsyn_ID].end(),
																											internal_time - STDP_WINDOW);
																							
										for (auto it = it_relevant_preEvents; it != delivered_preEvent_times_local[i][postsyn_ID].end(); it++)
										{
											double dt = internal_time - *it;
									
										
											//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
											
											if (dt <= synaptic_params.T_0)
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												LTD(weights_RA_RA_local[i][postsyn_ID], dt);
												
												//~ std::cout   << "LTD from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
											}
											else
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
											
												LTP(weights_RA_RA_local[i][postsyn_ID], dt);
												
												//~ std::cout   << "LTP from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
												
											}	
											//double w = weights_local[i][fired_ID];
											
											//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
											 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
											 //           dt, weights_local[i][fired_ID] - w);
											
											update_synapse(i, postsyn_ID);
											
										}
									}
								}
							}
                        }

                    }
                    // if not saturated apply LTP for all dendritic spikes
                    else
                    {
                        for (size_t j = 0; j < RA_neurons_new_postEvent_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_new_postEvent_global[j];
                            // don't allow self-to-self connections
                            if ( postsyn_ID != presyn_ID )
                            {
                                // find previous somatic spikes that are not too old
								if ( !delivered_preEvent_times_local[i][postsyn_ID].empty() )
								{
									std::vector<double>::iterator it_relevant_preEvents = std::lower_bound(delivered_preEvent_times_local[i][postsyn_ID].begin(),
																										delivered_preEvent_times_local[i][postsyn_ID].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it = it_relevant_preEvents; it != delivered_preEvent_times_local[i][postsyn_ID].end(); it++)
									{
										double dt = internal_time - *it;
								

										//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
										
									    if (dt <= synaptic_params.T_0)
									    {
											double w_before = weights_RA_RA_local[i][postsyn_ID];
												
											LTD(weights_RA_RA_local[i][postsyn_ID], dt);
											
											// show LTD on pool neurons only
											//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											//~ 
											//~ if (p.first == p.second)
											//~ {
												//~ std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											//~ }
											
											//~ std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
											
										}
										else
										{
											double w_before = weights_RA_RA_local[i][postsyn_ID];
											
											LTP(weights_RA_RA_local[i][postsyn_ID], dt);
											
											
											// show LTP on pool neurons only
											//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											//~ 
											//~ if (p.first == p.second)
											//~ {
												//~ std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											//~ }
											
											
											//~ std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
										}	
										//double w = weights_local[i][fired_ID];
										
										//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
										 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
										 //           dt, weights_local[i][fired_ID] - w);
										
										update_synapse(i, postsyn_ID);
									
									}
								}
                            }
                        }
                   }

                } // end for i -> N_RA_local

				       

               // check if we need axon remodeling

                for (int i = 0; i < N_RA_local; i++)
                {
                    if ( (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss) && (remodeled_local[i] == 0) )
					    this->axon_remodeling(i);

				   
                }
            
				// clear bursts
				RA_neurons_new_postEvent_local.clear();
				RA_neurons_new_postEvent_global.clear();
				
            } // end if some HVC-RA produced new postsynaptic event
            
            
            network_time += NETWORK_UPDATE_FREQUENCY;
        }
    }
    
    this->potentiation_decay();
    //printf("After potentiation decay")
    //this->update_all_synapses_sudden_maturation();
    this->update_all_synapses();
}


void NetworkGrowthSimulator::trial_1stSoma_pre_1stSoma_post_delays_fixedSpread(bool training, std::vector<double>& spread_times)
{
	double event_window = 30.0; // window in which all somatic spikes are considered as one event
	
	// indicators for conductance updates and bursting
	int some_RA_inh_conductance_was_updated_global = 0;
	int some_RA_inh_conductance_was_updated_local = 0;
	
	int some_RA_exc_conductance_was_updated_global = 0;
	int some_RA_exc_conductance_was_updated_local = 0;
	
	int some_I_exc_conductance_was_updated_global = 0;
	int some_I_exc_conductance_was_updated_local = 0;
	
	int some_RA_neuron_new_postEvent_global = 0;
	int some_RA_neuron_new_postEvent_local = 0;
	
	 // arrays with updates for conductances
    std::vector<double> update_Ge_RA_local(N_RA);
	std::vector<double> update_Ge_RA_global(N_RA);

	std::vector<double> update_Gi_RA_local(N_RA);
	std::vector<double> update_Gi_RA_global(N_RA);

	std::vector<double> update_Ge_I_local(N_I);
	std::vector<double> update_Ge_I_global(N_I);
	
	std::vector<int> RA_neurons_new_postEvent_local; // local array of HVC-RA neurons that produced new postsynaptic event
	std::vector<int> RA_neurons_new_postEvent_global; // global array of HVC-RA neurons that produced new postsynaptic event
	
	std::vector<std::vector<std::vector<double>>> delivered_preEvent_times_local(N_RA_local); // local array with the times when presynaptic events 
															 // source HVC-RA neuron were delivered to other HVC-RA neurons
		
	for (int i = 0; i < N_RA_local; i++)
		delivered_preEvent_times_local[i].resize(N_RA);
		
		
	std::vector<std::vector<std::pair<double,int>>> delivery_queue_RA_RA_soma_local(N_RA_local); // local queue with spike delivery times for HVC-RA -> HVC-RA interactions
	std::vector<std::vector<std::pair<double,int>>> delivery_queue_RA_I_local(N_RA_local); // local queue with spike delivery times for HVC-RA -> HVC-I interactions
	std::vector<std::vector<std::pair<double,int>>> delivery_queue_I_RA_local(N_I_local); // local queue with spike delivery times for HVC-I -> HVC-RA interactions
		
	std::vector<std::vector<double>> previous_postEvent_times_global(N_RA); // array with the previous somatic spike times of neurons that are not too old
		
	// sample training innervation time and send to all processes
  // training neurons innervation
	std::vector<bool> indicators_training(N_RA);
	std::fill(indicators_training.begin(), indicators_training.end(), false);
	
	for (int i = 0; i < N_TR; i++)
		indicators_training[training_neurons[i]] = true;
		
	std::vector<bool> indicators_current_injected(N_RA);
	std::fill(indicators_current_injected.begin(), indicators_current_injected.end(), false);
	
	// sample mean training innervation time and send to all processes
    double mean_injection_time;
    
	if (MPI_rank == 0)
	{
		mean_injection_time = WAITING_TIME + noise_generator.random(TRIAL_DURATION - 2*WAITING_TIME);
		std::cout << "mean_injection_time = " << mean_injection_time << std::endl;
    }
    MPI_Bcast(&mean_injection_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
    double internal_time = 0;
    double network_time = NETWORK_UPDATE_FREQUENCY;

   

    //printf("network time = %f\n", network_time);

    
	// initialize update arrays and fired indicators
	std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
	
	
    // evolve dynamics
    int num_steps = static_cast<int>(TRIAL_DURATION / TIMESTEP);
    
    for (int t = 1; t < num_steps; t++)
	{
		internal_time += TIMESTEP;
		
		//////////////////////
		// HVC-RA neurons
		//////////////////////
		
		for (int i = 0; i < N_RA_local; i++)
		{
			if ( (training) && ( indicators_training[Id_RA_local[i]] ) && ( !indicators_current_injected[Id_RA_local[i]] ) )
			{
				if ( internal_time > spread_times[Id_RA_local[i]] + mean_injection_time )
				{
					HVCRA_local[i].raiseE(G_TRAINING_KICK);
					indicators_current_injected[Id_RA_local[i]] = true;
				}
			}
			
            // Debraband step
            HVCRA_local[i].Debraband_step_no_target_update();
            
           
			///////////////////////////////
			// HVC-RA somatic spike
			///////////////////////////////
			// if some neuron produced somatic spike
            if (HVCRA_local[i].get_fired_soma())
            {
                spikes_in_trial_soma_local[i].push_back(internal_time);

				// update somatic spikes only if current spike occured later than previous 
				// somatic spike + event_window
				if ( previous_postEvent_times_global[Id_RA_local[i]].empty() )
				{
					RA_neurons_new_postEvent_local.push_back(Id_RA_local[i]);
					some_RA_neuron_new_postEvent_local = 1;
					//std::cout << "HVC-RA " << Id_RA_local[i] << " postsynaptic event at " << internal_time << std::endl;
				}
				else if ( internal_time > previous_postEvent_times_global[Id_RA_local[i]].back() + event_window)
				{
					RA_neurons_new_postEvent_local.push_back(Id_RA_local[i]);
					some_RA_neuron_new_postEvent_local = 1;
					//std::cout << "HVC-RA " << Id_RA_local[i] << " postsynaptic event at " << internal_time << std::endl;
				}
							
				// for inhibitory neurons
				// loop over all inhibitory targets of fired neurons
				size_t num_I_targets = syn_ID_RA_I_local[i].size();
				
				for (size_t j = 0; j < num_I_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_RA_I_local[i][j];
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_RA_I_local[i].empty() )
						delivery_queue_RA_I_local[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_RA_I_local[i].begin(), delivery_queue_RA_I_local[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_RA_I_local[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
					//std::cout << "neuron " << Id_RA_local[i] << " spike time = " << internal_time << " axonal delay = " << axonal_delays_RA_I[Id_RA_local[i]][syn_ID] << " delivery time RA to I: " << delivery_queue_RA_I[i].back().first << " delivery target id: " << syn_ID_RA_I_local[i][delivery_queue_RA_I[i].back().second] << std::endl;
				
				
				// if neuron is saturated, deliver spikes only to supersynaptic targets
				if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
				{
					for (size_t j = 0; j < supersynapses_local[i].size(); j++)
					{
						int target_id = supersynapses_local[i][j];
						
						double delivery_time = internal_time + axonal_delays_RA_RA_local[i][target_id];
						
						// if queue is empty, just add item to the queue
						if ( delivery_queue_RA_RA_soma_local[i].empty() )
							delivery_queue_RA_RA_soma_local[i].push_back(std::pair<double,int>(delivery_time, target_id));	
						// otherwise add item so that queue is sorted
						else
						{
							auto it = std::upper_bound(delivery_queue_RA_RA_soma_local[i].begin(), delivery_queue_RA_RA_soma_local[i].end(), std::pair<double,int>(delivery_time, target_id));
							delivery_queue_RA_RA_soma_local[i].insert(it, std::pair<double,int>(delivery_time, target_id));
						}
					}
				}
				else // deliver spikes to everyone except itself
				{	
					for (int j = 0; j < N_RA; j++)
					{
						if ( j != Id_RA_local[i])
						{
						
							double delivery_time = internal_time + axonal_delays_RA_RA_local[i][j];
						
							// if queue is empty, just add item to the queue
							if ( delivery_queue_RA_RA_soma_local[i].empty() )
								delivery_queue_RA_RA_soma_local[i].push_back(std::pair<double,int>(delivery_time, j));	
							// otherwise add item so that queue is sorted
							else
							{
								auto it = std::upper_bound(delivery_queue_RA_RA_soma_local[i].begin(), delivery_queue_RA_RA_soma_local[i].end(), std::pair<double,int>(delivery_time, j));
								delivery_queue_RA_RA_soma_local[i].insert(it, std::pair<double,int>(delivery_time, j));
							}
						}
					}
				}
			
            } // end if local HVC-RA neuron spiked  
            
             //////////////////////////////////
			// HVC-RA -> HVC-I delivery queue
			//////////////////////////////////
            // check if somatic spike was delivered to some interneuron
            // loop through the delivery queue to check if current time exceeds the spike delivery time
            std::vector<std::pair<double,int>>::iterator it = delivery_queue_RA_I_local[i].begin();
			
            for (; it != delivery_queue_RA_I_local[i].end(); it++)
            {
				if (internal_time >= it->first)
				{
					
					some_I_exc_conductance_was_updated_local = 1;
					
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_RA_I_local[i][pos_in_local_target_array];
					
					update_Ge_I_local[target_id] += weights_RA_I_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(RA) neuron " << Id_RA_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to HVC(I) neuron " << syn_ID_RA_I_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_I_local[i].erase(delivery_queue_RA_I_local[i].begin(), it);
			
			//////////////////////////////////////////////////
			// HVC-RA -> HVC-RA somatic spike delivery queue
			//////////////////////////////////////////////////
			// check if somatic spike was delivered to some HVC(RA) neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
            it = delivery_queue_RA_RA_soma_local[i].begin();
            
            for (; it != delivery_queue_RA_RA_soma_local[i].end(); it++)
            {
				double delivery_time = it->first;
				
				if ( internal_time >= delivery_time )
				{
					int target_id = it->second;
					
					// update conductance of target if synapse is active
					if ( active_indicators_local[i][target_id] == 1 )
					{
						some_RA_exc_conductance_was_updated_local = 1;
						update_Ge_RA_local[target_id] += weights_RA_RA_local[i][target_id];
					}
					
					//////////////
					//// LTD
					//////////////
					// update delivered spikes only if current spike occured later than previous delivered
					// spike + event_window
					bool new_preEvent_occured = false;
				
					if ( !delivered_preEvent_times_local[i][target_id].empty() ){
						if ( delivery_time > delivered_preEvent_times_local[i][target_id].back() + event_window ){
							delivered_preEvent_times_local[i][target_id].push_back(delivery_time);
							new_preEvent_occured = true;
						}
					}
					else{
						delivered_preEvent_times_local[i][target_id].push_back(delivery_time);
						new_preEvent_occured = true;
					}
					
					if ( new_preEvent_occured ){
						 // if neuron is saturated apply LTD only if target is among super synapses
						if (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss)
						{
							if ( supersynapses_indicators_local[i][target_id] == 1 )
							{	
								// find previous dendritic spikes that are not too old
								if ( !previous_postEvent_times_global[target_id].empty() )
								{
									std::vector<double>::iterator it_relevant_postEvent = std::lower_bound(previous_postEvent_times_global[target_id].begin(),
																										previous_postEvent_times_global[target_id].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it_postEvent = it_relevant_postEvent; it_postEvent != previous_postEvent_times_global[target_id].end(); it_postEvent++)
									{
										double dt = *it_postEvent - delivery_time;

										//std::cout << "dt in saturated LTD = " << dt << std::endl;

										
										double w_before = weights_RA_RA_local[i][target_id];     
							   
										LTD(weights_RA_RA_local[i][target_id], dt);
										
										
										//~ std::cout << "LTD from saturated " << Id_RA_local[i] << " -> " << supersynapse_id
												  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][supersynapse_id] - w_before
												  //~ << std::endl;
												  
										//printf("LTD from saturated neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n",
										//            Id_RA_local[i], supersynapse_id, spikes_in_trial_soma_local[i].back(), spike_times_dend_global[supersynapse_id],
										//            dt, weights_local[i][supersynapse_id] - w);
										
										update_synapse(i, target_id);	
										
									}   
								}
							}

							// if some supersynapse desaturated, update all synapses
							if (static_cast<int>(supersynapses_local[i].size()) < synaptic_params.Nss)
								for (int j = 0; j < N_RA; j++)
									this->update_synapse(i, j);
						}
						// if not saturated apply LTD rule 
						else
						{
							// find previous dendritic spikes that are not too old
							if ( !previous_postEvent_times_global[target_id].empty() )
							{
								std::vector<double>::iterator it_relevant_postEvent = std::lower_bound(previous_postEvent_times_global[target_id].begin(),
																									previous_postEvent_times_global[target_id].end(),
																									internal_time - STDP_WINDOW);
																					
								for (auto it_postEvent = it_relevant_postEvent; it_postEvent != previous_postEvent_times_global[target_id].end(); it_postEvent++)
								{
									double dt = *it_postEvent - delivery_time;

									//std::cout << "dt in saturated LTD = " << dt << std::endl;

									
									double w_before = weights_RA_RA_local[i][target_id];     
						   
									LTD(weights_RA_RA_local[i][target_id], dt);
									
									//~ // show LTD to pool neurons only
									//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), target_id);
									//~ 
									//~ if (p.first == p.second)
									//~ {
										//~ std::cout << "LTD from  " << Id_RA_local[i] << " -> " << target_id
											  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][target_id] - w_before
											  //~ << " w_after = " << weights_RA_RA_local[i][target_id] << std::endl;
									//~ }
								//~ 
									update_synapse(i, target_id);	
									
								}   
							}
						}
					}
				} // end if internal_time >= delivery spike time 
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_RA_soma_local[i].erase(delivery_queue_RA_RA_soma_local[i].begin(), it);
			
			
			///////////////////////////////
			// HVC-RA dendritic spike
			///////////////////////////////
			//if some neuron produced dendritic spike, store this neuron in array
			
			if ( HVCRA_local[i].get_fired_dend() )
			{
				//~ // show spikes of pool neurons only
				//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
				//~ 
				//~ if (p.first == p.second)
					//~ std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
				//~ 
				//std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
				
				spikes_in_trial_dend_local[i].push_back(internal_time); // dend spike time relative to trial onset
				
				
				//std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
			} // end if get_fired_dend
	
            
		} // end for i = 0 -> N_RA_local (loop through all local HVC-RA neurons)
		
		//////////////////////
		// HVC-I neurons
		//////////////////////
		for (int i = 0; i < N_I_local; i++)
		{
			HVCI_local[i].DP8_step_no_target_update();
			
			
			///////////////////////////////
			// HVC-I spike
			///////////////////////////////
			//  if some I neuron spikes, update delivery queue
			if (HVCI_local[i].get_fired())
			{
				//printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
				spikes_in_trial_interneuron_local[i].push_back(internal_time);

				size_t num_RA_targets = syn_ID_I_RA_local[i].size();
				// loop over all targets of fired neurons
				for (size_t j = 0; j < num_RA_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_I_RA_local[i][j];
					
					//std::cout << "HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << syn_ID_I_RA_local[i][j] << " spike time " << internal_time << " delivery time = " << delivery_time << std::endl; 
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_I_RA_local[i].empty() )
						delivery_queue_I_RA_local[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_I_RA_local[i].begin(), delivery_queue_I_RA_local[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_I_RA_local[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
			}
			
			//////////////////////////////////
			// HVC-I -> HVC-RA delivery queue
			//////////////////////////////////
			// check if interneuron spike was delivered to some HVC-RA neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
			std::vector<std::pair<double,int>>::iterator it = delivery_queue_I_RA_local[i].begin();
			
			for (; it != delivery_queue_I_RA_local[i].end(); it++)
			{
				if (internal_time >= it->first)
				{
					some_RA_inh_conductance_was_updated_local = 1;
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_I_RA_local[i][pos_in_local_target_array];
					
					//std::cout << "Delivered spike HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << target_id << " delivered_time " << internal_time << " delivery time in queue = " << it->first << std::endl; 
					
					
					update_Gi_RA_local[target_id] += weights_I_RA_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(I) neuron " << Id_I_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to neuron " << syn_ID_I_RA_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_I_RA_local[i].erase(delivery_queue_I_RA_local[i].begin(), it);
		
		}  // end for i = 0 -> N_I_local (loop through all local HVC-I neurons)
		
		/////////////////////////////
		// Network Synchronization //
		/////////////////////////////
		if ( ( internal_time > network_time ) || ( t == num_steps-1 ) )
        {
            //if (MPI_rank == 0)
            //{
            //    std::cout << "internal time = " << internal_time << " > network_time = " << network_time << std::endl; 
            //}

            MPI_Allreduce(&some_RA_exc_conductance_was_updated_local, &some_RA_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_RA_inh_conductance_was_updated_local, &some_RA_inh_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_I_exc_conductance_was_updated_local, &some_I_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_RA_neuron_new_postEvent_local, &some_RA_neuron_new_postEvent_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            
        
			
            ////////////////////////////////////////////
            //// Process HVC-I -> HVC-RA interactions
            ////////////////////////////////////////////
            if ( some_RA_inh_conductance_was_updated_global > 0 )
            {
            // sum update array and send to all processes

                MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                for (int i = 0; i < N_RA_local; i++)
                {
					//std::cout << "HVC-RA neuron " << Id_RA_local[i] << " raised inhibitory conductance by " << update_Gi_RA_global[Id_RA_local[i]] << std::endl;
                    HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
            	}
				// update conductance arrays and delivered indicators
				some_RA_inh_conductance_was_updated_global = 0;
				some_RA_inh_conductance_was_updated_local = 0;
				
				std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            }
            
            ////////////////////////////////////////////
            //// Process HVC-RA -> HVC-I interactions
            ////////////////////////////////////////////
			if ( some_I_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons

				for (int i = 0; i < N_I_local; i++)
					HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);

				// update conductance arrays and fired indicators
				some_I_exc_conductance_was_updated_local = 0;
				some_I_exc_conductance_was_updated_global = 0;

				std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
			}
			
			////////////////////////////////////////////
            //// Process HVC-RA -> HVC-RA interactions
            ////////////////////////////////////////////
			if ( some_RA_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons
				for (int i = 0; i < N_RA_local; i++)
					HVCRA_local[i].raiseExcWeight(update_Ge_RA_global[Id_RA_local[i]]); // update conductance

				// update conductance arrays and fired indicators
				some_RA_exc_conductance_was_updated_global = 0;
				some_RA_exc_conductance_was_updated_local = 0;
				
				std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
			}	

			////////////////////////////////////////////
            //// Process HVC-RA dendritic spikes
            ////////////////////////////////////////////
			if ( some_RA_neuron_new_postEvent_global > 0 )
            {
				// update conductance arrays and fired indicators
            	some_RA_neuron_new_postEvent_local = 0;
	        	some_RA_neuron_new_postEvent_global = 0;
	        	
	        	// gather all HVC-RA neurons that produced new postsynaptic event
				this->gather_spiked_or_bursted_neurons(RA_neurons_new_postEvent_local, RA_neurons_new_postEvent_global);
				
				// update postsynaptic event arrays
				for (size_t j = 0; j < RA_neurons_new_postEvent_global.size(); j++)
					previous_postEvent_times_global[RA_neurons_new_postEvent_global[j]].push_back(internal_time);
				
				////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// LTP or LTD of previously delivered somatic spike times on HVC-RA neuron with new postsynaptic events  ///
				////////////////////////////////////////////////////////////////////////////////////////////////////////////
				
				for (int i = 0; i < N_RA_local; i++)
                {
					int presyn_ID = Id_RA_local[i]; // real id of presynaptic neuron
					
                    // if neuron is saturated apply LTP only if spiked neurons are among supersynapse targets
                    if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
                    {
                        for (size_t j = 0; j < RA_neurons_new_postEvent_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_new_postEvent_global[j]; // id of postsynaptic neuron
                            
                            // do not allow self-to-self synapses to emerge
							if ( presyn_ID != postsyn_ID )
							{
								std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
											supersynapses_local[i].end(), postsyn_ID );

								if ( pos!=supersynapses_local[i].end() )
								{
									// find previous somatic spikes that are not too old
									if ( !delivered_preEvent_times_local[i][postsyn_ID].empty() )
									{
										std::vector<double>::iterator it_relevant_preEvents = std::lower_bound(delivered_preEvent_times_local[i][postsyn_ID].begin(),
																											delivered_preEvent_times_local[i][postsyn_ID].end(),
																											internal_time - STDP_WINDOW);
																							
										for (auto it = it_relevant_preEvents; it != delivered_preEvent_times_local[i][postsyn_ID].end(); it++)
										{
											double dt = internal_time - *it;
									
										
											//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
											
											if (dt <= synaptic_params.T_0)
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												LTD(weights_RA_RA_local[i][postsyn_ID], dt);
												
												//~ std::cout   << "LTD from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
											}
											else
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
											
												LTP(weights_RA_RA_local[i][postsyn_ID], dt);
												
												//~ std::cout   << "LTP from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
												
											}	
											//double w = weights_local[i][fired_ID];
											
											//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
											 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
											 //           dt, weights_local[i][fired_ID] - w);
											
											update_synapse(i, postsyn_ID);
											
										}
									}
								}
							}
                        }

                    }
                    // if not saturated apply LTP for all dendritic spikes
                    else
                    {
                        for (size_t j = 0; j < RA_neurons_new_postEvent_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_new_postEvent_global[j];
                            // don't allow self-to-self connections
                            if ( postsyn_ID != presyn_ID )
                            {
                                // find previous somatic spikes that are not too old
								if ( !delivered_preEvent_times_local[i][postsyn_ID].empty() )
								{
									std::vector<double>::iterator it_relevant_preEvents = std::lower_bound(delivered_preEvent_times_local[i][postsyn_ID].begin(),
																										delivered_preEvent_times_local[i][postsyn_ID].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it = it_relevant_preEvents; it != delivered_preEvent_times_local[i][postsyn_ID].end(); it++)
									{
										double dt = internal_time - *it;
								

										//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
										
									    if (dt <= synaptic_params.T_0)
									    {
											double w_before = weights_RA_RA_local[i][postsyn_ID];
												
											LTD(weights_RA_RA_local[i][postsyn_ID], dt);
											
											// show LTD on pool neurons only
											//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											//~ 
											//~ if (p.first == p.second)
											//~ {
												//~ std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											//~ }
											
											//~ std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
											
										}
										else
										{
											double w_before = weights_RA_RA_local[i][postsyn_ID];
											
											LTP(weights_RA_RA_local[i][postsyn_ID], dt);
											
											
											// show LTP on pool neurons only
											//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											//~ 
											//~ if (p.first == p.second)
											//~ {
												//~ std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											//~ }
											
											
											//~ std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
										}	
										//double w = weights_local[i][fired_ID];
										
										//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
										 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
										 //           dt, weights_local[i][fired_ID] - w);
										
										update_synapse(i, postsyn_ID);
									
									}
								}
                            }
                        }
                   }

                } // end for i -> N_RA_local

				       

               // check if we need axon remodeling

                for (int i = 0; i < N_RA_local; i++)
                {
                    if ( (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss) && (remodeled_local[i] == 0) )
					    this->axon_remodeling(i);

				   
                }
            
				// clear bursts
				RA_neurons_new_postEvent_local.clear();
				RA_neurons_new_postEvent_global.clear();
				
            } // end if some HVC-RA produced new postsynaptic event
            
            
            network_time += NETWORK_UPDATE_FREQUENCY;
        }
    }
    
    this->potentiation_decay();
    //printf("After potentiation decay")
    //this->update_all_synapses_sudden_maturation();
    this->update_all_synapses();
}

void NetworkGrowthSimulator::check_neuron_activity(std::vector<int>& neurons_to_replace)
{
	std::vector<int> neurons_to_replace_local;
	
	// check if some neurons became silent or driven
	for (int i = 0; i < N_RA_local; i++)
	{
		num_spikes_in_recent_trials_local[i].push_front(static_cast<int>(spikes_in_trial_soma_local[i].size()));
        
        //firing_rate_short_local[i] = std::accumulate(num_spikes_in_recent_trials[i].begin(), num_spikes_in_recent_trials[i].begin() + RATE_WINDOW_SHORT, 0.0)
          //                          / static_cast<double>(RATE_WINDOW_SHORT);

		// calculate firing rate in large window:
		int window = 1000;
		
		double driven_threshold = 0.50;
		
        firing_rate_long_local[i] = std::accumulate(num_spikes_in_recent_trials_local[i].begin(), num_spikes_in_recent_trials_local[i].end(), 0.0)
                                                                                            / static_cast<double>(maturation_params.RATE_WINDOW_LONG);
        	
		std::vector<double> recent_firing_robustness(window);
		
		for (int j = 0; j < window; j++ )
			if ( num_spikes_in_recent_trials_local[i][j] > 0 )
				recent_firing_robustness[j] = 1;
		
		double firing_robustness = std::accumulate(recent_firing_robustness.begin(), recent_firing_robustness.end(), 0.0) / static_cast<double>(window);
		
		if ( (firing_robustness >= driven_threshold) && (maturation_scale_local[i] != MATURATION_SCALE_DRIVEN) )
		{
			std::cout << "Neuron " << Id_RA_local[i] << " became driven!" << std::endl;
			maturation_scale_local[i] = MATURATION_SCALE_DRIVEN;
		}
		
		if ( ( firing_robustness < maturation_params.DEATH_RATE_THRESHOLD ) && ( num_trials_after_replacement_local[i] > window ) )
		{
			std::cout << "Neuron " << Id_RA_local[i] << " became silent!" << std::endl;
			
			neurons_to_replace_local.push_back(Id_RA_local[i]);
		}                                                                                                                          
			//if (std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].end(), 0.0) >= 1)
			//{
			//	std::cout << "Recent bursts of neuron " << Id_RA_local[i] << " :" << std::endl;
				
			//	for (size_t j = 0; j < num_bursts_in_recent_trials[i].size(); j++)
			//		std::cout << num_bursts_in_recent_trials[i][j] << "\t";
				
			//	std::cout << std::endl;
			//}
	}
	
	this->gather_localVectorsToGlobal(neurons_to_replace_local, neurons_to_replace);
	
	
}

void NetworkGrowthSimulator::trial_burst_pre_dend_event_post_delays_sudden_maturation_noImmatureOut_fixedSpread(bool training, std::vector<double>& spread_times)
{
	double event_window = 30.0; // window in which all somatic spikes are considered as one event
	
	// indicators for conductance updates and bursting
	int some_RA_inh_conductance_was_updated_global = 0;
	int some_RA_inh_conductance_was_updated_local = 0;
	
	int some_RA_exc_conductance_was_updated_global = 0;
	int some_RA_exc_conductance_was_updated_local = 0;
	
	int some_I_exc_conductance_was_updated_global = 0;
	int some_I_exc_conductance_was_updated_local = 0;
	
	int some_RA_neuron_bursted_global = 0;
	int some_RA_neuron_bursted_local = 0;
	
	 // arrays with updates for conductances
    std::vector<double> update_Ge_RA_local(N_RA);
	std::vector<double> update_Ge_RA_global(N_RA);

	std::vector<double> update_Gi_RA_local(N_RA);
	std::vector<double> update_Gi_RA_global(N_RA);

	std::vector<double> update_Ge_I_local(N_I);
	std::vector<double> update_Ge_I_global(N_I);
	
	std::vector<int> RA_neurons_bursted_local; // local array of HVC-RA neurons that bursted
	std::vector<int> RA_neurons_bursted_global; // global array of HVC-RA neurons that bursted
	
	// sample training innervation time and send to all processes
  // training neurons innervation
	std::vector<bool> indicators_training(N_RA);
	std::fill(indicators_training.begin(), indicators_training.end(), false);
	
	for (int i = 0; i < N_TR; i++)
		indicators_training[training_neurons[i]] = true;
		
	std::vector<bool> indicators_current_injected(N_RA);
	std::fill(indicators_current_injected.begin(), indicators_current_injected.end(), false);
	
	// sample mean training innervation time and send to all processes
    double mean_injection_time;
    
	if (MPI_rank == 0)
	{
		mean_injection_time = WAITING_TIME + noise_generator.random(TRIAL_DURATION - 2*WAITING_TIME);
		std::cout << "mean_injection_time = " << mean_injection_time << std::endl;
    }
    MPI_Bcast(&mean_injection_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
    double internal_time = 0;
    double network_time = NETWORK_UPDATE_FREQUENCY;

   

    //printf("network time = %f\n", network_time);

    
	// initialize update arrays and fired indicators
	std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
	
	
    // evolve dynamics
    int num_steps = static_cast<int>(TRIAL_DURATION / TIMESTEP);
    
    for (int t = 1; t < num_steps; t++)
	{
		internal_time += TIMESTEP;
		
		//////////////////////
		// HVC-RA neurons
		//////////////////////
		
		for (int i = 0; i < N_RA_local; i++)
		{
			if ( (training) && ( indicators_training[Id_RA_local[i]] ) && ( !indicators_current_injected[Id_RA_local[i]] ) )
			{
				if ( internal_time > spread_times[Id_RA_local[i]] + mean_injection_time )
				{
					HVCRA_local[i].raiseE(G_TRAINING_KICK);
					indicators_current_injected[Id_RA_local[i]] = true;
				}
			}
			
            // Debraband step
            HVCRA_local[i].Debraband_step_no_target_update();
            
           
			///////////////////////////////
			// HVC-RA somatic spike
			///////////////////////////////
			// if some neuron produced somatic spike
            if (HVCRA_local[i].get_fired_soma())
            {
                spikes_in_trial_soma_local[i].push_back(internal_time);
                
				// update delivery queues if neuron is mature (no output for immature neuron and no wiring mechanism)
				//if ( mature_global[Id_RA_local[i]] == 1 ){
				
					// for inhibitory neurons
					// loop over all inhibitory targets of fired neurons
					size_t num_I_targets = syn_ID_RA_I_local[i].size();
					
					for (size_t j = 0; j < num_I_targets; j++)
					{
						double delivery_time = internal_time + axonal_delays_RA_I_local[i][j];
						
						// if queue is empty, just add item to the queue
						if ( delivery_queue_RA_I[i].empty() )
							delivery_queue_RA_I[i].push_back(std::pair<double,int>(delivery_time, j));	
						// otherwise add item so that queue is sorted
						else
						{
							auto it = std::upper_bound(delivery_queue_RA_I[i].begin(), delivery_queue_RA_I[i].end(), std::pair<double,int>(delivery_time, j));
							delivery_queue_RA_I[i].insert(it, std::pair<double,int>(delivery_time, j));
						}
					}
						//std::cout << "neuron " << Id_RA_local[i] << " spike time = " << internal_time << " axonal delay = " << axonal_delays_RA_I[Id_RA_local[i]][syn_ID] << " delivery time RA to I: " << delivery_queue_RA_I[i].back().first << " delivery target id: " << syn_ID_RA_I_local[i][delivery_queue_RA_I[i].back().second] << std::endl;
					
					
					// if neuron is saturated, deliver spikes only to supersynaptic targets
					if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
					{
						for (size_t j = 0; j < supersynapses_local[i].size(); j++)
						{
							int target_id = supersynapses_local[i][j];
							
							double delivery_time = internal_time + axonal_delays_RA_RA_local[i][target_id];
							
							// if queue is empty, just add item to the queue
							if ( delivery_queue_RA_RA_soma[i].empty() )
								delivery_queue_RA_RA_soma[i].push_back(std::pair<double,int>(delivery_time, target_id));	
							// otherwise add item so that queue is sorted
							else
							{
								auto it = std::upper_bound(delivery_queue_RA_RA_soma[i].begin(), delivery_queue_RA_RA_soma[i].end(), std::pair<double,int>(delivery_time, target_id));
								delivery_queue_RA_RA_soma[i].insert(it, std::pair<double,int>(delivery_time, target_id));
							}
						}
					}
					else // deliver spikes to everyone except itself
					{	
						for (int j = 0; j < N_RA; j++)
						{
							if ( j != Id_RA_local[i])
							{
							
								double delivery_time = internal_time + axonal_delays_RA_RA_local[i][j];
							
								// if queue is empty, just add item to the queue
								if ( delivery_queue_RA_RA_soma[i].empty() )
									delivery_queue_RA_RA_soma[i].push_back(std::pair<double,int>(delivery_time, j));	
								// otherwise add item so that queue is sorted
								else
								{
									auto it = std::upper_bound(delivery_queue_RA_RA_soma[i].begin(), delivery_queue_RA_RA_soma[i].end(), std::pair<double,int>(delivery_time, j));
									delivery_queue_RA_RA_soma[i].insert(it, std::pair<double,int>(delivery_time, j));
								}
							}
						}
					}
				//}
                
            } // end if local HVC-RA neuron spiked  
            
             //////////////////////////////////
			// HVC-RA -> HVC-I delivery queue
			//////////////////////////////////
            // check if somatic spike was delivered to some interneuron
            // loop through the delivery queue to check if current time exceeds the spike delivery time
            std::vector<std::pair<double,int>>::iterator it = delivery_queue_RA_I[i].begin();
			
            for (; it != delivery_queue_RA_I[i].end(); it++)
            {
				if (internal_time >= it->first)
				{
					
					some_I_exc_conductance_was_updated_local = 1;
					
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_RA_I_local[i][pos_in_local_target_array];
					
					update_Ge_I_local[target_id] += weights_RA_I_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(RA) neuron " << Id_RA_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to HVC(I) neuron " << syn_ID_RA_I_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_I[i].erase(delivery_queue_RA_I[i].begin(), it);
			
			//////////////////////////////////////////////////
			// HVC-RA -> HVC-RA somatic spike delivery queue
			//////////////////////////////////////////////////
			// check if somatic spike was delivered to some HVC(RA) neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
            it = delivery_queue_RA_RA_soma[i].begin();
            
            for (; it != delivery_queue_RA_RA_soma[i].end(); it++)
            {
				double delivery_time = it->first;
				
				if ( internal_time >= delivery_time )
				{
					int target_id = it->second;
					
					
					
					// update conductance of target if synapse is active
					if ( active_indicators_local[i][target_id] == 1 )
					{
						some_RA_exc_conductance_was_updated_local = 1;
						update_Ge_RA_local[target_id] += weights_RA_RA_local[i][target_id];
					}
					
					
					
					//////////////
					//// LTD
					//////////////
					// do LTD on dendritic spike time of target neuron if it is not mature neuron
					//~ if ( mature_global[target_id] != 1 )
					//~ {
						// update delivered spikes only if current spike occured later than previous delivered
						// spike + event_window
						bool new_event_occured = false;
					
						if ( !delivered_spike_times[i][target_id].empty() ){
							if ( delivery_time > delivered_spike_times[i][target_id].back() + event_window ){
								delivered_spike_times[i][target_id].push_back(delivery_time);
								new_event_occured = true;
							}
						}
						else{
							delivered_spike_times[i][target_id].push_back(delivery_time);
							new_event_occured = true;
						}
						
						if ( new_event_occured ){
							 // if neuron is saturated apply LTD only if target is among super synapses
							if (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss)
							{
								if ( supersynapses_indicators_local[i][target_id] == 1 )
								{	
									// find previous dendritic spikes that are not too old
									if ( !previous_dendritic_spike_times_global[target_id].empty() )
									{
										std::vector<double>::iterator it_relevant_dend_spikes = std::lower_bound(previous_dendritic_spike_times_global[target_id].begin(),
																											previous_dendritic_spike_times_global[target_id].end(),
																											internal_time - STDP_WINDOW);
																							
										for (auto it_dend_spike = it_relevant_dend_spikes; it_dend_spike != previous_dendritic_spike_times_global[target_id].end(); it_dend_spike++)
										{
											double dt = *it_dend_spike - delivery_time;

											//std::cout << "dt in saturated LTD = " << dt << std::endl;

											
											double w_before = weights_RA_RA_local[i][target_id];     
								   
											LTD(weights_RA_RA_local[i][target_id], dt);
											
											
											//~ std::cout << "LTD from saturated " << Id_RA_local[i] << " -> " << supersynapse_id
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][supersynapse_id] - w_before
													  //~ << std::endl;
													  
											//printf("LTD from saturated neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n",
											//            Id_RA_local[i], supersynapse_id, spikes_in_trial_soma_local[i].back(), spike_times_dend_global[supersynapse_id],
											//            dt, weights_local[i][supersynapse_id] - w);
											
											update_synapse(i, target_id);	
											
										}   
									}
								}

								// if some supersynapse desaturated, update all synapses
								if (static_cast<int>(supersynapses_local[i].size()) < synaptic_params.Nss)
									for (int j = 0; j < N_RA; j++)
										this->update_synapse(i, j);
							}
							// if not saturated apply LTD rule 
							else
							{
								// find previous dendritic spikes that are not too old
								if ( !previous_dendritic_spike_times_global[target_id].empty() )
								{
									std::vector<double>::iterator it_relevant_dend_spikes = std::lower_bound(previous_dendritic_spike_times_global[target_id].begin(),
																										previous_dendritic_spike_times_global[target_id].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it_dend_spike = it_relevant_dend_spikes; it_dend_spike != previous_dendritic_spike_times_global[target_id].end(); it_dend_spike++)
									{
										double dt = *it_dend_spike - delivery_time;

										//std::cout << "dt in saturated LTD = " << dt << std::endl;

										
										double w_before = weights_RA_RA_local[i][target_id];     
							   
										LTD(weights_RA_RA_local[i][target_id], dt);
										
										//~ // show LTD to pool neurons only
										//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), target_id);
										//~ 
										//~ if (p.first == p.second)
										//~ {
											//~ std::cout << "LTD from  " << Id_RA_local[i] << " -> " << target_id
												  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][target_id] - w_before
												  //~ << " w_after = " << weights_RA_RA_local[i][target_id] << std::endl;
										//~ }
									//~ 
										update_synapse(i, target_id);	
										
									}   
								}
							}
						}
					//~ } // end if target neuron is not mature
				} // end if internal_time >= delivery spike time 
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_RA_soma[i].erase(delivery_queue_RA_RA_soma[i].begin(), it);
			
			
			///////////////////////////////
			// HVC-RA dendritic spike
			///////////////////////////////
			//if some neuron produced dendritic spike, store this neuron in array
			
			if ( HVCRA_local[i].get_fired_dend() )
			{
				//~ // show spikes of pool neurons only
				//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
				//~ 
				//~ if (p.first == p.second)
					//~ std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
				//~ 
				std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
				
				spikes_in_trial_dend_local[i].push_back(internal_time); // dend spike time relative to trial onset
				
				//if ( mature_global[Id_RA_local[i]] != 1 )
				//{
				
				// update dendritic spikes only if current spike occured later than previous 
				// dendritic spike + event_window
				if ( previous_dendritic_spike_times_global[Id_RA_local[i]].empty() )
				{
					RA_neurons_bursted_local.push_back(Id_RA_local[i]);
					some_RA_neuron_bursted_local = 1;
					std::cout << "HVC-RA " << Id_RA_local[i] << " dendritic event at " << internal_time << std::endl;
				}
				else if ( internal_time > previous_dendritic_spike_times_global[Id_RA_local[i]].back() + event_window)
				{
					RA_neurons_bursted_local.push_back(Id_RA_local[i]);
					some_RA_neuron_bursted_local = 1;
					std::cout << "HVC-RA " << Id_RA_local[i] << " dendritic event at " << internal_time << std::endl;
				}
				//}
				//std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
			} // end if get_fired_dend
	
            
		} // end for i = 0 -> N_RA_local (loop through all local HVC-RA neurons)
		
		//////////////////////
		// HVC-I neurons
		//////////////////////
		for (int i = 0; i < N_I_local; i++)
		{
			HVCI_local[i].DP8_step_no_target_update();
			
			
			///////////////////////////////
			// HVC-I spike
			///////////////////////////////
			//  if some I neuron spikes, update delivery queue
			if (HVCI_local[i].get_fired())
			{
				//printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
				spikes_in_trial_interneuron_local[i].push_back(internal_time);

				size_t num_RA_targets = syn_ID_I_RA_local[i].size();
				// loop over all targets of fired neurons
				for (size_t j = 0; j < num_RA_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_I_RA_local[i][j];
					
					//std::cout << "HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << syn_ID_I_RA_local[i][j] << " spike time " << internal_time << " delivery time = " << delivery_time << std::endl; 
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_I_RA[i].empty() )
						delivery_queue_I_RA[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_I_RA[i].begin(), delivery_queue_I_RA[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_I_RA[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
			}
			
			//////////////////////////////////
			// HVC-I -> HVC-RA delivery queue
			//////////////////////////////////
			// check if interneuron spike was delivered to some HVC-RA neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
			std::vector<std::pair<double,int>>::iterator it = delivery_queue_I_RA[i].begin();
			
			for (; it != delivery_queue_I_RA[i].end(); it++)
			{
				if (internal_time >= it->first)
				{
					some_RA_inh_conductance_was_updated_local = 1;
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_I_RA_local[i][pos_in_local_target_array];
					
					//std::cout << "Delivered spike HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << target_id << " delivered_time " << internal_time << " delivery time in queue = " << it->first << std::endl; 
					
					
					update_Gi_RA_local[target_id] += weights_I_RA_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(I) neuron " << Id_I_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to neuron " << syn_ID_I_RA_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_I_RA[i].erase(delivery_queue_I_RA[i].begin(), it);
		
		}  // end for i = 0 -> N_I_local (loop through all local HVC-I neurons)
		
		/////////////////////////////
		// Network Synchronization //
		/////////////////////////////
		if ( ( internal_time > network_time ) || ( t == num_steps-1 ) )
        {
            //if (MPI_rank == 0)
            //{
            //    std::cout << "internal time = " << internal_time << " > network_time = " << network_time << std::endl; 
            //}

            MPI_Allreduce(&some_RA_exc_conductance_was_updated_local, &some_RA_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_RA_inh_conductance_was_updated_local, &some_RA_inh_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_I_exc_conductance_was_updated_local, &some_I_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_RA_neuron_bursted_local, &some_RA_neuron_bursted_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            
        
			
            ////////////////////////////////////////////
            //// Process HVC-I -> HVC-RA interactions
            ////////////////////////////////////////////
            if ( some_RA_inh_conductance_was_updated_global > 0 )
            {
            // sum update array and send to all processes

                MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                for (int i = 0; i < N_RA_local; i++)
                {
					//std::cout << "HVC-RA neuron " << Id_RA_local[i] << " raised inhibitory conductance by " << update_Gi_RA_global[Id_RA_local[i]] << std::endl;
                    HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
            	}
				// update conductance arrays and delivered indicators
				some_RA_inh_conductance_was_updated_global = 0;
				some_RA_inh_conductance_was_updated_local = 0;
				
				std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            }
            
            ////////////////////////////////////////////
            //// Process HVC-RA -> HVC-I interactions
            ////////////////////////////////////////////
			if ( some_I_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons

				for (int i = 0; i < N_I_local; i++)
					HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);

				// update conductance arrays and fired indicators
				some_I_exc_conductance_was_updated_local = 0;
				some_I_exc_conductance_was_updated_global = 0;

				std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
			}
			
			////////////////////////////////////////////
            //// Process HVC-RA -> HVC-RA interactions
            ////////////////////////////////////////////
			if ( some_RA_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons
				for (int i = 0; i < N_RA_local; i++)
					HVCRA_local[i].raiseExcWeight(update_Ge_RA_global[Id_RA_local[i]]); // update conductance

				// update conductance arrays and fired indicators
				some_RA_exc_conductance_was_updated_global = 0;
				some_RA_exc_conductance_was_updated_local = 0;
				
				std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
			}	

			////////////////////////////////////////////
            //// Process HVC-RA dendritic spikes
            ////////////////////////////////////////////
			if ( some_RA_neuron_bursted_global > 0 )
            {
				// update conductance arrays and fired indicators
            	some_RA_neuron_bursted_local = 0;
	        	some_RA_neuron_bursted_global = 0;
	        	
	        	// gather all bursted HVC-RA neurons
				this->gather_spiked_or_bursted_neurons(RA_neurons_bursted_local, RA_neurons_bursted_global);
				
				// add dendritic spikes for bursted neurons
				for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
					previous_dendritic_spike_times_global[RA_neurons_bursted_global[j]].push_back(internal_time);
				
				////////////////////////////////////////////////////////////////////////////////////
				// LTP or LTD of previously delivered somatic spike times on bursted HVC-RA neuron
				////////////////////////////////////////////////////////////////////////////////////
				
				for (int i = 0; i < N_RA_local; i++)
                {
					int presyn_ID = Id_RA_local[i]; // real id of presynaptic neuron
					
                    // if neuron is saturated apply LTP only if spiked neurons are among supersynapse targets
                    if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
                    {
                        for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_bursted_global[j]; // id of postsynaptic neuron
                            
                            // do not allow self-to-self synapses to emerge
							if ( presyn_ID != postsyn_ID )
							{
								std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
											supersynapses_local[i].end(), postsyn_ID );

								if ( pos!=supersynapses_local[i].end() )
								{
									// find previous somatic spikes that are not too old
									if ( !delivered_spike_times[i][postsyn_ID].empty() )
									{
										std::vector<double>::iterator it_relevant_spikes = std::lower_bound(delivered_spike_times[i][postsyn_ID].begin(),
																											delivered_spike_times[i][postsyn_ID].end(),
																											internal_time - STDP_WINDOW);
																							
										for (auto it = it_relevant_spikes; it != delivered_spike_times[i][postsyn_ID].end(); it++)
										{
											double dt = internal_time - *it;
									
										
											//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
											
											if (dt <= synaptic_params.T_0)
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												LTD(weights_RA_RA_local[i][postsyn_ID], dt);
												
												//~ std::cout   << "LTD from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
											}
											else
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												if ( rescaled_indicators_local[i][postsyn_ID] == 0 )
													LTP(weights_RA_RA_local[i][postsyn_ID], dt);
												else
													LTP_toRescaled(weights_RA_RA_local[i][postsyn_ID], dt);
												//~ std::cout   << "LTP from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
												
											}	
											//double w = weights_local[i][fired_ID];
											
											//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
											 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
											 //           dt, weights_local[i][fired_ID] - w);
											
											update_synapse(i, postsyn_ID);
											
										}
									}
								}
							}
                        }

                    }
                    // if not saturated apply LTP for all dendritic spikes
                    else
                    {
                        for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_bursted_global[j];
                            // don't allow self-to-self connections
                            if ( postsyn_ID != presyn_ID )
                            {
                                // find previous somatic spikes that are not too old
								if ( !delivered_spike_times[i][postsyn_ID].empty() )
								{
									std::vector<double>::iterator it_relevant_spikes = std::lower_bound(delivered_spike_times[i][postsyn_ID].begin(),
																										delivered_spike_times[i][postsyn_ID].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it = it_relevant_spikes; it != delivered_spike_times[i][postsyn_ID].end(); it++)
									{
										double dt = internal_time - *it;
								

										//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
										
									    if (dt <= synaptic_params.T_0)
									    {
											double w_before = weights_RA_RA_local[i][postsyn_ID];
												
											LTD(weights_RA_RA_local[i][postsyn_ID], dt);
											
											// show LTD on pool neurons only
											//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											//~ 
											//~ if (p.first == p.second)
											//~ {
												//~ std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											//~ }
											
											//~ std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
											
										}
										else
										{
											double w_before = weights_RA_RA_local[i][postsyn_ID];
											
											if ( rescaled_indicators_local[i][postsyn_ID] == 0 )
												LTP(weights_RA_RA_local[i][postsyn_ID], dt);
											else
												LTP_toRescaled(weights_RA_RA_local[i][postsyn_ID], dt);
										
											
											// show LTP on pool neurons only
											//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											//~ 
											//~ if (p.first == p.second)
											//~ {
												//~ std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											//~ }
											
											
											//~ std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
										}	
										//double w = weights_local[i][fired_ID];
										
										//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
										 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
										 //           dt, weights_local[i][fired_ID] - w);
										
										update_synapse(i, postsyn_ID);
									
									}
									
									
								}
                                
                            }
                        }
                   }

                } // end for i -> N_RA_local

				       

               // check if we need axon remodeling

                for (int i = 0; i < N_RA_local; i++)
                {
                    if ( (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss) && (remodeled_local[i] == 0) )
                    {
					    this->axon_remodeling(i);

				    }
                }
            
				// clear bursts
				RA_neurons_bursted_local.clear();
				RA_neurons_bursted_global.clear();
				
            } // end if some HVC-RA bursted
            
            
            network_time += NETWORK_UPDATE_FREQUENCY;
        }
    }
    
    this->potentiation_decay_sudden_maturation();
    //printf("After potentiation decay")
    //this->update_all_synapses_sudden_maturation();
    this->update_all_synapses();
	
	// update maturation info
	std::vector<int> RA_matured_local;
	std::vector<int> RA_matured_global;
	
	for (int i = 0; i < N_RA_local; i++)
	{
		num_spikes_in_recent_trials_local[i].push_front(static_cast<int>(spikes_in_trial_soma_local[i].size()));
        
        //firing_rate_short_local[i] = std::accumulate(num_spikes_in_recent_trials[i].begin(), num_spikes_in_recent_trials[i].begin() + RATE_WINDOW_SHORT, 0.0)
          //                          / static_cast<double>(RATE_WINDOW_SHORT);

		// calculate firing rate in large window:
		int maturation_window = 1000;
		double maturation_threshold = 0.7;
		
        firing_rate_long_local[i] = std::accumulate(num_spikes_in_recent_trials_local[i].begin(), num_spikes_in_recent_trials_local[i].end(), 0.0)
                                                                                            / static_cast<double>(maturation_params.RATE_WINDOW_LONG);
        if ( mature_global[Id_RA_local[i]] != 1)
        {
			
			std::vector<double> recent_firing_robustness(maturation_window);
			
			for (int j = 0; j < maturation_window; j++ )
				if ( num_spikes_in_recent_trials_local[i][j] > 0 )
					recent_firing_robustness[j] = 1;
			
			double firing_robustness = std::accumulate(recent_firing_robustness.begin(), recent_firing_robustness.end(), 0.0) / static_cast<double>(maturation_window);
			
			if ( firing_robustness > maturation_threshold )
			{
				RA_matured_local.push_back(Id_RA_local[i]);
					
				this->set_neuron_mature(i);
			}                                                                                                                          
			//if (std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].end(), 0.0) >= 1)
			//{
			//	std::cout << "Recent bursts of neuron " << Id_RA_local[i] << " :" << std::endl;
				
			//	for (size_t j = 0; j < num_bursts_in_recent_trials[i].size(); j++)
			//		std::cout << num_bursts_in_recent_trials[i][j] << "\t";
				
			//	std::cout << std::endl;
			//}
		}
	}
	
	this->gather_spiked_or_bursted_neurons(RA_matured_local, RA_matured_global);
	
	for (size_t i = 0; i < RA_matured_global.size(); i++)
	{
		int neuron_id = RA_matured_global[i];
		
		mature_global[neuron_id] = 1;
		this->rescale_synapses_to_mature(neuron_id);
	}
}

void NetworkGrowthSimulator::trial_noImmatureOut_fixedSpread(bool training, std::vector<double>& spread_times)
{
	double event_window = 10.0; // window in which all somatic spikes are considered as one event
	
	// indicators for conductance updates and bursting
	int some_RA_inh_conductance_was_updated_global = 0;
	int some_RA_inh_conductance_was_updated_local = 0;
	
	int some_RA_exc_conductance_was_updated_global = 0;
	int some_RA_exc_conductance_was_updated_local = 0;
	
	int some_I_exc_conductance_was_updated_global = 0;
	int some_I_exc_conductance_was_updated_local = 0;
	
	int some_RA_neuron_bursted_global = 0;
	int some_RA_neuron_bursted_local = 0;
	
	 // arrays with updates for conductances
    std::vector<double> update_Ge_RA_local(N_RA);
	std::vector<double> update_Ge_RA_global(N_RA);

	std::vector<double> update_Gi_RA_local(N_RA);
	std::vector<double> update_Gi_RA_global(N_RA);

	std::vector<double> update_Ge_I_local(N_I);
	std::vector<double> update_Ge_I_global(N_I);
	
	std::vector<int> RA_neurons_bursted_local; // local array of HVC-RA neurons that bursted
	std::vector<int> RA_neurons_bursted_global; // global array of HVC-RA neurons that bursted
	
	// sample training innervation time and send to all processes
  // training neurons innervation
	std::vector<bool> indicators_training(N_RA);
	std::fill(indicators_training.begin(), indicators_training.end(), false);
	
	for (int i = 0; i < N_TR; i++)
		indicators_training[training_neurons[i]] = true;
		
	std::vector<bool> indicators_current_injected(N_RA);
	std::fill(indicators_current_injected.begin(), indicators_current_injected.end(), false);
	
	// sample mean training innervation time and send to all processes
    double mean_injection_time;
    
	if (MPI_rank == 0)
	{
		mean_injection_time = WAITING_TIME + noise_generator.random(TRIAL_DURATION - 2*WAITING_TIME);
		std::cout << "mean_injection_time = " << mean_injection_time << std::endl;
    }
    MPI_Bcast(&mean_injection_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
    double internal_time = 0;
    double network_time = NETWORK_UPDATE_FREQUENCY;

   

    //printf("network time = %f\n", network_time);

    
	// initialize update arrays and fired indicators
	std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
	
	
    // evolve dynamics
    int num_steps = static_cast<int>(TRIAL_DURATION / TIMESTEP);
    
    for (int t = 1; t < num_steps; t++)
	{
		internal_time += TIMESTEP;
		
		//////////////////////
		// HVC-RA neurons
		//////////////////////
		
		for (int i = 0; i < N_RA_local; i++)
		{
			if ( ( indicators_training[Id_RA_local[i]] ) && ( !indicators_current_injected[Id_RA_local[i]] ) )
			{
				if ( internal_time > spread_times[Id_RA_local[i]] + mean_injection_time )
				{
					HVCRA_local[i].raiseE(G_TRAINING_KICK);
					indicators_current_injected[Id_RA_local[i]] = true;
				}
			}
			
            // Debraband step
            HVCRA_local[i].Debraband_step_no_target_update();
            
           
			///////////////////////////////
			// HVC-RA somatic spike
			///////////////////////////////
			// if some neuron produced somatic spike
            if (HVCRA_local[i].get_fired_soma())
            {
                spikes_in_trial_soma_local[i].push_back(internal_time);
                
				// update delivery queues if neuron is mature (no output for immature neuron and no wiring mechanism)
				if ( mature_global[Id_RA_local[i]] == 1 ){
				
					// for inhibitory neurons
					// loop over all inhibitory targets of fired neurons
					size_t num_I_targets = syn_ID_RA_I_local[i].size();
					
					for (size_t j = 0; j < num_I_targets; j++)
					{
						double delivery_time = internal_time + axonal_delays_RA_I_local[i][j];
						
						// if queue is empty, just add item to the queue
						if ( delivery_queue_RA_I[i].empty() )
							delivery_queue_RA_I[i].push_back(std::pair<double,int>(delivery_time, j));	
						// otherwise add item so that queue is sorted
						else
						{
							auto it = std::upper_bound(delivery_queue_RA_I[i].begin(), delivery_queue_RA_I[i].end(), std::pair<double,int>(delivery_time, j));
							delivery_queue_RA_I[i].insert(it, std::pair<double,int>(delivery_time, j));
						}
					}
						//std::cout << "neuron " << Id_RA_local[i] << " spike time = " << internal_time << " axonal delay = " << axonal_delays_RA_I[Id_RA_local[i]][syn_ID] << " delivery time RA to I: " << delivery_queue_RA_I[i].back().first << " delivery target id: " << syn_ID_RA_I_local[i][delivery_queue_RA_I[i].back().second] << std::endl;
					
					
					// if neuron is saturated, deliver spikes only to supersynaptic targets
					if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
					{
						for (size_t j = 0; j < supersynapses_local[i].size(); j++)
						{
							int target_id = supersynapses_local[i][j];
							
							double delivery_time = internal_time + axonal_delays_RA_RA_local[i][target_id];
							
							// if queue is empty, just add item to the queue
							if ( delivery_queue_RA_RA_soma[i].empty() )
								delivery_queue_RA_RA_soma[i].push_back(std::pair<double,int>(delivery_time, target_id));	
							// otherwise add item so that queue is sorted
							else
							{
								auto it = std::upper_bound(delivery_queue_RA_RA_soma[i].begin(), delivery_queue_RA_RA_soma[i].end(), std::pair<double,int>(delivery_time, target_id));
								delivery_queue_RA_RA_soma[i].insert(it, std::pair<double,int>(delivery_time, target_id));
							}
						}
					}
					else // deliver spikes to everyone except itself
					{	
						for (int j = 0; j < N_RA; j++)
						{
							if ( j != Id_RA_local[i])
							{
							
								double delivery_time = internal_time + axonal_delays_RA_RA_local[i][j];
							
								// if queue is empty, just add item to the queue
								if ( delivery_queue_RA_RA_soma[i].empty() )
									delivery_queue_RA_RA_soma[i].push_back(std::pair<double,int>(delivery_time, j));	
								// otherwise add item so that queue is sorted
								else
								{
									auto it = std::upper_bound(delivery_queue_RA_RA_soma[i].begin(), delivery_queue_RA_RA_soma[i].end(), std::pair<double,int>(delivery_time, j));
									delivery_queue_RA_RA_soma[i].insert(it, std::pair<double,int>(delivery_time, j));
								}
							}
						}
					}
				}
                
            } // end if local HVC-RA neuron spiked  
            
             //////////////////////////////////
			// HVC-RA -> HVC-I delivery queue
			//////////////////////////////////
            // check if somatic spike was delivered to some interneuron
            // loop through the delivery queue to check if current time exceeds the spike delivery time
            std::vector<std::pair<double,int>>::iterator it = delivery_queue_RA_I[i].begin();
			
            for (; it != delivery_queue_RA_I[i].end(); it++)
            {
				if (internal_time >= it->first)
				{
					
					some_I_exc_conductance_was_updated_local = 1;
					
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_RA_I_local[i][pos_in_local_target_array];
					
					update_Ge_I_local[target_id] += weights_RA_I_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(RA) neuron " << Id_RA_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to HVC(I) neuron " << syn_ID_RA_I_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_I[i].erase(delivery_queue_RA_I[i].begin(), it);
			
			//////////////////////////////////////////////////
			// HVC-RA -> HVC-RA somatic spike delivery queue
			//////////////////////////////////////////////////
			// check if somatic spike was delivered to some HVC(RA) neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
            it = delivery_queue_RA_RA_soma[i].begin();
            
            for (; it != delivery_queue_RA_RA_soma[i].end(); it++)
            {
				double delivery_time = it->first;
				
				if ( internal_time >= delivery_time )
				{
					int target_id = it->second;
					
					
					
					// update conductance of target if synapse is active
					if ( active_indicators_local[i][target_id] == 1 )
					{
						some_RA_exc_conductance_was_updated_local = 1;
						update_Ge_RA_local[target_id] += weights_RA_RA_local[i][target_id];
					}
					
					
					
					//////////////
					//// LTD
					//////////////
					// do LTD on dendritic spike time of target neuron if it is not mature neuron
					//~ if ( mature_global[target_id] != 1 )
					//~ {
						// update delivered spikes only if current spike occured later than previous delivered
						// spike + event_window
						bool new_event_occured = false;
					
						if ( !delivered_spike_times[i][target_id].empty() ){
							if ( delivery_time > delivered_spike_times[i][target_id].back() + event_window ){
								delivered_spike_times[i][target_id].push_back(delivery_time);
								new_event_occured = true;
							}
						}
						else{
							delivered_spike_times[i][target_id].push_back(delivery_time);
							new_event_occured = true;
						}
						
						if ( new_event_occured ){
							 // if neuron is saturated apply LTD only if target is among super synapses
							if (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss)
							{
								if ( supersynapses_indicators_local[i][target_id] == 1 )
								{	
									// find previous dendritic spikes that are not too old
									if ( !previous_dendritic_spike_times_global[target_id].empty() )
									{
										std::vector<double>::iterator it_relevant_dend_spikes = std::lower_bound(previous_dendritic_spike_times_global[target_id].begin(),
																											previous_dendritic_spike_times_global[target_id].end(),
																											internal_time - STDP_WINDOW);
																							
										for (auto it_dend_spike = it_relevant_dend_spikes; it_dend_spike != previous_dendritic_spike_times_global[target_id].end(); it_dend_spike++)
										{
											double dt = *it_dend_spike - delivery_time;

											//std::cout << "dt in saturated LTD = " << dt << std::endl;

											
											double w_before = weights_RA_RA_local[i][target_id];     
								   
											LTD(weights_RA_RA_local[i][target_id], dt);
											
											
											//~ std::cout << "LTD from saturated " << Id_RA_local[i] << " -> " << supersynapse_id
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][supersynapse_id] - w_before
													  //~ << std::endl;
													  
											//printf("LTD from saturated neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n",
											//            Id_RA_local[i], supersynapse_id, spikes_in_trial_soma_local[i].back(), spike_times_dend_global[supersynapse_id],
											//            dt, weights_local[i][supersynapse_id] - w);
											
											update_synapse(i, target_id);	
											
										}   
									}
								}

								// if some supersynapse desaturated, update all synapses
								if (static_cast<int>(supersynapses_local[i].size()) < synaptic_params.Nss)
									for (int j = 0; j < N_RA; j++)
										this->update_synapse(i, j);
							}
							// if not saturated apply LTD rule 
							else
							{
								// find previous dendritic spikes that are not too old
								if ( !previous_dendritic_spike_times_global[target_id].empty() )
								{
									std::vector<double>::iterator it_relevant_dend_spikes = std::lower_bound(previous_dendritic_spike_times_global[target_id].begin(),
																										previous_dendritic_spike_times_global[target_id].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it_dend_spike = it_relevant_dend_spikes; it_dend_spike != previous_dendritic_spike_times_global[target_id].end(); it_dend_spike++)
									{
										double dt = *it_dend_spike - delivery_time;

										//std::cout << "dt in saturated LTD = " << dt << std::endl;

										
										double w_before = weights_RA_RA_local[i][target_id];     
							   
										LTD(weights_RA_RA_local[i][target_id], dt);
										
										//~ // show LTD to pool neurons only
										//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), target_id);
										//~ 
										//~ if (p.first == p.second)
										//~ {
											//~ std::cout << "LTD from  " << Id_RA_local[i] << " -> " << target_id
												  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][target_id] - w_before
												  //~ << " w_after = " << weights_RA_RA_local[i][target_id] << std::endl;
										//~ }
									//~ 
										update_synapse(i, target_id);	
										
									}   
								}
							}
						}
					//~ } // end if target neuron is not mature
				} // end if internal_time >= delivery spike time 
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_RA_soma[i].erase(delivery_queue_RA_RA_soma[i].begin(), it);
			
			
			///////////////////////////////
			// HVC-RA dendritic spike
			///////////////////////////////
			//if some neuron produced dendritic spike, store this neuron in array
			
			if ( HVCRA_local[i].get_fired_dend() )
			{
				//~ // show spikes of pool neurons only
				//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
				//~ 
				//~ if (p.first == p.second)
					//~ std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
				//~ 
				std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
				
				spikes_in_trial_dend_local[i].push_back(internal_time); // dend spike time relative to trial onset
				
				//if ( mature_global[Id_RA_local[i]] != 1 )
				//{
					RA_neurons_bursted_local.push_back(Id_RA_local[i]);
				
					some_RA_neuron_bursted_local = 1;
				//}
				//std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
			} // end if get_fired_dend
	
            
		} // end for i = 0 -> N_RA_local (loop through all local HVC-RA neurons)
		
		//////////////////////
		// HVC-I neurons
		//////////////////////
		for (int i = 0; i < N_I_local; i++)
		{
			HVCI_local[i].DP8_step_no_target_update();
			
			
			///////////////////////////////
			// HVC-I spike
			///////////////////////////////
			//  if some I neuron spikes, update delivery queue
			if (HVCI_local[i].get_fired())
			{
				//printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
				spikes_in_trial_interneuron_local[i].push_back(internal_time);

				size_t num_RA_targets = syn_ID_I_RA_local[i].size();
				// loop over all targets of fired neurons
				for (size_t j = 0; j < num_RA_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_I_RA_local[i][j];
					
					//std::cout << "HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << syn_ID_I_RA_local[i][j] << " spike time " << internal_time << " delivery time = " << delivery_time << std::endl; 
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_I_RA[i].empty() )
						delivery_queue_I_RA[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_I_RA[i].begin(), delivery_queue_I_RA[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_I_RA[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
			}
			
			//////////////////////////////////
			// HVC-I -> HVC-RA delivery queue
			//////////////////////////////////
			// check if interneuron spike was delivered to some HVC-RA neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
			std::vector<std::pair<double,int>>::iterator it = delivery_queue_I_RA[i].begin();
			
			for (; it != delivery_queue_I_RA[i].end(); it++)
			{
				if (internal_time >= it->first)
				{
					some_RA_inh_conductance_was_updated_local = 1;
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_I_RA_local[i][pos_in_local_target_array];
					
					//std::cout << "Delivered spike HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << target_id << " delivered_time " << internal_time << " delivery time in queue = " << it->first << std::endl; 
					
					
					update_Gi_RA_local[target_id] += weights_I_RA_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(I) neuron " << Id_I_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to neuron " << syn_ID_I_RA_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_I_RA[i].erase(delivery_queue_I_RA[i].begin(), it);
		
		}  // end for i = 0 -> N_I_local (loop through all local HVC-I neurons)
		
		/////////////////////////////
		// Network Synchronization //
		/////////////////////////////
		if ( ( internal_time > network_time ) || ( t == num_steps-1 ) )
        {
            //if (MPI_rank == 0)
            //{
            //    std::cout << "internal time = " << internal_time << " > network_time = " << network_time << std::endl; 
            //}

            MPI_Allreduce(&some_RA_exc_conductance_was_updated_local, &some_RA_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_RA_inh_conductance_was_updated_local, &some_RA_inh_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_I_exc_conductance_was_updated_local, &some_I_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_RA_neuron_bursted_local, &some_RA_neuron_bursted_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            
        
			
            ////////////////////////////////////////////
            //// Process HVC-I -> HVC-RA interactions
            ////////////////////////////////////////////
            if ( some_RA_inh_conductance_was_updated_global > 0 )
            {
            // sum update array and send to all processes

                MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                for (int i = 0; i < N_RA_local; i++)
                {
					//std::cout << "HVC-RA neuron " << Id_RA_local[i] << " raised inhibitory conductance by " << update_Gi_RA_global[Id_RA_local[i]] << std::endl;
                    HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
            	}
				// update conductance arrays and delivered indicators
				some_RA_inh_conductance_was_updated_global = 0;
				some_RA_inh_conductance_was_updated_local = 0;
				
				std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            }
            
            ////////////////////////////////////////////
            //// Process HVC-RA -> HVC-I interactions
            ////////////////////////////////////////////
			if ( some_I_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons

				for (int i = 0; i < N_I_local; i++)
					HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);

				// update conductance arrays and fired indicators
				some_I_exc_conductance_was_updated_local = 0;
				some_I_exc_conductance_was_updated_global = 0;

				std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
			}
			
			////////////////////////////////////////////
            //// Process HVC-RA -> HVC-RA interactions
            ////////////////////////////////////////////
			if ( some_RA_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons
				for (int i = 0; i < N_RA_local; i++)
					HVCRA_local[i].raiseExcWeight(update_Ge_RA_global[Id_RA_local[i]]); // update conductance

				// update conductance arrays and fired indicators
				some_RA_exc_conductance_was_updated_global = 0;
				some_RA_exc_conductance_was_updated_local = 0;
				
				std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
			}	

			////////////////////////////////////////////
            //// Process HVC-RA dendritic spikes
            ////////////////////////////////////////////
			if ( some_RA_neuron_bursted_global > 0 )
            {
				// update conductance arrays and fired indicators
            	some_RA_neuron_bursted_local = 0;
	        	some_RA_neuron_bursted_global = 0;
	        	
	        	// gather all bursted HVC-RA neurons
				this->gather_spiked_or_bursted_neurons(RA_neurons_bursted_local, RA_neurons_bursted_global);
				
				// add dendritic spikes for bursted neurons
				for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
					previous_dendritic_spike_times_global[RA_neurons_bursted_global[j]].push_back(internal_time);
				
				////////////////////////////////////////////////////////////////////////////////////
				// LTP or LTD of previously delivered somatic spike times on bursted HVC-RA neuron
				////////////////////////////////////////////////////////////////////////////////////
				
				for (int i = 0; i < N_RA_local; i++)
                {
					int presyn_ID = Id_RA_local[i]; // real id of presynaptic neuron
					
                    // if neuron is saturated apply LTP only if spiked neurons are among supersynapse targets
                    if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
                    {
                        for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_bursted_global[j]; // id of postsynaptic neuron
                            
                            // do not allow self-to-self synapses to emerge
							if ( presyn_ID != postsyn_ID )
							{
								std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
											supersynapses_local[i].end(), postsyn_ID );

								if ( pos!=supersynapses_local[i].end() )
								{
									// find previous somatic spikes that are not too old
									if ( !delivered_spike_times[i][postsyn_ID].empty() )
									{
										std::vector<double>::iterator it_relevant_spikes = std::lower_bound(delivered_spike_times[i][postsyn_ID].begin(),
																											delivered_spike_times[i][postsyn_ID].end(),
																											internal_time - STDP_WINDOW);
																							
										for (auto it = it_relevant_spikes; it != delivered_spike_times[i][postsyn_ID].end(); it++)
										{
											double dt = internal_time - *it;
									
										
											//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
											
											if (dt <= synaptic_params.T_0)
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												LTD(weights_RA_RA_local[i][postsyn_ID], dt);
												
												//~ std::cout   << "LTD from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
											}
											else
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												if ( rescaled_indicators_local[i][postsyn_ID] == 0 )
													LTP(weights_RA_RA_local[i][postsyn_ID], dt);
												else
													LTP_toRescaled(weights_RA_RA_local[i][postsyn_ID], dt);
												//~ std::cout   << "LTP from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
												
											}	
											//double w = weights_local[i][fired_ID];
											
											//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
											 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
											 //           dt, weights_local[i][fired_ID] - w);
											
											update_synapse(i, postsyn_ID);
											
										}
									}
								}
							}
                        }

                    }
                    // if not saturated apply LTP for all dendritic spikes
                    else
                    {
                        for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_bursted_global[j];
                            // don't allow self-to-self connections
                            if ( postsyn_ID != presyn_ID )
                            {
                                // find previous somatic spikes that are not too old
								if ( !delivered_spike_times[i][postsyn_ID].empty() )
								{
									std::vector<double>::iterator it_relevant_spikes = std::lower_bound(delivered_spike_times[i][postsyn_ID].begin(),
																										delivered_spike_times[i][postsyn_ID].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it = it_relevant_spikes; it != delivered_spike_times[i][postsyn_ID].end(); it++)
									{
										double dt = internal_time - *it;
								

										//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
										
									    if (dt <= synaptic_params.T_0)
									    {
											double w_before = weights_RA_RA_local[i][postsyn_ID];
												
											LTD(weights_RA_RA_local[i][postsyn_ID], dt);
											
											// show LTD on pool neurons only
											//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											//~ 
											//~ if (p.first == p.second)
											//~ {
												//~ std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											//~ }
											
											//~ std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
											
										}
										else
										{
											double w_before = weights_RA_RA_local[i][postsyn_ID];
											
											if ( rescaled_indicators_local[i][postsyn_ID] == 0 )
												LTP(weights_RA_RA_local[i][postsyn_ID], dt);
											else
												LTP_toRescaled(weights_RA_RA_local[i][postsyn_ID], dt);
										
											
											// show LTP on pool neurons only
											//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											//~ 
											//~ if (p.first == p.second)
											//~ {
												//~ std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											//~ }
											
											
											//~ std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
										}	
										//double w = weights_local[i][fired_ID];
										
										//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
										 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
										 //           dt, weights_local[i][fired_ID] - w);
										
										update_synapse(i, postsyn_ID);
									
									}
									
									
								}
                                
                            }
                        }
                   }

                } // end for i -> N_RA_local

				       

               // check if we need axon remodeling

                for (int i = 0; i < N_RA_local; i++)
                {
                    if ( (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss) && (remodeled_local[i] == 0) )
                    {
					    this->axon_remodeling(i);

				    }
                }
            
				// clear bursts
				RA_neurons_bursted_local.clear();
				RA_neurons_bursted_global.clear();
				
            } // end if some HVC-RA bursted
            
            
            network_time += NETWORK_UPDATE_FREQUENCY;
        }
    }
    
    this->potentiation_decay_sudden_maturation();
    //printf("After potentiation decay")
    //this->update_all_synapses_sudden_maturation();
    this->update_all_synapses();
	
	// update maturation info
	std::vector<int> RA_matured_local;
	std::vector<int> RA_matured_global;
	
	for (int i = 0; i < N_RA_local; i++)
	{
		num_spikes_in_recent_trials_local[i].push_front(static_cast<int>(spikes_in_trial_soma_local[i].size()));
        
        //firing_rate_short_local[i] = std::accumulate(num_spikes_in_recent_trials[i].begin(), num_spikes_in_recent_trials[i].begin() + RATE_WINDOW_SHORT, 0.0)
          //                          / static_cast<double>(RATE_WINDOW_SHORT);

		// calculate firing rate in large window:
		int maturation_window = 300;
		double maturation_threshold = 0.7;
		
        firing_rate_long_local[i] = std::accumulate(num_spikes_in_recent_trials_local[i].begin(), num_spikes_in_recent_trials_local[i].end(), 0.0)
                                                                                            / static_cast<double>(maturation_params.RATE_WINDOW_LONG);
        if ( mature_global[Id_RA_local[i]] != 1)
        {
			
			std::vector<double> recent_firing_robustness(maturation_window);
			
			for (int j = 0; j < maturation_window; j++ )
				if ( num_spikes_in_recent_trials_local[i][j] > 0 )
					recent_firing_robustness[j] = 1;
			
			double firing_robustness = std::accumulate(recent_firing_robustness.begin(), recent_firing_robustness.end(), 0.0) / static_cast<double>(maturation_window);
			
			if ( firing_robustness > maturation_threshold )
			{
				RA_matured_local.push_back(Id_RA_local[i]);
					
				this->set_neuron_mature(i);
			}                                                                                                                          
			//if (std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].end(), 0.0) >= 1)
			//{
			//	std::cout << "Recent bursts of neuron " << Id_RA_local[i] << " :" << std::endl;
				
			//	for (size_t j = 0; j < num_bursts_in_recent_trials[i].size(); j++)
			//		std::cout << num_bursts_in_recent_trials[i][j] << "\t";
				
			//	std::cout << std::endl;
			//}
		}
	}
	
	this->gather_spiked_or_bursted_neurons(RA_matured_local, RA_matured_global);
	
	for (size_t i = 0; i < RA_matured_global.size(); i++)
	{
		int neuron_id = RA_matured_global[i];
		
		mature_global[neuron_id] = 1;
		this->rescale_synapses_to_mature(neuron_id);
	}
}

void NetworkGrowthSimulator::trial_event_pre_dend_post_delays_sudden_maturation(bool training)
{
	double event_window = 10.0; // window in which all somatic spikes are considered as one event
	
	// indicators for conductance updates and bursting
	int some_RA_inh_conductance_was_updated_global = 0;
	int some_RA_inh_conductance_was_updated_local = 0;
	
	int some_RA_exc_conductance_was_updated_global = 0;
	int some_RA_exc_conductance_was_updated_local = 0;
	
	int some_I_exc_conductance_was_updated_global = 0;
	int some_I_exc_conductance_was_updated_local = 0;
	
	int some_RA_neuron_bursted_global = 0;
	int some_RA_neuron_bursted_local = 0;
	
	 // arrays with updates for conductances
    std::vector<double> update_Ge_RA_local(N_RA);
	std::vector<double> update_Ge_RA_global(N_RA);

	std::vector<double> update_Gi_RA_local(N_RA);
	std::vector<double> update_Gi_RA_global(N_RA);

	std::vector<double> update_Ge_I_local(N_I);
	std::vector<double> update_Ge_I_global(N_I);
	
	std::vector<int> RA_neurons_bursted_local; // local array of HVC-RA neurons that bursted
	std::vector<int> RA_neurons_bursted_global; // global array of HVC-RA neurons that bursted
	
	// sample training innervation time and send to all processes
    double training_kick_time;
    
	if (MPI_rank == 0)
	{
		training_kick_time = WAITING_TIME + noise_generator.random(TRIAL_DURATION - 2*WAITING_TIME);
		std::cout << "training_kick_time = " << training_kick_time << std::endl;
    }
    MPI_Bcast(&training_kick_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    bool training_excited = false; // indicator that training neurons were already excited
    
    double internal_time = 0;
    double network_time = NETWORK_UPDATE_FREQUENCY;

   

    //printf("network time = %f\n", network_time);

    
	// initialize update arrays and fired indicators
	std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
	
	
    // evolve dynamics
    int num_steps = static_cast<int>(TRIAL_DURATION / TIMESTEP);
    
    for (int t = 1; t < num_steps; t++)
	{
		internal_time += TIMESTEP;
	
		if ( ( training ) && ( !training_excited ) && ( internal_time >= training_kick_time ) )
		{
			for (int i = 0; i < N_TR; i++)
			{
				int rank;
				int shift;
				
				this->get_neuronRA_location(training_neurons[i], &rank, &shift);
				
				if (MPI_rank == rank)
					HVCRA_local[shift].raiseE(G_TRAINING_KICK);
			}
			
			training_excited = true;
		}
		
		//////////////////////
		// HVC-RA neurons
		//////////////////////
		
		for (int i = 0; i < N_RA_local; i++)
		{
            // Debraband step
            HVCRA_local[i].Debraband_step_no_target_update();
            
           
			///////////////////////////////
			// HVC-RA somatic spike
			///////////////////////////////
			// if some neuron produced somatic spike
            if (HVCRA_local[i].get_fired_soma())
            {
                spikes_in_trial_soma_local[i].push_back(internal_time);
                
				// update delivery queues
				
				// for inhibitory neurons
				// loop over all inhibitory targets of fired neurons
				size_t num_I_targets = syn_ID_RA_I_local[i].size();
				
				for (size_t j = 0; j < num_I_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_RA_I_local[i][j];
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_RA_I[i].empty() )
						delivery_queue_RA_I[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_RA_I[i].begin(), delivery_queue_RA_I[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_RA_I[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
					//std::cout << "neuron " << Id_RA_local[i] << " spike time = " << internal_time << " axonal delay = " << axonal_delays_RA_I[Id_RA_local[i]][syn_ID] << " delivery time RA to I: " << delivery_queue_RA_I[i].back().first << " delivery target id: " << syn_ID_RA_I_local[i][delivery_queue_RA_I[i].back().second] << std::endl;
				
				
				// if neuron is saturated, deliver spikes only to supersynaptic targets
				if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
				{
					for (size_t j = 0; j < supersynapses_local[i].size(); j++)
					{
						int target_id = supersynapses_local[i][j];
						
						double delivery_time = internal_time + axonal_delays_RA_RA_local[i][target_id];
						
						// if queue is empty, just add item to the queue
						if ( delivery_queue_RA_RA_soma[i].empty() )
							delivery_queue_RA_RA_soma[i].push_back(std::pair<double,int>(delivery_time, target_id));	
						// otherwise add item so that queue is sorted
						else
						{
							auto it = std::upper_bound(delivery_queue_RA_RA_soma[i].begin(), delivery_queue_RA_RA_soma[i].end(), std::pair<double,int>(delivery_time, target_id));
							delivery_queue_RA_RA_soma[i].insert(it, std::pair<double,int>(delivery_time, target_id));
						}
					}
				}
				else // deliver spikes to everyone except itself
				{	
					for (int j = 0; j < N_RA; j++)
					{
						if ( j != Id_RA_local[i])
						{
						
							double delivery_time = internal_time + axonal_delays_RA_RA_local[i][j];
						
							// if queue is empty, just add item to the queue
							if ( delivery_queue_RA_RA_soma[i].empty() )
								delivery_queue_RA_RA_soma[i].push_back(std::pair<double,int>(delivery_time, j));	
							// otherwise add item so that queue is sorted
							else
							{
								auto it = std::upper_bound(delivery_queue_RA_RA_soma[i].begin(), delivery_queue_RA_RA_soma[i].end(), std::pair<double,int>(delivery_time, j));
								delivery_queue_RA_RA_soma[i].insert(it, std::pair<double,int>(delivery_time, j));
							}
						}
					}
				}
                
            } // end if local HVC-RA neuron spiked  
            
             //////////////////////////////////
			// HVC-RA -> HVC-I delivery queue
			//////////////////////////////////
            // check if somatic spike was delivered to some interneuron
            // loop through the delivery queue to check if current time exceeds the spike delivery time
            std::vector<std::pair<double,int>>::iterator it = delivery_queue_RA_I[i].begin();
			
            for (; it != delivery_queue_RA_I[i].end(); it++)
            {
				if (internal_time >= it->first)
				{
					
					some_I_exc_conductance_was_updated_local = 1;
					
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_RA_I_local[i][pos_in_local_target_array];
					
					update_Ge_I_local[target_id] += weights_RA_I_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(RA) neuron " << Id_RA_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to HVC(I) neuron " << syn_ID_RA_I_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_I[i].erase(delivery_queue_RA_I[i].begin(), it);
			
			//////////////////////////////////////////////////
			// HVC-RA -> HVC-RA somatic spike delivery queue
			//////////////////////////////////////////////////
			// check if somatic spike was delivered to some HVC(RA) neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
            it = delivery_queue_RA_RA_soma[i].begin();
            
            for (; it != delivery_queue_RA_RA_soma[i].end(); it++)
            {
				double delivery_time = it->first;
				
				if ( internal_time >= delivery_time )
				{
					int target_id = it->second;
					
					
					
					// update conductance of target if synapse is active
					if ( active_indicators_local[i][target_id] == 1 )
					{
						some_RA_exc_conductance_was_updated_local = 1;
						update_Ge_RA_local[target_id] += weights_RA_RA_local[i][target_id];
					}
					
					
					
					//////////////
					//// LTD
					//////////////
					// do LTD on dendritic spike time of target neuron if it is not mature neuron
					if ( mature_global[target_id] != 1 )
					{
						// update delivered spikes only if current spike occured later than previous delivered
						// spike + event_window
						bool new_event_occured = false;
					
						if ( !delivered_spike_times[i][target_id].empty() ){
							if ( delivery_time > delivered_spike_times[i][target_id].back() + event_window ){
								delivered_spike_times[i][target_id].push_back(delivery_time);
								new_event_occured = true;
							}
						}
						else{
							delivered_spike_times[i][target_id].push_back(delivery_time);
							new_event_occured = true;
						}
						
						if ( new_event_occured ){
							 // if neuron is saturated apply LTD only if target is among super synapses
							if (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss)
							{
								if ( supersynapses_indicators_local[i][target_id] == 1 )
								{	
									// find previous dendritic spikes that are not too old
									if ( !previous_dendritic_spike_times_global[target_id].empty() )
									{
										std::vector<double>::iterator it_relevant_dend_spikes = std::lower_bound(previous_dendritic_spike_times_global[target_id].begin(),
																											previous_dendritic_spike_times_global[target_id].end(),
																											internal_time - STDP_WINDOW);
																							
										for (auto it_dend_spike = it_relevant_dend_spikes; it_dend_spike != previous_dendritic_spike_times_global[target_id].end(); it_dend_spike++)
										{
											double dt = *it_dend_spike - delivery_time;

											//std::cout << "dt in saturated LTD = " << dt << std::endl;

											
											double w_before = weights_RA_RA_local[i][target_id];     
								   
											LTD(weights_RA_RA_local[i][target_id], dt);
											
											
											//~ std::cout << "LTD from saturated " << Id_RA_local[i] << " -> " << supersynapse_id
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][supersynapse_id] - w_before
													  //~ << std::endl;
													  
											//printf("LTD from saturated neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n",
											//            Id_RA_local[i], supersynapse_id, spikes_in_trial_soma_local[i].back(), spike_times_dend_global[supersynapse_id],
											//            dt, weights_local[i][supersynapse_id] - w);
											
											update_synapse(i, target_id);	
											
										}   
									}
								}

								// if some supersynapse desaturated, update all synapses
								if (static_cast<int>(supersynapses_local[i].size()) < synaptic_params.Nss)
									for (int j = 0; j < N_RA; j++)
										this->update_synapse(i, j);
							}
							// if not saturated apply LTD rule 
							else
							{
								// find previous dendritic spikes that are not too old
								if ( !previous_dendritic_spike_times_global[target_id].empty() )
								{
									std::vector<double>::iterator it_relevant_dend_spikes = std::lower_bound(previous_dendritic_spike_times_global[target_id].begin(),
																										previous_dendritic_spike_times_global[target_id].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it_dend_spike = it_relevant_dend_spikes; it_dend_spike != previous_dendritic_spike_times_global[target_id].end(); it_dend_spike++)
									{
										double dt = *it_dend_spike - delivery_time;

										//std::cout << "dt in saturated LTD = " << dt << std::endl;

										
										double w_before = weights_RA_RA_local[i][target_id];     
							   
										LTD(weights_RA_RA_local[i][target_id], dt);
										
										// show LTD of pool neurons only
										//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
										//~ 
										//~ if (p.first == p.second)
										//~ {
											//~ std::cout << "LTD from  " << Id_RA_local[i] << " -> " << target_id
												  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][target_id] - w_before
												  //~ << " w_after = " << weights_RA_RA_local[i][target_id] << std::endl;
										//~ }
									
										update_synapse(i, target_id);	
										
									}   
								}
							}
						}
					} // end if target neuron is not mature
				} // end if internal_time >= delivery spike time 
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_RA_soma[i].erase(delivery_queue_RA_RA_soma[i].begin(), it);
			
			
			///////////////////////////////
			// HVC-RA dendritic spike
			///////////////////////////////
			//if some neuron produced dendritic spike, store this neuron in array
			
			if ( HVCRA_local[i].get_fired_dend() )
			{
				// show spikes of pool neurons only
				auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
				
				if (p.first == p.second)
					std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
				
				
				spikes_in_trial_dend_local[i].push_back(internal_time); // dend spike time relative to trial onset
				
				if ( mature_global[Id_RA_local[i]] != 1 )
				{
					RA_neurons_bursted_local.push_back(Id_RA_local[i]);
				
					some_RA_neuron_bursted_local = 1;
				}
				//std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
			} // end if get_fired_dend
	
            
		} // end for i = 0 -> N_RA_local (loop through all local HVC-RA neurons)
		
		//////////////////////
		// HVC-I neurons
		//////////////////////
		for (int i = 0; i < N_I_local; i++)
		{
			HVCI_local[i].DP8_step_no_target_update();
			
			
			///////////////////////////////
			// HVC-I spike
			///////////////////////////////
			//  if some I neuron spikes, update delivery queue
			if (HVCI_local[i].get_fired())
			{
				//printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
				spikes_in_trial_interneuron_local[i].push_back(internal_time);

				size_t num_RA_targets = syn_ID_I_RA_local[i].size();
				// loop over all targets of fired neurons
				for (size_t j = 0; j < num_RA_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_I_RA_local[i][j];
					
					//std::cout << "HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << syn_ID_I_RA_local[i][j] << " spike time " << internal_time << " delivery time = " << delivery_time << std::endl; 
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_I_RA[i].empty() )
						delivery_queue_I_RA[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_I_RA[i].begin(), delivery_queue_I_RA[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_I_RA[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
			}
			
			//////////////////////////////////
			// HVC-I -> HVC-RA delivery queue
			//////////////////////////////////
			// check if interneuron spike was delivered to some HVC-RA neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
			std::vector<std::pair<double,int>>::iterator it = delivery_queue_I_RA[i].begin();
			
			for (; it != delivery_queue_I_RA[i].end(); it++)
			{
				if (internal_time >= it->first)
				{
					some_RA_inh_conductance_was_updated_local = 1;
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_I_RA_local[i][pos_in_local_target_array];
					
					//std::cout << "Delivered spike HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << target_id << " delivered_time " << internal_time << " delivery time in queue = " << it->first << std::endl; 
					
					
					update_Gi_RA_local[target_id] += weights_I_RA_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(I) neuron " << Id_I_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to neuron " << syn_ID_I_RA_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_I_RA[i].erase(delivery_queue_I_RA[i].begin(), it);
		
		}  // end for i = 0 -> N_I_local (loop through all local HVC-I neurons)
		
		/////////////////////////////
		// Network Synchronization //
		/////////////////////////////
		if ( ( internal_time > network_time ) || ( t == num_steps-1 ) )
        {
            //if (MPI_rank == 0)
            //{
            //    std::cout << "internal time = " << internal_time << " > network_time = " << network_time << std::endl; 
            //}

            MPI_Allreduce(&some_RA_exc_conductance_was_updated_local, &some_RA_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_RA_inh_conductance_was_updated_local, &some_RA_inh_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_I_exc_conductance_was_updated_local, &some_I_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_RA_neuron_bursted_local, &some_RA_neuron_bursted_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            
        
			
            ////////////////////////////////////////////
            //// Process HVC-I -> HVC-RA interactions
            ////////////////////////////////////////////
            if ( some_RA_inh_conductance_was_updated_global > 0 )
            {
            // sum update array and send to all processes

                MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                for (int i = 0; i < N_RA_local; i++)
                {
					//std::cout << "HVC-RA neuron " << Id_RA_local[i] << " raised inhibitory conductance by " << update_Gi_RA_global[Id_RA_local[i]] << std::endl;
                    HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
            	}
				// update conductance arrays and delivered indicators
				some_RA_inh_conductance_was_updated_global = 0;
				some_RA_inh_conductance_was_updated_local = 0;
				
				std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            }
            
            ////////////////////////////////////////////
            //// Process HVC-RA -> HVC-I interactions
            ////////////////////////////////////////////
			if ( some_I_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons

				for (int i = 0; i < N_I_local; i++)
					HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);

				// update conductance arrays and fired indicators
				some_I_exc_conductance_was_updated_local = 0;
				some_I_exc_conductance_was_updated_global = 0;

				std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
			}
			
			////////////////////////////////////////////
            //// Process HVC-RA -> HVC-RA interactions
            ////////////////////////////////////////////
			if ( some_RA_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons
				for (int i = 0; i < N_RA_local; i++)
					HVCRA_local[i].raiseExcWeight(update_Ge_RA_global[Id_RA_local[i]]); // update conductance

				// update conductance arrays and fired indicators
				some_RA_exc_conductance_was_updated_global = 0;
				some_RA_exc_conductance_was_updated_local = 0;
				
				std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
			}	

			////////////////////////////////////////////
            //// Process HVC-RA dendritic spikes
            ////////////////////////////////////////////
			if ( some_RA_neuron_bursted_global > 0 )
            {
				// update conductance arrays and fired indicators
            	some_RA_neuron_bursted_local = 0;
	        	some_RA_neuron_bursted_global = 0;
	        	
	        	// gather all bursted HVC-RA neurons
				this->gather_spiked_or_bursted_neurons(RA_neurons_bursted_local, RA_neurons_bursted_global);
				
				// add dendritic spikes for bursted neurons
				for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
					previous_dendritic_spike_times_global[RA_neurons_bursted_global[j]].push_back(internal_time);
				
				////////////////////////////////////////////////////////////////////////////////////
				// LTP or LTD of previously delivered somatic spike times on bursted HVC-RA neuron
				////////////////////////////////////////////////////////////////////////////////////
				
				for (int i = 0; i < N_RA_local; i++)
                {
					int presyn_ID = Id_RA_local[i]; // real id of presynaptic neuron
					
                    // if neuron is saturated apply LTP only if spiked neurons are among supersynapse targets
                    if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
                    {
                        for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_bursted_global[j]; // id of postsynaptic neuron
                            
                            // do not allow self-to-self synapses to emerge
							if ( presyn_ID != postsyn_ID )
							{
								std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
											supersynapses_local[i].end(), postsyn_ID );

								if ( pos!=supersynapses_local[i].end() )
								{
									// find previous somatic spikes that are not too old
									if ( !delivered_spike_times[i][postsyn_ID].empty() )
									{
										std::vector<double>::iterator it_relevant_spikes = std::lower_bound(delivered_spike_times[i][postsyn_ID].begin(),
																											delivered_spike_times[i][postsyn_ID].end(),
																											internal_time - STDP_WINDOW);
																							
										for (auto it = it_relevant_spikes; it != delivered_spike_times[i][postsyn_ID].end(); it++)
										{
											double dt = internal_time - *it;
									
										
											//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
											
											if (dt <= synaptic_params.T_0)
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												LTD(weights_RA_RA_local[i][postsyn_ID], dt);
												
												//~ std::cout   << "LTD from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
											}
											else
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												
												LTP(weights_RA_RA_local[i][postsyn_ID], dt);
												
												//~ std::cout   << "LTP from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
												
											}	
											//double w = weights_local[i][fired_ID];
											
											//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
											 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
											 //           dt, weights_local[i][fired_ID] - w);
											
											update_synapse(i, postsyn_ID);
											
										}
									}
								}
							}
                        }

                    }
                    // if not saturated apply LTP for all dendritic spikes
                    else
                    {
                        for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_bursted_global[j];
                            // don't allow self-to-self connections
                            if ( postsyn_ID != presyn_ID )
                            {
                                // find previous somatic spikes that are not too old
								if ( !delivered_spike_times[i][postsyn_ID].empty() )
								{
									std::vector<double>::iterator it_relevant_spikes = std::lower_bound(delivered_spike_times[i][postsyn_ID].begin(),
																										delivered_spike_times[i][postsyn_ID].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it = it_relevant_spikes; it != delivered_spike_times[i][postsyn_ID].end(); it++)
									{
										double dt = internal_time - *it;
								

										//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
										
									    if (dt <= synaptic_params.T_0)
									    {
											double w_before = weights_RA_RA_local[i][postsyn_ID];
												
											LTD(weights_RA_RA_local[i][postsyn_ID], dt);
											
											//~ // show LTD on pool neurons only
											//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											//~ 
											//~ if (p.first == p.second)
											//~ {
												//~ std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											//~ }
											//~ 
											//~ std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
											
										}
										else
										{
											double w_before = weights_RA_RA_local[i][postsyn_ID];
											
											LTP(weights_RA_RA_local[i][postsyn_ID], dt);
											
											//~ // show LTP on pool neurons only
											//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											//~ 
											//~ if (p.first == p.second)
											//~ {
												//~ std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											//~ }
											
											
											//~ std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
										}	
										//double w = weights_local[i][fired_ID];
										
										//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
										 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
										 //           dt, weights_local[i][fired_ID] - w);
										
										update_synapse(i, postsyn_ID);
									
									}
									
									
								}
                                
                            }
                        }
                   }

                } // end for i -> N_RA_local

				       

               // check if we need axon remodeling

                for (int i = 0; i < N_RA_local; i++)
                {
                    if ( (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss) && (remodeled_local[i] == 0) )
                    {
					    this->axon_remodeling(i);

				    }
                }
            
				// clear bursts
				RA_neurons_bursted_local.clear();
				RA_neurons_bursted_global.clear();
				
            } // end if some HVC-RA bursted
            
            
            network_time += NETWORK_UPDATE_FREQUENCY;
        }
    }
    
    this->potentiation_decay_sudden_maturation();
    //printf("After potentiation decay")
    this->update_all_synapses_sudden_maturation();
	
	// update maturation info
	std::vector<int> RA_matured_local;
	std::vector<int> RA_matured_global;
	
	for (int i = 0; i < N_RA_local; i++)
	{
		num_spikes_in_recent_trials_local[i].push_front(static_cast<int>(spikes_in_trial_soma_local[i].size()));
        
        //firing_rate_short_local[i] = std::accumulate(num_spikes_in_recent_trials[i].begin(), num_spikes_in_recent_trials[i].begin() + RATE_WINDOW_SHORT, 0.0)
          //                          / static_cast<double>(RATE_WINDOW_SHORT);

		// calculate firing rate in large window:
		int maturation_window = 300;
		double maturation_threshold = 0.95;
		
        firing_rate_long_local[i] = std::accumulate(num_spikes_in_recent_trials_local[i].begin(), num_spikes_in_recent_trials_local[i].end(), 0.0)
                                                                                            / static_cast<double>(maturation_params.RATE_WINDOW_LONG);
        if ( mature_global[Id_RA_local[i]] != 1)
        {
			
			std::vector<double> recent_firing_robustness(maturation_window);
			
			for (int j = 0; j < maturation_window; j++ )
				if ( num_spikes_in_recent_trials_local[i][j] > 0 )
					recent_firing_robustness[j] = 1;
			
			double firing_robustness = std::accumulate(recent_firing_robustness.begin(), recent_firing_robustness.end(), 0.0) / static_cast<double>(maturation_window);
			
			if ( firing_robustness > maturation_threshold )
			{
				RA_matured_local.push_back(Id_RA_local[i]);
					
				this->set_neuron_mature(i);
			}                                                                                                                          
			//if (std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].end(), 0.0) >= 1)
			//{
			//	std::cout << "Recent bursts of neuron " << Id_RA_local[i] << " :" << std::endl;
				
			//	for (size_t j = 0; j < num_bursts_in_recent_trials[i].size(); j++)
			//		std::cout << num_bursts_in_recent_trials[i][j] << "\t";
				
			//	std::cout << std::endl;
			//}
		}
	}
	
	this->gather_spiked_or_bursted_neurons(RA_matured_local, RA_matured_global);
	
	for (size_t i = 0; i < RA_matured_global.size(); i++)
	{
		int neuron_id = RA_matured_global[i];
		
		mature_global[neuron_id] = 1;
		this->rescale_synapses_to_mature(neuron_id);
	}
}

void NetworkGrowthSimulator::trial_soma_pre_dend_post_stdp_delays_sudden_maturation(bool training)
{
	// indicators for conductance updates and bursting
	int some_RA_inh_conductance_was_updated_global = 0;
	int some_RA_inh_conductance_was_updated_local = 0;
	
	int some_RA_exc_conductance_was_updated_global = 0;
	int some_RA_exc_conductance_was_updated_local = 0;
	
	int some_I_exc_conductance_was_updated_global = 0;
	int some_I_exc_conductance_was_updated_local = 0;
	
	int some_RA_neuron_bursted_global = 0;
	int some_RA_neuron_bursted_local = 0;
	
	 // arrays with updates for conductances
    std::vector<double> update_Ge_RA_local(N_RA);
	std::vector<double> update_Ge_RA_global(N_RA);

	std::vector<double> update_Gi_RA_local(N_RA);
	std::vector<double> update_Gi_RA_global(N_RA);

	std::vector<double> update_Ge_I_local(N_I);
	std::vector<double> update_Ge_I_global(N_I);
	
	std::vector<int> RA_neurons_bursted_local; // local array of HVC-RA neurons that bursted
	std::vector<int> RA_neurons_bursted_global; // global array of HVC-RA neurons that bursted
	
	// sample training innervation time and send to all processes
    double training_kick_time;
    
	if (MPI_rank == 0)
	{
		training_kick_time = WAITING_TIME + noise_generator.random(TRIAL_DURATION - 2*WAITING_TIME);
		std::cout << "training_kick_time = " << training_kick_time << std::endl;
    }
    MPI_Bcast(&training_kick_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    bool training_excited = false; // indicator that training neurons were already excited
    
    double internal_time = 0;
    double network_time = NETWORK_UPDATE_FREQUENCY;

   

    //printf("network time = %f\n", network_time);

    
	// initialize update arrays and fired indicators
	std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
	
	
    // evolve dynamics
    int num_steps = static_cast<int>(TRIAL_DURATION / TIMESTEP);
    
    for (int t = 1; t < num_steps; t++)
	{
		internal_time += TIMESTEP;
	
		if ( ( training ) && ( !training_excited ) && ( internal_time >= training_kick_time ) )
		{
			for (int i = 0; i < N_TR; i++)
			{
				int rank;
				int shift;
				
				this->get_neuronRA_location(training_neurons[i], &rank, &shift);
				
				if (MPI_rank == rank)
					HVCRA_local[shift].raiseE(G_TRAINING_KICK);
			}
			
			training_excited = true;
		}
		
		//////////////////////
		// HVC-RA neurons
		//////////////////////
		
		for (int i = 0; i < N_RA_local; i++)
		{
            // Debraband step
            HVCRA_local[i].Debraband_step_no_target_update();
            
           
			///////////////////////////////
			// HVC-RA somatic spike
			///////////////////////////////
			// if some neuron produced somatic spike
            if (HVCRA_local[i].get_fired_soma())
            {
                spikes_in_trial_soma_local[i].push_back(internal_time);
                
				// update delivery queues
				
				// for inhibitory neurons
				// loop over all inhibitory targets of fired neurons
				size_t num_I_targets = syn_ID_RA_I_local[i].size();
				
				for (size_t j = 0; j < num_I_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_RA_I_local[i][j];
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_RA_I[i].empty() )
						delivery_queue_RA_I[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_RA_I[i].begin(), delivery_queue_RA_I[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_RA_I[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
					//std::cout << "neuron " << Id_RA_local[i] << " spike time = " << internal_time << " axonal delay = " << axonal_delays_RA_I[Id_RA_local[i]][syn_ID] << " delivery time RA to I: " << delivery_queue_RA_I[i].back().first << " delivery target id: " << syn_ID_RA_I_local[i][delivery_queue_RA_I[i].back().second] << std::endl;
				
				
				// if neuron is saturated, deliver spikes only to supersynaptic targets
				if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
				{
					for (size_t j = 0; j < supersynapses_local[i].size(); j++)
					{
						int target_id = supersynapses_local[i][j];
						
						double delivery_time = internal_time + axonal_delays_RA_RA_local[i][target_id];
						
						// if queue is empty, just add item to the queue
						if ( delivery_queue_RA_RA_soma[i].empty() )
							delivery_queue_RA_RA_soma[i].push_back(std::pair<double,int>(delivery_time, target_id));	
						// otherwise add item so that queue is sorted
						else
						{
							auto it = std::upper_bound(delivery_queue_RA_RA_soma[i].begin(), delivery_queue_RA_RA_soma[i].end(), std::pair<double,int>(delivery_time, target_id));
							delivery_queue_RA_RA_soma[i].insert(it, std::pair<double,int>(delivery_time, target_id));
						}
					}
				}
				else // deliver spikes to everyone except itself
				{	
					for (int j = 0; j < N_RA; j++)
					{
						if ( j != Id_RA_local[i])
						{
						
							double delivery_time = internal_time + axonal_delays_RA_RA_local[i][j];
						
							// if queue is empty, just add item to the queue
							if ( delivery_queue_RA_RA_soma[i].empty() )
								delivery_queue_RA_RA_soma[i].push_back(std::pair<double,int>(delivery_time, j));	
							// otherwise add item so that queue is sorted
							else
							{
								auto it = std::upper_bound(delivery_queue_RA_RA_soma[i].begin(), delivery_queue_RA_RA_soma[i].end(), std::pair<double,int>(delivery_time, j));
								delivery_queue_RA_RA_soma[i].insert(it, std::pair<double,int>(delivery_time, j));
							}
						}
					}
				}
                
            } // end if local HVC-RA neuron spiked  
            
             //////////////////////////////////
			// HVC-RA -> HVC-I delivery queue
			//////////////////////////////////
            // check if somatic spike was delivered to some interneuron
            // loop through the delivery queue to check if current time exceeds the spike delivery time
            std::vector<std::pair<double,int>>::iterator it = delivery_queue_RA_I[i].begin();
			
            for (; it != delivery_queue_RA_I[i].end(); it++)
            {
				if (internal_time >= it->first)
				{
					
					some_I_exc_conductance_was_updated_local = 1;
					
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_RA_I_local[i][pos_in_local_target_array];
					
					update_Ge_I_local[target_id] += weights_RA_I_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(RA) neuron " << Id_RA_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to HVC(I) neuron " << syn_ID_RA_I_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_I[i].erase(delivery_queue_RA_I[i].begin(), it);
			
			//////////////////////////////////////////////////
			// HVC-RA -> HVC-RA somatic spike delivery queue
			//////////////////////////////////////////////////
			// check if somatic spike was delivered to some HVC(RA) neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
            it = delivery_queue_RA_RA_soma[i].begin();
            
            for (; it != delivery_queue_RA_RA_soma[i].end(); it++)
            {
				double delivery_time = it->first;
				
				if ( internal_time >= delivery_time )
				{
					int target_id = it->second;
					
					delivered_spike_times[i][target_id].push_back(delivery_time);
					
					// update conductance of target if synapse is active
					if ( active_indicators_local[i][target_id] == 1 )
					{
						some_RA_exc_conductance_was_updated_local = 1;
						update_Ge_RA_local[target_id] += weights_RA_RA_local[i][target_id];
					}
					
					//////////////
					//// LTD
					//////////////
					// do LTD on dendritic spike time of target neuron if it is not mature neuron
					if ( mature_global[target_id] != 1 )
					{
						
						 // if neuron is saturated apply LTD only if target is among super synapses
						if (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss)
						{
							if ( supersynapses_indicators_local[i][target_id] == 1 )
							{	
								// find previous dendritic spikes that are not too old
								if ( !previous_dendritic_spike_times_global[target_id].empty() )
								{
									std::vector<double>::iterator it_relevant_dend_spikes = std::lower_bound(previous_dendritic_spike_times_global[target_id].begin(),
																										previous_dendritic_spike_times_global[target_id].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it_dend_spike = it_relevant_dend_spikes; it_dend_spike != previous_dendritic_spike_times_global[target_id].end(); it_dend_spike++)
									{
										double dt = *it_dend_spike - delivery_time;

										//std::cout << "dt in saturated LTD = " << dt << std::endl;

										
										double w_before = weights_RA_RA_local[i][target_id];     
							   
										LTD(weights_RA_RA_local[i][target_id], dt);
										
										
										//~ std::cout << "LTD from saturated " << Id_RA_local[i] << " -> " << supersynapse_id
												  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][supersynapse_id] - w_before
												  //~ << std::endl;
												  
										//printf("LTD from saturated neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n",
										//            Id_RA_local[i], supersynapse_id, spikes_in_trial_soma_local[i].back(), spike_times_dend_global[supersynapse_id],
										//            dt, weights_local[i][supersynapse_id] - w);
										
										update_synapse(i, target_id);	
										
									}   
								}
							}

							// if some supersynapse desaturated, update all synapses
							if (static_cast<int>(supersynapses_local[i].size()) < synaptic_params.Nss)
								for (int j = 0; j < N_RA; j++)
									this->update_synapse(i, j);
						}
						// if not saturated apply LTD rule 
						else
						{
							// find previous dendritic spikes that are not too old
							if ( !previous_dendritic_spike_times_global[target_id].empty() )
							{
								std::vector<double>::iterator it_relevant_dend_spikes = std::lower_bound(previous_dendritic_spike_times_global[target_id].begin(),
																									previous_dendritic_spike_times_global[target_id].end(),
																									internal_time - STDP_WINDOW);
																					
								for (auto it_dend_spike = it_relevant_dend_spikes; it_dend_spike != previous_dendritic_spike_times_global[target_id].end(); it_dend_spike++)
								{
									double dt = *it_dend_spike - delivery_time;

									//std::cout << "dt in saturated LTD = " << dt << std::endl;

									
									double w_before = weights_RA_RA_local[i][target_id];     
						   
									LTD(weights_RA_RA_local[i][target_id], dt);
									
									//~ // show LTD of pool neurons only
									//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
									//~ 
									//~ if (p.first == p.second)
									//~ {
										//~ std::cout << "LTD from  " << Id_RA_local[i] << " -> " << target_id
											  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][target_id] - w_before
											  //~ << " w_after = " << weights_RA_RA_local[i][target_id] << std::endl;
									//~ }
								//~ 
									update_synapse(i, target_id);	
									
								}   
							}
						}
					} // end if target neuron is not mature
				} // end if internal_time >= delivery spike time 
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_RA_soma[i].erase(delivery_queue_RA_RA_soma[i].begin(), it);
			
			
			///////////////////////////////
			// HVC-RA dendritic spike
			///////////////////////////////
			//if some neuron produced dendritic spike, store this neuron in array
			
			if ( HVCRA_local[i].get_fired_dend() )
			{
				// show spikes of pool neurons only
				auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
				
				if (p.first == p.second)
					std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
				
				
				spikes_in_trial_dend_local[i].push_back(internal_time); // dend spike time relative to trial onset
				
				if ( mature_global[Id_RA_local[i]] != 1 )
				{
					RA_neurons_bursted_local.push_back(Id_RA_local[i]);
				
					some_RA_neuron_bursted_local = 1;
				}
				//std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
			} // end if get_fired_dend
	
            
		} // end for i = 0 -> N_RA_local (loop through all local HVC-RA neurons)
		
		//////////////////////
		// HVC-I neurons
		//////////////////////
		for (int i = 0; i < N_I_local; i++)
		{
			HVCI_local[i].DP8_step_no_target_update();
			
			
			///////////////////////////////
			// HVC-I spike
			///////////////////////////////
			//  if some I neuron spikes, update delivery queue
			if (HVCI_local[i].get_fired())
			{
				//printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
				spikes_in_trial_interneuron_local[i].push_back(internal_time);

				size_t num_RA_targets = syn_ID_I_RA_local[i].size();
				// loop over all targets of fired neurons
				for (size_t j = 0; j < num_RA_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_I_RA_local[i][j];
					
					//std::cout << "HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << syn_ID_I_RA_local[i][j] << " spike time " << internal_time << " delivery time = " << delivery_time << std::endl; 
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_I_RA[i].empty() )
						delivery_queue_I_RA[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_I_RA[i].begin(), delivery_queue_I_RA[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_I_RA[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
			}
			
			//////////////////////////////////
			// HVC-I -> HVC-RA delivery queue
			//////////////////////////////////
			// check if interneuron spike was delivered to some HVC-RA neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
			std::vector<std::pair<double,int>>::iterator it = delivery_queue_I_RA[i].begin();
			
			for (; it != delivery_queue_I_RA[i].end(); it++)
			{
				if (internal_time >= it->first)
				{
					some_RA_inh_conductance_was_updated_local = 1;
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_I_RA_local[i][pos_in_local_target_array];
					
					//std::cout << "Delivered spike HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << target_id << " delivered_time " << internal_time << " delivery time in queue = " << it->first << std::endl; 
					
					
					update_Gi_RA_local[target_id] += weights_I_RA_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(I) neuron " << Id_I_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to neuron " << syn_ID_I_RA_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_I_RA[i].erase(delivery_queue_I_RA[i].begin(), it);
		
		}  // end for i = 0 -> N_I_local (loop through all local HVC-I neurons)
		
		/////////////////////////////
		// Network Synchronization //
		/////////////////////////////
		if ( ( internal_time > network_time ) || ( t == num_steps-1 ) )
        {
            //if (MPI_rank == 0)
            //{
            //    std::cout << "internal time = " << internal_time << " > network_time = " << network_time << std::endl; 
            //}

            MPI_Allreduce(&some_RA_exc_conductance_was_updated_local, &some_RA_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_RA_inh_conductance_was_updated_local, &some_RA_inh_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_I_exc_conductance_was_updated_local, &some_I_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_RA_neuron_bursted_local, &some_RA_neuron_bursted_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            
        
			
            ////////////////////////////////////////////
            //// Process HVC-I -> HVC-RA interactions
            ////////////////////////////////////////////
            if ( some_RA_inh_conductance_was_updated_global > 0 )
            {
            // sum update array and send to all processes

                MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                for (int i = 0; i < N_RA_local; i++)
                {
					//std::cout << "HVC-RA neuron " << Id_RA_local[i] << " raised inhibitory conductance by " << update_Gi_RA_global[Id_RA_local[i]] << std::endl;
                    HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
            	}
				// update conductance arrays and delivered indicators
				some_RA_inh_conductance_was_updated_global = 0;
				some_RA_inh_conductance_was_updated_local = 0;
				
				std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            }
            
            ////////////////////////////////////////////
            //// Process HVC-RA -> HVC-I interactions
            ////////////////////////////////////////////
			if ( some_I_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons

				for (int i = 0; i < N_I_local; i++)
					HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);

				// update conductance arrays and fired indicators
				some_I_exc_conductance_was_updated_local = 0;
				some_I_exc_conductance_was_updated_global = 0;

				std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
			}
			
			////////////////////////////////////////////
            //// Process HVC-RA -> HVC-RA interactions
            ////////////////////////////////////////////
			if ( some_RA_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons
				for (int i = 0; i < N_RA_local; i++)
					HVCRA_local[i].raiseExcWeight(update_Ge_RA_global[Id_RA_local[i]]); // update conductance

				// update conductance arrays and fired indicators
				some_RA_exc_conductance_was_updated_global = 0;
				some_RA_exc_conductance_was_updated_local = 0;
				
				std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
			}	

			////////////////////////////////////////////
            //// Process HVC-RA dendritic spikes
            ////////////////////////////////////////////
			if ( some_RA_neuron_bursted_global > 0 )
            {
				// update conductance arrays and fired indicators
            	some_RA_neuron_bursted_local = 0;
	        	some_RA_neuron_bursted_global = 0;
	        	
	        	// gather all bursted HVC-RA neurons
				this->gather_spiked_or_bursted_neurons(RA_neurons_bursted_local, RA_neurons_bursted_global);
				
				// add dendritic spikes for bursted neurons
				for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
					previous_dendritic_spike_times_global[RA_neurons_bursted_global[j]].push_back(internal_time);
				
				////////////////////////////////////////////////////////////////////////////////////
				// LTP or LTD of previously delivered somatic spike times on bursted HVC-RA neuron
				////////////////////////////////////////////////////////////////////////////////////
				
				for (int i = 0; i < N_RA_local; i++)
                {
					int presyn_ID = Id_RA_local[i]; // real id of presynaptic neuron
					
                    // if neuron is saturated apply LTP only if spiked neurons are among supersynapse targets
                    if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
                    {
                        for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_bursted_global[j]; // id of postsynaptic neuron
                            
                            // do not allow self-to-self synapses to emerge
							if ( presyn_ID != postsyn_ID )
							{
								std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
											supersynapses_local[i].end(), postsyn_ID );

								if ( pos!=supersynapses_local[i].end() )
								{
									// find previous somatic spikes that are not too old
									if ( !delivered_spike_times[i][postsyn_ID].empty() )
									{
										std::vector<double>::iterator it_relevant_spikes = std::lower_bound(delivered_spike_times[i][postsyn_ID].begin(),
																											delivered_spike_times[i][postsyn_ID].end(),
																											internal_time - STDP_WINDOW);
																							
										for (auto it = it_relevant_spikes; it != delivered_spike_times[i][postsyn_ID].end(); it++)
										{
											double dt = internal_time - *it;
									
										
											//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
											
											if (dt <= synaptic_params.T_0)
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												LTD(weights_RA_RA_local[i][postsyn_ID], dt);
												
												//~ std::cout   << "LTD from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
											}
											else
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												
												LTP(weights_RA_RA_local[i][postsyn_ID], dt);
												
												//~ std::cout   << "LTP from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
												
											}	
											//double w = weights_local[i][fired_ID];
											
											//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
											 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
											 //           dt, weights_local[i][fired_ID] - w);
											
											update_synapse(i, postsyn_ID);
											
										}
									}
								}
							}
                        }

                    }
                    // if not saturated apply LTP for all dendritic spikes
                    else
                    {
                        for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_bursted_global[j];
                            // don't allow self-to-self connections
                            if ( postsyn_ID != presyn_ID )
                            {
                                // find previous somatic spikes that are not too old
								if ( !delivered_spike_times[i][postsyn_ID].empty() )
								{
									std::vector<double>::iterator it_relevant_spikes = std::lower_bound(delivered_spike_times[i][postsyn_ID].begin(),
																										delivered_spike_times[i][postsyn_ID].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it = it_relevant_spikes; it != delivered_spike_times[i][postsyn_ID].end(); it++)
									{
										double dt = internal_time - *it;
								

										//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
										
									    if (dt <= synaptic_params.T_0)
									    {
											double w_before = weights_RA_RA_local[i][postsyn_ID];
												
											LTD(weights_RA_RA_local[i][postsyn_ID], dt);
											
											//~ // show LTD on pool neurons only
											//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											//~ 
											//~ if (p.first == p.second)
											//~ {
												//~ std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											//~ }
											
											//~ std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
											
										}
										else
										{
											double w_before = weights_RA_RA_local[i][postsyn_ID];
											
											LTP(weights_RA_RA_local[i][postsyn_ID], dt);
											
											//~ // show LTP on pool neurons only
											//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											//~ 
											//~ if (p.first == p.second)
											//~ {
												//~ std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											//~ }
											
											
											//~ std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
										}	
										//double w = weights_local[i][fired_ID];
										
										//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
										 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
										 //           dt, weights_local[i][fired_ID] - w);
										
										update_synapse(i, postsyn_ID);
									
									}
									
									
								}
                                
                            }
                        }
                   }

                } // end for i -> N_RA_local

				       

               // check if we need axon remodeling

                for (int i = 0; i < N_RA_local; i++)
                {
                    if ( (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss) && (remodeled_local[i] == 0) )
                    {
					    this->axon_remodeling(i);

				    }
                }
            
				// clear bursts
				RA_neurons_bursted_local.clear();
				RA_neurons_bursted_global.clear();
				
            } // end if some HVC-RA bursted
            
            
            network_time += NETWORK_UPDATE_FREQUENCY;
        }
    }
    
    this->potentiation_decay_sudden_maturation();
    //printf("After potentiation decay")
    this->update_all_synapses_sudden_maturation();
	
	// update maturation info
	std::vector<int> RA_matured_local;
	std::vector<int> RA_matured_global;
	
	for (int i = 0; i < N_RA_local; i++)
	{
		num_spikes_in_recent_trials_local[i].push_front(static_cast<int>(spikes_in_trial_soma_local[i].size()));
        
        //firing_rate_short_local[i] = std::accumulate(num_spikes_in_recent_trials[i].begin(), num_spikes_in_recent_trials[i].begin() + RATE_WINDOW_SHORT, 0.0)
          //                          / static_cast<double>(RATE_WINDOW_SHORT);

		// calculate firing rate in large window:
		int maturation_window = 300;
		double maturation_threshold = 0.95;
		
        firing_rate_long_local[i] = std::accumulate(num_spikes_in_recent_trials_local[i].begin(), num_spikes_in_recent_trials_local[i].end(), 0.0)
                                                                                            / static_cast<double>(maturation_params.RATE_WINDOW_LONG);
        if ( mature_global[Id_RA_local[i]] != 1)
        {
			
			std::vector<double> recent_firing_robustness(maturation_window);
			
			for (int j = 0; j < maturation_window; j++ )
				if ( num_spikes_in_recent_trials_local[i][j] > 0 )
					recent_firing_robustness[j] = 1;
			
			double firing_robustness = std::accumulate(recent_firing_robustness.begin(), recent_firing_robustness.end(), 0.0) / static_cast<double>(maturation_window);
			
			if ( firing_robustness > maturation_threshold )
			{
				RA_matured_local.push_back(Id_RA_local[i]);
					
				this->set_neuron_mature(i);
			}                                                                                                                          
			//if (std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].end(), 0.0) >= 1)
			//{
			//	std::cout << "Recent bursts of neuron " << Id_RA_local[i] << " :" << std::endl;
				
			//	for (size_t j = 0; j < num_bursts_in_recent_trials[i].size(); j++)
			//		std::cout << num_bursts_in_recent_trials[i][j] << "\t";
				
			//	std::cout << std::endl;
			//}
		}
	}
	
	this->gather_spiked_or_bursted_neurons(RA_matured_local, RA_matured_global);
	
	for (size_t i = 0; i < RA_matured_global.size(); i++)
	{
		int neuron_id = RA_matured_global[i];
		
		mature_global[neuron_id] = 1;
		this->rescale_synapses_to_mature(neuron_id);
	}
}

void NetworkGrowthSimulator::trial_no_stdp_fixedSpread(std::vector<double>& spread_times)
{
	// indicators for conductance updates and bursting
	int some_RA_inh_conductance_was_updated_global = 0;
	int some_RA_inh_conductance_was_updated_local = 0;
	
	int some_RA_exc_conductance_was_updated_global = 0;
	int some_RA_exc_conductance_was_updated_local = 0;
	
	int some_I_exc_conductance_was_updated_global = 0;
	int some_I_exc_conductance_was_updated_local = 0;
	
	 // arrays with updates for conductances
    std::vector<double> update_Ge_RA_local(N_RA);
	std::vector<double> update_Ge_RA_global(N_RA);

	std::vector<double> update_Gi_RA_local(N_RA);
	std::vector<double> update_Gi_RA_global(N_RA);

	std::vector<double> update_Ge_I_local(N_I);
	std::vector<double> update_Ge_I_global(N_I);
	
	std::vector<int> RA_neurons_bursted_local; // local array of HVC-RA neurons that bursted
	std::vector<int> RA_neurons_bursted_global; // global array of HVC-RA neurons that bursted
	
	
	// sample training innervation time and send to all processes
  // training neurons innervation
	std::vector<bool> indicators_training(N_RA);
	std::fill(indicators_training.begin(), indicators_training.end(), false);
	
	for (int i = 0; i < N_TR; i++)
		indicators_training[training_neurons[i]] = true;
		
	std::vector<bool> indicators_current_injected(N_RA);
	std::fill(indicators_current_injected.begin(), indicators_current_injected.end(), false);
	
	// sample mean training innervation time and send to all processes
    double mean_injection_time = 100;
    
    
    double internal_time = 0;
    double network_time = NETWORK_UPDATE_FREQUENCY;

   

    //printf("network time = %f\n", network_time);

    
	// initialize update arrays and fired indicators
	std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
	
	
    // evolve dynamics
    int num_steps = static_cast<int>(TRIAL_DURATION / TIMESTEP);
    
    for (int t = 1; t < num_steps; t++)
	{
		internal_time += TIMESTEP;
		
		//////////////////////
		// HVC-RA neurons
		//////////////////////
		
		for (int i = 0; i < N_RA_local; i++)
		{
			
			if ( ( indicators_training[Id_RA_local[i]] ) && ( !indicators_current_injected[Id_RA_local[i]] ) )
			{
				if ( internal_time > spread_times[Id_RA_local[i]] + mean_injection_time )
				{
					HVCRA_local[i].raiseE(G_TRAINING_KICK);
					indicators_current_injected[Id_RA_local[i]] = true;
				}
			}
            // Debraband step
            HVCRA_local[i].Debraband_step_no_target_update();
            
           
			///////////////////////////////
			// HVC-RA somatic spike
			///////////////////////////////
			// if some neuron produced somatic spike
            if (HVCRA_local[i].get_fired_soma())
            {
                spikes_in_trial_soma_local[i].push_back(internal_time);
                
				// update delivery queues
				
				// for inhibitory neurons
				// loop over all inhibitory targets of fired neurons
				size_t num_I_targets = syn_ID_RA_I_local[i].size();
				
				for (size_t j = 0; j < num_I_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_RA_I_local[i][j];
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_RA_I[i].empty() )
						delivery_queue_RA_I[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_RA_I[i].begin(), delivery_queue_RA_I[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_RA_I[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
					//std::cout << "neuron " << Id_RA_local[i] << " spike time = " << internal_time << " axonal delay = " << axonal_delays_RA_I[Id_RA_local[i]][syn_ID] << " delivery time RA to I: " << delivery_queue_RA_I[i].back().first << " delivery target id: " << syn_ID_RA_I_local[i][delivery_queue_RA_I[i].back().second] << std::endl;
				
				// deliver spikes to active targets
				size_t num_RA_targets = active_synapses_local[i].size();
				
				for (size_t j = 0; j < num_RA_targets; j++)
				{	
					int target_id = active_synapses_local[i][j];
					double delivery_time = internal_time + axonal_delays_RA_RA_local[i][target_id];
				
					// if queue is empty, just add item to the queue
					if ( delivery_queue_RA_RA_soma[i].empty() )
						delivery_queue_RA_RA_soma[i].push_back(std::pair<double,int>(delivery_time, target_id));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_RA_RA_soma[i].begin(), delivery_queue_RA_RA_soma[i].end(), std::pair<double,int>(delivery_time, target_id));
						delivery_queue_RA_RA_soma[i].insert(it, std::pair<double,int>(delivery_time, target_id));
					}
				}
                
            } // end if local HVC-RA neuron spiked  
            
             //////////////////////////////////
			// HVC-RA -> HVC-I delivery queue
			//////////////////////////////////
            // check if somatic spike was delivered to some interneuron
            // loop through the delivery queue to check if current time exceeds the spike delivery time
            std::vector<std::pair<double,int>>::iterator it = delivery_queue_RA_I[i].begin();
			
            for (; it != delivery_queue_RA_I[i].end(); it++)
            {
				if (internal_time >= it->first)
				{
					
					some_I_exc_conductance_was_updated_local = 1;
					
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_RA_I_local[i][pos_in_local_target_array];
					
					update_Ge_I_local[target_id] += weights_RA_I_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(RA) neuron " << Id_RA_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to HVC(I) neuron " << syn_ID_RA_I_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_I[i].erase(delivery_queue_RA_I[i].begin(), it);
			
			//////////////////////////////////////////////////
			// HVC-RA -> HVC-RA somatic spike delivery queue
			//////////////////////////////////////////////////
			// check if somatic spike was delivered to some HVC(RA) neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
            it = delivery_queue_RA_RA_soma[i].begin();
            
            for (; it != delivery_queue_RA_RA_soma[i].end(); it++)
            {
				double delivery_time = it->first;
				
				if ( internal_time >= delivery_time )
				{
					int target_id = it->second;
					
					
					some_RA_exc_conductance_was_updated_local = 1;
					update_Ge_RA_local[target_id] += weights_RA_RA_local[i][target_id];
					
					
				} // end if internal_time >= delivery spike time 
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_RA_soma[i].erase(delivery_queue_RA_RA_soma[i].begin(), it);
			
			
			///////////////////////////////
			// HVC-RA dendritic spike
			///////////////////////////////
			//if some neuron produced dendritic spike, store this neuron in array
			
			if ( HVCRA_local[i].get_fired_dend() )
			{
				//~ // show spikes of pool neurons only
				//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
				//~ 
				//~ if (p.first == p.second)
					//~ std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
				//~ 
				//~ 
				spikes_in_trial_dend_local[i].push_back(internal_time); // dend spike time relative to trial onset
				
				//std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
			} // end if get_fired_dend
	
            
		} // end for i = 0 -> N_RA_local (loop through all local HVC-RA neurons)
		
		//////////////////////
		// HVC-I neurons
		//////////////////////
		for (int i = 0; i < N_I_local; i++)
		{
			HVCI_local[i].DP8_step_no_target_update();
			
			
			///////////////////////////////
			// HVC-I spike
			///////////////////////////////
			//  if some I neuron spikes, update delivery queue
			if (HVCI_local[i].get_fired())
			{
				//printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
				spikes_in_trial_interneuron_local[i].push_back(internal_time);

				size_t num_RA_targets = syn_ID_I_RA_local[i].size();
				// loop over all targets of fired neurons
				for (size_t j = 0; j < num_RA_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_I_RA_local[i][j];
					
					//std::cout << "HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << syn_ID_I_RA_local[i][j] << " spike time " << internal_time << " delivery time = " << delivery_time << std::endl; 
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_I_RA[i].empty() )
						delivery_queue_I_RA[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_I_RA[i].begin(), delivery_queue_I_RA[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_I_RA[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
			}
			
			//////////////////////////////////
			// HVC-I -> HVC-RA delivery queue
			//////////////////////////////////
			// check if interneuron spike was delivered to some HVC-RA neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
			std::vector<std::pair<double,int>>::iterator it = delivery_queue_I_RA[i].begin();
			
			for (; it != delivery_queue_I_RA[i].end(); it++)
			{
				if (internal_time >= it->first)
				{
					some_RA_inh_conductance_was_updated_local = 1;
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_I_RA_local[i][pos_in_local_target_array];
					
					//std::cout << "Delivered spike HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << target_id << " delivered_time " << internal_time << " delivery time in queue = " << it->first << std::endl; 
					
					
					update_Gi_RA_local[target_id] += weights_I_RA_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(I) neuron " << Id_I_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to neuron " << syn_ID_I_RA_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_I_RA[i].erase(delivery_queue_I_RA[i].begin(), it);
		
		}  // end for i = 0 -> N_I_local (loop through all local HVC-I neurons)
		
		/////////////////////////////
		// Network Synchronization //
		/////////////////////////////
		if ( ( internal_time > network_time ) || ( t == num_steps-1 ) )
        {
            //if (MPI_rank == 0)
            //{
            //    std::cout << "internal time = " << internal_time << " > network_time = " << network_time << std::endl; 
            //}

            MPI_Allreduce(&some_RA_exc_conductance_was_updated_local, &some_RA_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_RA_inh_conductance_was_updated_local, &some_RA_inh_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_I_exc_conductance_was_updated_local, &some_I_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
			
            ////////////////////////////////////////////
            //// Process HVC-I -> HVC-RA interactions
            ////////////////////////////////////////////
            if ( some_RA_inh_conductance_was_updated_global > 0 )
            {
            // sum update array and send to all processes

                MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                for (int i = 0; i < N_RA_local; i++)
                {
					//std::cout << "HVC-RA neuron " << Id_RA_local[i] << " raised inhibitory conductance by " << update_Gi_RA_global[Id_RA_local[i]] << std::endl;
                    HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
            	}
				// update conductance arrays and delivered indicators
				some_RA_inh_conductance_was_updated_global = 0;
				some_RA_inh_conductance_was_updated_local = 0;
				
				std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            }
            
            ////////////////////////////////////////////
            //// Process HVC-RA -> HVC-I interactions
            ////////////////////////////////////////////
			if ( some_I_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons

				for (int i = 0; i < N_I_local; i++)
					HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);

				// update conductance arrays and fired indicators
				some_I_exc_conductance_was_updated_local = 0;
				some_I_exc_conductance_was_updated_global = 0;

				std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
			}
			
			////////////////////////////////////////////
            //// Process HVC-RA -> HVC-RA interactions
            ////////////////////////////////////////////
			if ( some_RA_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons
				for (int i = 0; i < N_RA_local; i++)
					HVCRA_local[i].raiseExcWeight(update_Ge_RA_global[Id_RA_local[i]]); // update conductance

				// update conductance arrays and fired indicators
				some_RA_exc_conductance_was_updated_global = 0;
				some_RA_exc_conductance_was_updated_local = 0;
				
				std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
			}	
            
            network_time += NETWORK_UPDATE_FREQUENCY;
        }
    }
}

void NetworkGrowthSimulator::trial_no_stdp(double training_kick_time)
{
	// indicators for conductance updates and bursting
	int some_RA_inh_conductance_was_updated_global = 0;
	int some_RA_inh_conductance_was_updated_local = 0;
	
	int some_RA_exc_conductance_was_updated_global = 0;
	int some_RA_exc_conductance_was_updated_local = 0;
	
	int some_I_exc_conductance_was_updated_global = 0;
	int some_I_exc_conductance_was_updated_local = 0;
	
	 // arrays with updates for conductances
    std::vector<double> update_Ge_RA_local(N_RA);
	std::vector<double> update_Ge_RA_global(N_RA);

	std::vector<double> update_Gi_RA_local(N_RA);
	std::vector<double> update_Gi_RA_global(N_RA);

	std::vector<double> update_Ge_I_local(N_I);
	std::vector<double> update_Ge_I_global(N_I);
	
	std::vector<int> RA_neurons_bursted_local; // local array of HVC-RA neurons that bursted
	std::vector<int> RA_neurons_bursted_global; // global array of HVC-RA neurons that bursted
	
	
	if (MPI_rank == 0)
	{
		std::cout << "training_kick_time = " << training_kick_time << std::endl;
    }
    
    
    bool training_excited = false; // indicator that training neurons were already excited
    
    double internal_time = 0;
    double network_time = NETWORK_UPDATE_FREQUENCY;

   

    //printf("network time = %f\n", network_time);

    
	// initialize update arrays and fired indicators
	std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
	
	
    // evolve dynamics
    int num_steps = static_cast<int>(TRIAL_DURATION / TIMESTEP);
    
    for (int t = 1; t < num_steps; t++)
	{
		internal_time += TIMESTEP;
	
		if ( ( !training_excited ) && ( internal_time >= training_kick_time ) )
		{
			for (int i = 0; i < N_TR; i++)
			{
				int rank;
				int shift;
				
				this->get_neuronRA_location(training_neurons[i], &rank, &shift);
				
				if (MPI_rank == rank)
					HVCRA_local[shift].raiseE(G_TRAINING_KICK);
			}
			
			training_excited = true;
		}
		
		//////////////////////
		// HVC-RA neurons
		//////////////////////
		
		for (int i = 0; i < N_RA_local; i++)
		{
            // Debraband step
            HVCRA_local[i].Debraband_step_no_target_update();
            
           
			///////////////////////////////
			// HVC-RA somatic spike
			///////////////////////////////
			// if some neuron produced somatic spike
            if (HVCRA_local[i].get_fired_soma())
            {
                spikes_in_trial_soma_local[i].push_back(internal_time);
                
				// update delivery queues
				
				// for inhibitory neurons
				// loop over all inhibitory targets of fired neurons
				size_t num_I_targets = syn_ID_RA_I_local[i].size();
				
				for (size_t j = 0; j < num_I_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_RA_I_local[i][j];
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_RA_I[i].empty() )
						delivery_queue_RA_I[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_RA_I[i].begin(), delivery_queue_RA_I[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_RA_I[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
					//std::cout << "neuron " << Id_RA_local[i] << " spike time = " << internal_time << " axonal delay = " << axonal_delays_RA_I[Id_RA_local[i]][syn_ID] << " delivery time RA to I: " << delivery_queue_RA_I[i].back().first << " delivery target id: " << syn_ID_RA_I_local[i][delivery_queue_RA_I[i].back().second] << std::endl;
				
				// deliver spikes to active targets
				size_t num_RA_targets = active_synapses_local[i].size();
				
				for (size_t j = 0; j < num_RA_targets; j++)
				{	
					int target_id = active_synapses_local[i][j];
					double delivery_time = internal_time + axonal_delays_RA_RA_local[i][target_id];
				
					// if queue is empty, just add item to the queue
					if ( delivery_queue_RA_RA_soma[i].empty() )
						delivery_queue_RA_RA_soma[i].push_back(std::pair<double,int>(delivery_time, target_id));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_RA_RA_soma[i].begin(), delivery_queue_RA_RA_soma[i].end(), std::pair<double,int>(delivery_time, target_id));
						delivery_queue_RA_RA_soma[i].insert(it, std::pair<double,int>(delivery_time, target_id));
					}
				}
                
            } // end if local HVC-RA neuron spiked  
            
             //////////////////////////////////
			// HVC-RA -> HVC-I delivery queue
			//////////////////////////////////
            // check if somatic spike was delivered to some interneuron
            // loop through the delivery queue to check if current time exceeds the spike delivery time
            std::vector<std::pair<double,int>>::iterator it = delivery_queue_RA_I[i].begin();
			
            for (; it != delivery_queue_RA_I[i].end(); it++)
            {
				if (internal_time >= it->first)
				{
					
					some_I_exc_conductance_was_updated_local = 1;
					
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_RA_I_local[i][pos_in_local_target_array];
					
					update_Ge_I_local[target_id] += weights_RA_I_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(RA) neuron " << Id_RA_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to HVC(I) neuron " << syn_ID_RA_I_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_I[i].erase(delivery_queue_RA_I[i].begin(), it);
			
			//////////////////////////////////////////////////
			// HVC-RA -> HVC-RA somatic spike delivery queue
			//////////////////////////////////////////////////
			// check if somatic spike was delivered to some HVC(RA) neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
            it = delivery_queue_RA_RA_soma[i].begin();
            
            for (; it != delivery_queue_RA_RA_soma[i].end(); it++)
            {
				double delivery_time = it->first;
				
				if ( internal_time >= delivery_time )
				{
					int target_id = it->second;
					
					
					some_RA_exc_conductance_was_updated_local = 1;
					update_Ge_RA_local[target_id] += weights_RA_RA_local[i][target_id];
					
					
				} // end if internal_time >= delivery spike time 
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_RA_soma[i].erase(delivery_queue_RA_RA_soma[i].begin(), it);
			
			
			///////////////////////////////
			// HVC-RA dendritic spike
			///////////////////////////////
			//if some neuron produced dendritic spike, store this neuron in array
			
			if ( HVCRA_local[i].get_fired_dend() )
			{
				//~ // show spikes of pool neurons only
				//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
				//~ 
				//~ if (p.first == p.second)
					//~ std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
				//~ 
				//~ 
				spikes_in_trial_dend_local[i].push_back(internal_time); // dend spike time relative to trial onset
				
				//std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
			} // end if get_fired_dend
	
            
		} // end for i = 0 -> N_RA_local (loop through all local HVC-RA neurons)
		
		//////////////////////
		// HVC-I neurons
		//////////////////////
		for (int i = 0; i < N_I_local; i++)
		{
			HVCI_local[i].DP8_step_no_target_update();
			
			
			///////////////////////////////
			// HVC-I spike
			///////////////////////////////
			//  if some I neuron spikes, update delivery queue
			if (HVCI_local[i].get_fired())
			{
				//printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
				spikes_in_trial_interneuron_local[i].push_back(internal_time);

				size_t num_RA_targets = syn_ID_I_RA_local[i].size();
				// loop over all targets of fired neurons
				for (size_t j = 0; j < num_RA_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_I_RA_local[i][j];
					
					//std::cout << "HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << syn_ID_I_RA_local[i][j] << " spike time " << internal_time << " delivery time = " << delivery_time << std::endl; 
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_I_RA[i].empty() )
						delivery_queue_I_RA[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_I_RA[i].begin(), delivery_queue_I_RA[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_I_RA[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
			}
			
			//////////////////////////////////
			// HVC-I -> HVC-RA delivery queue
			//////////////////////////////////
			// check if interneuron spike was delivered to some HVC-RA neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
			std::vector<std::pair<double,int>>::iterator it = delivery_queue_I_RA[i].begin();
			
			for (; it != delivery_queue_I_RA[i].end(); it++)
			{
				if (internal_time >= it->first)
				{
					some_RA_inh_conductance_was_updated_local = 1;
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_I_RA_local[i][pos_in_local_target_array];
					
					//std::cout << "Delivered spike HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << target_id << " delivered_time " << internal_time << " delivery time in queue = " << it->first << std::endl; 
					
					
					update_Gi_RA_local[target_id] += weights_I_RA_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(I) neuron " << Id_I_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to neuron " << syn_ID_I_RA_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_I_RA[i].erase(delivery_queue_I_RA[i].begin(), it);
		
		}  // end for i = 0 -> N_I_local (loop through all local HVC-I neurons)
		
		/////////////////////////////
		// Network Synchronization //
		/////////////////////////////
		if ( ( internal_time > network_time ) || ( t == num_steps-1 ) )
        {
            //if (MPI_rank == 0)
            //{
            //    std::cout << "internal time = " << internal_time << " > network_time = " << network_time << std::endl; 
            //}

            MPI_Allreduce(&some_RA_exc_conductance_was_updated_local, &some_RA_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_RA_inh_conductance_was_updated_local, &some_RA_inh_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_I_exc_conductance_was_updated_local, &some_I_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
			
            ////////////////////////////////////////////
            //// Process HVC-I -> HVC-RA interactions
            ////////////////////////////////////////////
            if ( some_RA_inh_conductance_was_updated_global > 0 )
            {
            // sum update array and send to all processes

                MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                for (int i = 0; i < N_RA_local; i++)
                {
					//std::cout << "HVC-RA neuron " << Id_RA_local[i] << " raised inhibitory conductance by " << update_Gi_RA_global[Id_RA_local[i]] << std::endl;
                    HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
            	}
				// update conductance arrays and delivered indicators
				some_RA_inh_conductance_was_updated_global = 0;
				some_RA_inh_conductance_was_updated_local = 0;
				
				std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            }
            
            ////////////////////////////////////////////
            //// Process HVC-RA -> HVC-I interactions
            ////////////////////////////////////////////
			if ( some_I_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons

				for (int i = 0; i < N_I_local; i++)
					HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);

				// update conductance arrays and fired indicators
				some_I_exc_conductance_was_updated_local = 0;
				some_I_exc_conductance_was_updated_global = 0;

				std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
			}
			
			////////////////////////////////////////////
            //// Process HVC-RA -> HVC-RA interactions
            ////////////////////////////////////////////
			if ( some_RA_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons
				for (int i = 0; i < N_RA_local; i++)
					HVCRA_local[i].raiseExcWeight(update_Ge_RA_global[Id_RA_local[i]]); // update conductance

				// update conductance arrays and fired indicators
				some_RA_exc_conductance_was_updated_global = 0;
				some_RA_exc_conductance_was_updated_local = 0;
				
				std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
			}	
            
            network_time += NETWORK_UPDATE_FREQUENCY;
        }
    }
}

void NetworkGrowthSimulator::trial_soma_pre_dend_post_stdp_delays(bool training)
{
	// indicators for conductance updates and bursting
	int some_RA_inh_conductance_was_updated_global = 0;
	int some_RA_inh_conductance_was_updated_local = 0;
	
	int some_RA_exc_conductance_was_updated_global = 0;
	int some_RA_exc_conductance_was_updated_local = 0;
	
	int some_I_exc_conductance_was_updated_global = 0;
	int some_I_exc_conductance_was_updated_local = 0;
	
	int some_RA_neuron_bursted_global = 0;
	int some_RA_neuron_bursted_local = 0;
	
	 // arrays with updates for conductances
    std::vector<double> update_Ge_RA_local(N_RA);
	std::vector<double> update_Ge_RA_global(N_RA);

	std::vector<double> update_Gi_RA_local(N_RA);
	std::vector<double> update_Gi_RA_global(N_RA);

	std::vector<double> update_Ge_I_local(N_I);
	std::vector<double> update_Ge_I_global(N_I);
	
	std::vector<int> RA_neurons_bursted_local; // local array of HVC-RA neurons that bursted
	std::vector<int> RA_neurons_bursted_global; // global array of HVC-RA neurons that bursted
	
	// sample training innervation time and send to all processes
    double training_kick_time;
    
	if (MPI_rank == 0)
	{
		training_kick_time = WAITING_TIME + noise_generator.random(TRIAL_DURATION - 2*WAITING_TIME);
		std::cout << "training_kick_time = " << training_kick_time << std::endl;
    }
    MPI_Bcast(&training_kick_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    bool training_excited = false; // indicator that training neurons were already excited
    
    double internal_time = 0;
    double network_time = NETWORK_UPDATE_FREQUENCY;

   

    //printf("network time = %f\n", network_time);

    
	// initialize update arrays and fired indicators
	std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
	
	
    // evolve dynamics
    int num_steps = static_cast<int>(TRIAL_DURATION / TIMESTEP);
    
    for (int t = 1; t < num_steps; t++)
	{
		internal_time += TIMESTEP;
	
		if ( ( training ) && ( !training_excited ) && ( internal_time >= training_kick_time ) )
		{
			for (int i = 0; i < N_TR; i++)
			{
				int rank;
				int shift;
				
				this->get_neuronRA_location(training_neurons[i], &rank, &shift);
				
				if (MPI_rank == rank)
					HVCRA_local[shift].raiseE(G_TRAINING_KICK);
			}
			
			training_excited = true;
		}
		
		//////////////////////
		// HVC-RA neurons
		//////////////////////
		
		for (int i = 0; i < N_RA_local; i++)
		{
            // Debraband step
            HVCRA_local[i].Debraband_step_no_target_update();
            
           
			///////////////////////////////
			// HVC-RA somatic spike
			///////////////////////////////
			// if some neuron produced somatic spike
            if (HVCRA_local[i].get_fired_soma())
            {
                spikes_in_trial_soma_local[i].push_back(internal_time);
                
				// update delivery queues
				
				// for inhibitory neurons
				// loop over all inhibitory targets of fired neurons
				size_t num_I_targets = syn_ID_RA_I_local[i].size();
				
				for (size_t j = 0; j < num_I_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_RA_I_local[i][j];
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_RA_I[i].empty() )
						delivery_queue_RA_I[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_RA_I[i].begin(), delivery_queue_RA_I[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_RA_I[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
					//std::cout << "neuron " << Id_RA_local[i] << " spike time = " << internal_time << " axonal delay = " << axonal_delays_RA_I[Id_RA_local[i]][syn_ID] << " delivery time RA to I: " << delivery_queue_RA_I[i].back().first << " delivery target id: " << syn_ID_RA_I_local[i][delivery_queue_RA_I[i].back().second] << std::endl;
				
				
				// if neuron is saturated, deliver spikes only to supersynaptic targets
				if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
				{
					for (size_t j = 0; j < supersynapses_local[i].size(); j++)
					{
						int target_id = supersynapses_local[i][j];
						
						double delivery_time = internal_time + axonal_delays_RA_RA_local[i][target_id];
						
						// if queue is empty, just add item to the queue
						if ( delivery_queue_RA_RA_soma[i].empty() )
							delivery_queue_RA_RA_soma[i].push_back(std::pair<double,int>(delivery_time, target_id));	
						// otherwise add item so that queue is sorted
						else
						{
							auto it = std::upper_bound(delivery_queue_RA_RA_soma[i].begin(), delivery_queue_RA_RA_soma[i].end(), std::pair<double,int>(delivery_time, target_id));
							delivery_queue_RA_RA_soma[i].insert(it, std::pair<double,int>(delivery_time, target_id));
						}
					}
				}
				else // deliver spikes to everyone except itself
				{	
					for (int j = 0; j < N_RA; j++)
					{
						if ( j != Id_RA_local[i])
						{
						
							double delivery_time = internal_time + axonal_delays_RA_RA_local[i][j];
						
							// if queue is empty, just add item to the queue
							if ( delivery_queue_RA_RA_soma[i].empty() )
								delivery_queue_RA_RA_soma[i].push_back(std::pair<double,int>(delivery_time, j));	
							// otherwise add item so that queue is sorted
							else
							{
								auto it = std::upper_bound(delivery_queue_RA_RA_soma[i].begin(), delivery_queue_RA_RA_soma[i].end(), std::pair<double,int>(delivery_time, j));
								delivery_queue_RA_RA_soma[i].insert(it, std::pair<double,int>(delivery_time, j));
							}
						}
					}
				}
                
            } // end if local HVC-RA neuron spiked  
            
             //////////////////////////////////
			// HVC-RA -> HVC-I delivery queue
			//////////////////////////////////
            // check if somatic spike was delivered to some interneuron
            // loop through the delivery queue to check if current time exceeds the spike delivery time
            std::vector<std::pair<double,int>>::iterator it = delivery_queue_RA_I[i].begin();
			
            for (; it != delivery_queue_RA_I[i].end(); it++)
            {
				if (internal_time >= it->first)
				{
					
					some_I_exc_conductance_was_updated_local = 1;
					
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_RA_I_local[i][pos_in_local_target_array];
					
					update_Ge_I_local[target_id] += weights_RA_I_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(RA) neuron " << Id_RA_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to HVC(I) neuron " << syn_ID_RA_I_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_I[i].erase(delivery_queue_RA_I[i].begin(), it);
			
			//////////////////////////////////////////////////
			// HVC-RA -> HVC-RA somatic spike delivery queue
			//////////////////////////////////////////////////
			// check if somatic spike was delivered to some HVC(RA) neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
            it = delivery_queue_RA_RA_soma[i].begin();
            
            for (; it != delivery_queue_RA_RA_soma[i].end(); it++)
            {
				double delivery_time = it->first;
				
				if ( internal_time >= delivery_time )
				{
					int target_id = it->second;
					
					delivered_spike_times[i][target_id].push_back(delivery_time);
					
					// update conductance of target if synapse is active
					if ( active_indicators_local[i][target_id] == 1 )
					{
						some_RA_exc_conductance_was_updated_local = 1;
						update_Ge_RA_local[target_id] += weights_RA_RA_local[i][target_id];
					}
					
					//////////////
					//// LTD
					//////////////
					// do LTD on dendritic spike time of target neuron
					 // if neuron is saturated apply LTD only if target is among super synapses
					if (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss)
					{
						if ( supersynapses_indicators_local[i][target_id] == 1 )
						{	
							// find previous dendritic spikes that are not too old
							if ( !previous_dendritic_spike_times_global[target_id].empty() )
							{
								std::vector<double>::iterator it_relevant_dend_spikes = std::lower_bound(previous_dendritic_spike_times_global[target_id].begin(),
																									previous_dendritic_spike_times_global[target_id].end(),
																									internal_time - STDP_WINDOW);
																					
								for (auto it_dend_spike = it_relevant_dend_spikes; it_dend_spike != previous_dendritic_spike_times_global[target_id].end(); it_dend_spike++)
								{
									double dt = *it_dend_spike - delivery_time;

									//std::cout << "dt in saturated LTD = " << dt << std::endl;

									
									double w_before = weights_RA_RA_local[i][target_id];     
						   
									LTD(weights_RA_RA_local[i][target_id], dt);
									
									
									//~ std::cout << "LTD from saturated " << Id_RA_local[i] << " -> " << supersynapse_id
											  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][supersynapse_id] - w_before
											  //~ << std::endl;
											  
									//printf("LTD from saturated neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n",
									//            Id_RA_local[i], supersynapse_id, spikes_in_trial_soma_local[i].back(), spike_times_dend_global[supersynapse_id],
									//            dt, weights_local[i][supersynapse_id] - w);
									
									update_synapse(i, target_id);	
									
								}   
							}
						}

						// if some supersynapse desaturated, update all synapses
						if (static_cast<int>(supersynapses_local[i].size()) < synaptic_params.Nss)
							for (int j = 0; j < N_RA; j++)
								this->update_synapse(i, j);
					}
					// if not saturated apply LTD rule 
					else
					{
						// find previous dendritic spikes that are not too old
						if ( !previous_dendritic_spike_times_global[target_id].empty() )
						{
							std::vector<double>::iterator it_relevant_dend_spikes = std::lower_bound(previous_dendritic_spike_times_global[target_id].begin(),
																								previous_dendritic_spike_times_global[target_id].end(),
																								internal_time - STDP_WINDOW);
																				
							for (auto it_dend_spike = it_relevant_dend_spikes; it_dend_spike != previous_dendritic_spike_times_global[target_id].end(); it_dend_spike++)
							{
								double dt = *it_dend_spike - delivery_time;

								//std::cout << "dt in saturated LTD = " << dt << std::endl;

								
								double w_before = weights_RA_RA_local[i][target_id];     
					   
								LTD(weights_RA_RA_local[i][target_id], dt);
								
								//~ // show LTD of pool neurons only
								//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
								//~ 
								//~ if (p.first == p.second)
								//~ {
									//~ std::cout << "LTD from  " << Id_RA_local[i] << " -> " << target_id
										  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][target_id] - w_before
										  //~ << " w_after = " << weights_RA_RA_local[i][target_id] << std::endl;
								//~ }
							//~ 
								update_synapse(i, target_id);	
								
							}   
						}
					}
				} // end if internal_time >= delivery spike time 
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_RA_RA_soma[i].erase(delivery_queue_RA_RA_soma[i].begin(), it);
			
			
			///////////////////////////////
			// HVC-RA dendritic spike
			///////////////////////////////
			//if some neuron produced dendritic spike, store this neuron in array
			
			if ( HVCRA_local[i].get_fired_dend() )
			{
				// show spikes of pool neurons only
				auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
				
				if (p.first == p.second)
					std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
				
				
				spikes_in_trial_dend_local[i].push_back(internal_time); // dend spike time relative to trial onset
				
				RA_neurons_bursted_local.push_back(Id_RA_local[i]);
				
				some_RA_neuron_bursted_local = 1;
				
				//std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
			} // end if get_fired_dend
	
            
		} // end for i = 0 -> N_RA_local (loop through all local HVC-RA neurons)
		
		//////////////////////
		// HVC-I neurons
		//////////////////////
		for (int i = 0; i < N_I_local; i++)
		{
			HVCI_local[i].DP8_step_no_target_update();
			
			
			///////////////////////////////
			// HVC-I spike
			///////////////////////////////
			//  if some I neuron spikes, update delivery queue
			if (HVCI_local[i].get_fired())
			{
				//printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
				spikes_in_trial_interneuron_local[i].push_back(internal_time);

				size_t num_RA_targets = syn_ID_I_RA_local[i].size();
				// loop over all targets of fired neurons
				for (size_t j = 0; j < num_RA_targets; j++)
				{
					double delivery_time = internal_time + axonal_delays_I_RA_local[i][j];
					
					//std::cout << "HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << syn_ID_I_RA_local[i][j] << " spike time " << internal_time << " delivery time = " << delivery_time << std::endl; 
					
					// if queue is empty, just add item to the queue
					if ( delivery_queue_I_RA[i].empty() )
						delivery_queue_I_RA[i].push_back(std::pair<double,int>(delivery_time, j));	
					// otherwise add item so that queue is sorted
					else
					{
						auto it = std::upper_bound(delivery_queue_I_RA[i].begin(), delivery_queue_I_RA[i].end(), std::pair<double,int>(delivery_time, j));
						delivery_queue_I_RA[i].insert(it, std::pair<double,int>(delivery_time, j));
					}
				}
			}
			
			//////////////////////////////////
			// HVC-I -> HVC-RA delivery queue
			//////////////////////////////////
			// check if interneuron spike was delivered to some HVC-RA neuron
			// loop through the delivery queue to check if current time exceeds the spike delivery time
			std::vector<std::pair<double,int>>::iterator it = delivery_queue_I_RA[i].begin();
			
			for (; it != delivery_queue_I_RA[i].end(); it++)
			{
				if (internal_time >= it->first)
				{
					some_RA_inh_conductance_was_updated_local = 1;
					int pos_in_local_target_array = it->second;
					int target_id = syn_ID_I_RA_local[i][pos_in_local_target_array];
					
					//std::cout << "Delivered spike HVC-I -> HVC-RA: " << Id_I_local[i] << " -> " << target_id << " delivered_time " << internal_time << " delivery time in queue = " << it->first << std::endl; 
					
					
					update_Gi_RA_local[target_id] += weights_I_RA_local[i][pos_in_local_target_array];
					
					//~ std::cout << "HVC(I) neuron " << Id_I_local[i] << " at time " << internal_time 
							  //~ << " delivered spike with delivery time " << it->first 
							  //~ << " to neuron " << syn_ID_I_RA_local[i][pos_in_local_target_array] << std::endl;
				}
				else
					break;
			}
			
			// delete all processed spikes
			delivery_queue_I_RA[i].erase(delivery_queue_I_RA[i].begin(), it);
		
		}  // end for i = 0 -> N_I_local (loop through all local HVC-I neurons)
		
		/////////////////////////////
		// Network Synchronization //
		/////////////////////////////
		if ( ( internal_time > network_time ) || ( t == num_steps-1 ) )
        {
            //if (MPI_rank == 0)
            //{
            //    std::cout << "internal time = " << internal_time << " > network_time = " << network_time << std::endl; 
            //}

            MPI_Allreduce(&some_RA_exc_conductance_was_updated_local, &some_RA_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_RA_inh_conductance_was_updated_local, &some_RA_inh_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_I_exc_conductance_was_updated_local, &some_I_exc_conductance_was_updated_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            MPI_Allreduce(&some_RA_neuron_bursted_local, &some_RA_neuron_bursted_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            
        
			
            ////////////////////////////////////////////
            //// Process HVC-I -> HVC-RA interactions
            ////////////////////////////////////////////
            if ( some_RA_inh_conductance_was_updated_global > 0 )
            {
            // sum update array and send to all processes

                MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                for (int i = 0; i < N_RA_local; i++)
                {
					//std::cout << "HVC-RA neuron " << Id_RA_local[i] << " raised inhibitory conductance by " << update_Gi_RA_global[Id_RA_local[i]] << std::endl;
                    HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
            	}
				// update conductance arrays and delivered indicators
				some_RA_inh_conductance_was_updated_global = 0;
				some_RA_inh_conductance_was_updated_local = 0;
				
				std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            }
            
            ////////////////////////////////////////////
            //// Process HVC-RA -> HVC-I interactions
            ////////////////////////////////////////////
			if ( some_I_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons

				for (int i = 0; i < N_I_local; i++)
					HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);

				// update conductance arrays and fired indicators
				some_I_exc_conductance_was_updated_local = 0;
				some_I_exc_conductance_was_updated_global = 0;

				std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
			}
			
			////////////////////////////////////////////
            //// Process HVC-RA -> HVC-RA interactions
            ////////////////////////////////////////////
			if ( some_RA_exc_conductance_was_updated_global > 0 )
			{
				// sum all update arrays and send to all processes

				MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// now update excitatory conductances of all neurons
				for (int i = 0; i < N_RA_local; i++)
					HVCRA_local[i].raiseExcWeight(update_Ge_RA_global[Id_RA_local[i]]); // update conductance

				// update conductance arrays and fired indicators
				some_RA_exc_conductance_was_updated_global = 0;
				some_RA_exc_conductance_was_updated_local = 0;
				
				std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
			}	

			////////////////////////////////////////////
            //// Process HVC-RA dendritic spikes
            ////////////////////////////////////////////
			if ( some_RA_neuron_bursted_global > 0 )
            {
				// update conductance arrays and fired indicators
            	some_RA_neuron_bursted_local = 0;
	        	some_RA_neuron_bursted_global = 0;
	        	
	        	// gather all bursted HVC-RA neurons
				this->gather_spiked_or_bursted_neurons(RA_neurons_bursted_local, RA_neurons_bursted_global);
				
				// add dendritic spikes for bursted neurons
				for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
					previous_dendritic_spike_times_global[RA_neurons_bursted_global[j]].push_back(internal_time);
				
				////////////////////////////////////////////////////////////////////////////////////
				// LTP or LTD of previously delivered somatic spike times on bursted HVC-RA neuron
				////////////////////////////////////////////////////////////////////////////////////
				
				for (int i = 0; i < N_RA_local; i++)
                {
					int presyn_ID = Id_RA_local[i]; // real id of presynaptic neuron
					
                    // if neuron is saturated apply LTP only if spiked neurons are among supersynapse targets
                    if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
                    {
                        for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_bursted_global[j]; // id of postsynaptic neuron
                            
                            // do not allow self-to-self synapses to emerge
							if ( presyn_ID != postsyn_ID )
							{
								std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
											supersynapses_local[i].end(), postsyn_ID );

								if ( pos!=supersynapses_local[i].end() )
								{
									// find previous somatic spikes that are not too old
									if ( !delivered_spike_times[i][postsyn_ID].empty() )
									{
										std::vector<double>::iterator it_relevant_spikes = std::lower_bound(delivered_spike_times[i][postsyn_ID].begin(),
																											delivered_spike_times[i][postsyn_ID].end(),
																											internal_time - STDP_WINDOW);
																							
										for (auto it = it_relevant_spikes; it != delivered_spike_times[i][postsyn_ID].end(); it++)
										{
											double dt = internal_time - *it;
									
										
											//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
											
											if (dt <= synaptic_params.T_0)
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												LTD(weights_RA_RA_local[i][postsyn_ID], dt);
												
												//~ std::cout   << "LTD from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
											}
											else
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												
												LTP(weights_RA_RA_local[i][postsyn_ID], dt);
												
												//~ std::cout   << "LTP from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
												
											}	
											//double w = weights_local[i][fired_ID];
											
											//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
											 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
											 //           dt, weights_local[i][fired_ID] - w);
											
											update_synapse(i, postsyn_ID);
											
										}
									}
								}
							}
                        }

                    }
                    // if not saturated apply LTP for all dendritic spikes
                    else
                    {
                        for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_bursted_global[j];
                            // don't allow self-to-self connections
                            if ( postsyn_ID != presyn_ID )
                            {
                                // find previous somatic spikes that are not too old
								if ( !delivered_spike_times[i][postsyn_ID].empty() )
								{
									std::vector<double>::iterator it_relevant_spikes = std::lower_bound(delivered_spike_times[i][postsyn_ID].begin(),
																										delivered_spike_times[i][postsyn_ID].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it = it_relevant_spikes; it != delivered_spike_times[i][postsyn_ID].end(); it++)
									{
										double dt = internal_time - *it;
								

										//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
										
									    if (dt <= synaptic_params.T_0)
									    {
											double w_before = weights_RA_RA_local[i][postsyn_ID];
												
											LTD(weights_RA_RA_local[i][postsyn_ID], dt);
											
											//~ // show LTD on pool neurons only
											//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											//~ 
											//~ if (p.first == p.second)
											//~ {
												//~ std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											//~ }
											
											//~ std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
											
										}
										else
										{
											double w_before = weights_RA_RA_local[i][postsyn_ID];
											
											LTP(weights_RA_RA_local[i][postsyn_ID], dt);
											
											//~ // show LTP on pool neurons only
											//~ auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											//~ 
											//~ if (p.first == p.second)
											//~ {
												//~ std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << " w_after = " << weights_RA_RA_local[i][postsyn_ID] << std::endl;
											//~ }
											
											
											//~ std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
										}	
										//double w = weights_local[i][fired_ID];
										
										//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
										 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
										 //           dt, weights_local[i][fired_ID] - w);
										
										update_synapse(i, postsyn_ID);
									
									}
									
									
								}
                                
                            }
                        }
                   }

                } // end for i -> N_RA_local

				       

               // check if we need axon remodeling

                for (int i = 0; i < N_RA_local; i++)
                {
                    if ( (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss) && (remodeled_local[i] == 0) )
                    {
					    this->axon_remodeling(i);

				    }
                }
            
				// clear bursts
				RA_neurons_bursted_local.clear();
				RA_neurons_bursted_global.clear();
				
            } // end if some HVC-RA bursted
            
            
            network_time += NETWORK_UPDATE_FREQUENCY;
        }
    }
    
    this->potentiation_decay();
    //printf("After potentiation decay")
    this->update_all_synapses();
	
	// update rate
	for (int i = 0; i < N_RA_local; i++)
	{
		num_spikes_in_recent_trials_local[i].push_front(static_cast<int>(spikes_in_trial_soma_local[i].size()));
        
        //firing_rate_short_local[i] = std::accumulate(num_spikes_in_recent_trials[i].begin(), num_spikes_in_recent_trials[i].begin() + RATE_WINDOW_SHORT, 0.0)
          //                          / static_cast<double>(RATE_WINDOW_SHORT);

		// calculate firing rate in large window:
        firing_rate_long_local[i] = std::accumulate(num_spikes_in_recent_trials_local[i].begin(), num_spikes_in_recent_trials_local[i].end(), 0.0)
                                                                                            / static_cast<double>(maturation_params.RATE_WINDOW_LONG);
                                                                                            
		//if (std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].end(), 0.0) >= 1)
		//{
		//	std::cout << "Recent bursts of neuron " << Id_RA_local[i] << " :" << std::endl;
			
		//	for (size_t j = 0; j < num_bursts_in_recent_trials[i].size(); j++)
		//		std::cout << num_bursts_in_recent_trials[i][j] << "\t";
			
		//	std::cout << std::endl;
		//}
	}
}

void NetworkGrowthSimulator::trial_soma_pre_dend_post_stdp_no_delays(bool training)
{
	// spike indicators
	int some_RA_neuron_spiked_local = 0;
	int some_RA_neuron_bursted_local = 0;
	int some_I_neuron_fired_local = 0;
	
	int some_RA_neuron_spiked_global = 0;
	int some_RA_neuron_bursted_global = 0;
	int some_I_neuron_fired_global = 0;

    std::vector<double> spike_times_fired_dend_local;
    
    std::vector<int> RA_neurons_spiked_global;
    std::vector<int> RA_neurons_spiked_local;
    
    std::vector<int> RA_neurons_bursted_global;
    std::vector<int> RA_neurons_bursted_local;
    
    
    // arrays with updates for conductances
    std::vector<double> update_Ge_RA_local(N_RA);
	std::vector<double> update_Ge_RA_global(N_RA);

	std::vector<double> update_Gi_RA_local(N_RA);
	std::vector<double> update_Gi_RA_global(N_RA);

	std::vector<double> update_Ge_I_local(N_I);
	std::vector<double> update_Ge_I_global(N_I);

	// sample training innervation time and send to all processes
    double training_kick_time;
    
	if (MPI_rank == 0)
	{
		training_kick_time = WAITING_TIME + noise_generator.random(TRIAL_DURATION - 2*WAITING_TIME);
		std::cout << "training_kick_time = " << training_kick_time << std::endl;
    }
    MPI_Bcast(&training_kick_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    bool training_excited = false; // indicator that training neurons were already excited
    
    double internal_time = 0;
    double network_time = NETWORK_UPDATE_FREQUENCY;

   

    //printf("network time = %f\n", network_time);

    
	// initialize update arrays and fired indicators
	std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
	
	int num_steps = static_cast<int>(TRIAL_DURATION / TIMESTEP);
	
    // evolve dynamics
    for (int t = 1; t < num_steps; t++)
	{
		internal_time += TIMESTEP;
	
		if ( ( training ) && ( !training_excited ) && ( internal_time >= training_kick_time ) )
		{
			for (int i = 0; i < N_TR; i++)
			{
				int rank;
				int shift;
				
				this->get_neuronRA_location(training_neurons[i], &rank, &shift);
				
				if (MPI_rank == rank)
					HVCRA_local[shift].raiseE(G_TRAINING_KICK);
			}
			
			training_excited = true;
		}
		
		for (int i = 0; i < N_RA_local; i++)
		{
            // Debraband step
            HVCRA_local[i].Debraband_step_no_target_update();
            
            
            
            // if some neuron produced somatic spike, do LTD for all previous dendritic spikes
            if ( HVCRA_local[i].get_fired_soma() )
            {
				
				// show spikes of pool neurons only
                auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
                
                if (p.first == p.second)
					std::cout << "HVC-RA " << Id_RA_local[i] << " soma spike at " << internal_time << std::endl;
                	
				//std::cout << "HVC-RA " << Id_RA_local[i] << " soma spike at " << internal_time << std::endl;
                	
				
                spikes_in_trial_soma_local[i].push_back(internal_time);
				previous_somatic_spike_times_local[i].push_back(internal_time);
				
				//RA_neurons_spiked_local.push_back(Id_RA_local[i]);
                
				some_RA_neuron_spiked_local = 1;
				
				// loop over all inhibitory targets of fired neurons
				size_t num_I_targets = syn_ID_RA_I_local[i].size();
				for (size_t j = 0; j < num_I_targets; j++)
				{
					int syn_ID = syn_ID_RA_I_local[i][j];
					update_Ge_I_local[syn_ID] += weights_RA_I_local[i][j];

				}
				
				
				// loop over all excitatory targets
				size_t num_RA_targets = active_synapses_local[i].size();
				for (size_t j = 0; j < num_RA_targets; j++)
				{
					int syn_ID = active_synapses_local[i][j];
					update_Ge_RA_local[syn_ID] += weights_RA_RA_local[i][syn_ID];
				}

				
                //for (int j = 0; j < last_soma_spikes_local[i].size(); j++)
                //{
                  //  printf("My rank = %d; Soma spikes of neuron %d:  spike_time = %f\n", MPI_rank, Id_RA_local[i],
                    //    last_soma_spikes_local[i][j]);

                //}

				//std::cout << "HVC-RA " << Id_RA_local[i] << " soma spike at " << internal_time << std::endl;
                
                 // if neuron is saturated apply LTD only to supersynapses
                if (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss)
                {
                    for (size_t k = 0; k < supersynapses_local[i].size(); k++)
                    {
                        int supersynapse_id = supersynapses_local[i][k];
                        
                        if ( Id_RA_local[i] != supersynapse_id )
                        {
							// find previous dendritic spikes that are not too old
							if ( !previous_dendritic_spike_times_global[supersynapse_id].empty() )
							{
								std::vector<double>::iterator it_relevant_dend_spikes = std::lower_bound(previous_dendritic_spike_times_global[supersynapse_id].begin(),
																									previous_dendritic_spike_times_global[supersynapse_id].end(),
																									internal_time - STDP_WINDOW);
																					
								for (auto it = it_relevant_dend_spikes; it != previous_dendritic_spike_times_global[supersynapse_id].end(); it++)
								{
									double dt = *it - internal_time;

									//std::cout << "dt in saturated LTD = " << dt << std::endl;

									
									double w_before = weights_RA_RA_local[i][supersynapse_id];     
						   
									LTD(weights_RA_RA_local[i][supersynapse_id], dt);
									
									
									//~ std::cout << "LTD from saturated " << Id_RA_local[i] << " -> " << supersynapse_id
											  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][supersynapse_id] - w_before
											  //~ << std::endl;
											  
									//printf("LTD from saturated neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n",
									//            Id_RA_local[i], supersynapse_id, spikes_in_trial_soma_local[i].back(), spike_times_dend_global[supersynapse_id],
									//            dt, weights_local[i][supersynapse_id] - w);
									
									update_synapse(i, supersynapse_id);	
									
								}
							}
                           
                        }
                    }

					// if some supersynapse desaturated, update all synapses
                	if (static_cast<int>(supersynapses_local[i].size()) < synaptic_params.Nss)
						for (int j = 0; j < N_RA; j++)
							this->update_synapse(i, j);
                }
                // if not saturated apply LTD rule to recent dendritic spikes of all neurons 
                else
                {
                    for (int j = 0; j < N_RA; j++)
                    {
                        if ( Id_RA_local[i] != j )
                        {
							// find previous dendritic spikes that are not too old
							if ( !previous_dendritic_spike_times_global[j].empty() )
							{
								std::vector<double>::iterator it_relevant_spikes = std::lower_bound(previous_dendritic_spike_times_global[j].begin(),
																									previous_dendritic_spike_times_global[j].end(),
																									internal_time - STDP_WINDOW);
																					
								for (auto it = it_relevant_spikes; it != previous_dendritic_spike_times_global[j].end(); it++)
								{
									double dt = *it - internal_time;
                            
									//std::cout << "dt in LTD = " << dt << std::endl;
									double w_before = weights_RA_RA_local[i][j];

									LTD(weights_RA_RA_local[i][j], dt);
									
									// show LTD of pool neurons only
									auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
									
									if (p.first == p.second)
									{
										std::cout << "LTD from  " << Id_RA_local[i] << " -> " << j
											  << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][j] - w_before
											  << std::endl;
									}
									
									
									//~ std::cout << "LTD from  " << Id_RA_local[i] << " -> " << j
											  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][j] - w_before
											  //~ << std::endl;
									//printf("LTD from neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n", Id_RA_local[i], j,
									//	spikes_in_trial_soma_local[i].back(), spike_times_dend_global[j], dt, weights_local[i][j] - w);
									
									update_synapse(i, j);
											
									
								}
							}
                        }
                    }
                }

              
                

            } // end if get fired soma
            
            // if some neuron produced dendritic spike, store this neuron in array
            if ( HVCRA_local[i].get_fired_dend() )
            {
				// show spikes of pool neurons only
                auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
                
                if (p.first == p.second)
					std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
                
				
                spikes_in_trial_dend_local[i].push_back(internal_time); // dend spike time relative to trial onset
                
                RA_neurons_bursted_local.push_back(Id_RA_local[i]);
                
                some_RA_neuron_bursted_local = 1;
                
               
                //std::cout << "HVC-RA " << Id_RA_local[i] << " dend spike at " << internal_time << std::endl;
                
                
                //printf("My rank = %d; RA neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], internal_time);
               
                //printf("My rank = %d; Dend neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], spike_times_dend_local[i]);

               
            } // end if get_fired_dend

		} // end for i -> N_RA_local

		for (int i = 0; i < N_I_local; i++)
		{
            HVCI_local[i].DP8_step_no_target_update();
            
            //  if some I neuron spikes, change conductance update array
            if (HVCI_local[i].get_fired())
            {
                //printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
                some_I_neuron_fired_local = 1;
                spikes_in_trial_interneuron_local[i].push_back(internal_time);

                size_t num_RA_targets = syn_ID_I_RA_local[i].size();
                // loop over all targets of fired neurons
                for (size_t j = 0; j < num_RA_targets; j++)
                {
                    int syn_ID = syn_ID_I_RA_local[i][j];
                    update_Gi_RA_local[syn_ID] += weights_I_RA_local[i][j];
                    //printf("Rank = %d; i = %d; update_Gi_RA_local[%d] = %f; weights_I_RA_local[%d][%d] = %f\n", MPI_rank, i, syn_ID,
                     //   update_Gi_RA_local[syn_ID], weights_I_RA_local[fired_ID][j], fired_ID, j);
                }
                
                //std::cout << "HVC-I " << Id_I_local[i] << " spike at " << internal_time << std::endl;
                
            }
		} // end if i -> N_I_local

        // if we need to update network state or if we reached the end of the trial
        // get if any neurons fired in some process

        if ( ( internal_time > network_time ) || ( t == num_steps-1 ) )
        {
            //if (MPI_rank == 0)
            //{
            //    std::cout << "internal time = " << internal_time << " > network_time = " << network_time << std::endl; 
            //}

            MPI_Allreduce(&some_RA_neuron_spiked_local, &some_RA_neuron_spiked_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_RA_neuron_bursted_local, &some_RA_neuron_bursted_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&some_I_neuron_fired_local, &some_I_neuron_fired_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
            if (some_I_neuron_fired_global > 0)
            {
            // sum update array and send to all processes

                MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                for (int i = 0; i < N_RA_local; i++)
                    HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
            	
				// update conductance arrays and fired indicators
				some_I_neuron_fired_local = 0;
            	some_I_neuron_fired_global = 0;
				
				std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            }


            //if (some_RA_neuron_fired_global == 1)
            //    printf("Rank %d; some_RA_neuron_fired_global: %d\n", MPI_rank, some_RA_neuron_fired_global);

            // if somatic compartment of any neuron in the pool fired, update synaptic conductances
             if ( some_RA_neuron_spiked_global > 0 )
             {
            // sum all update arrays and send to all processes

                MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                // now update excitatory conductances of all neurons
                for (int i = 0; i < N_RA_local; i++)
                    HVCRA_local[i].raiseExcWeight(update_Ge_RA_global[Id_RA_local[i]]); // update conductance
				
                for (int i = 0; i < N_I_local; i++)
                    HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);
                
				// update conductance arrays and fired indicators
            	some_RA_neuron_spiked_local = 0;
	        	some_RA_neuron_spiked_global = 0;
				
				std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
				
				std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
            
				//~ // gather all spiked HVC-RA neurons
				//~ this->gather_spiked_or_bursted_neurons(RA_neurons_spiked_local, RA_neurons_spiked_global);
				
				//~ if (MPI_rank == 0)
				//~ {
					//~ std::cout << "Spiked neurons:\n";
					//~ 
					//~ for (size_t j = 0; j < RA_neurons_spiked_global.size(); j++)
						//~ std::cout << RA_neurons_spiked_global[j] << " ";
					//~ std::cout << std::endl;
				//~ }
				
				//~ // add somatic spikes for spiked neurons
				//~ for (size_t j = 0; j < RA_neurons_spiked_global.size(); j++)
					//~ previous_somatic_spike_times_global[RA_neurons_spiked_global[j]].push_back(internal_time);
                  
              	/* 
                if (MPI_rank == 0)
                {
					for (size_t i = 0; i < RA_neurons_fired_dend_global.size(); i++)
                        printf("neuron %d bursted; spike_time_dend = %f\n", RA_neurons_fired_dend_global[i], 
								spikes_in_trial_dend_global[RA_neurons_fired_dend_global[i]].back());
                }
                */
                
                // clear arrays
				//~ RA_neurons_spiked_local.clear();
				//~ RA_neurons_spiked_global.clear();
			
			} // end if some HVC-RA neuron spiked
			
			if ( some_RA_neuron_bursted_global > 0 )
            {
				// update conductance arrays and fired indicators
            	some_RA_neuron_bursted_local = 0;
	        	some_RA_neuron_bursted_global = 0;
	        	
	        	// gather all bursted HVC-RA neurons
				this->gather_spiked_or_bursted_neurons(RA_neurons_bursted_local, RA_neurons_bursted_global);
				
				// add dendritic spikes for bursted neurons
				for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
					previous_dendritic_spike_times_global[RA_neurons_bursted_global[j]].push_back(internal_time);
               
				
                // apply LTP rule to the previous somatic spikes of pool neurons and dendritic spike of bursted neurons

                for (int i = 0; i < N_RA_local; i++)
                {
					int presyn_ID = Id_RA_local[i]; // real id of presynaptic neuron
					
                    // if neuron is saturated apply LTP only if spiked neurons are among supersynapse targets
                    if ( static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss )
                    {
                        for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_bursted_global[j]; // id of postsynaptic neuron
                            
                            // do not allow self-to-self synapses to emerge
							if ( presyn_ID != postsyn_ID )
							{
								std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
											supersynapses_local[i].end(), postsyn_ID );

								if ( pos!=supersynapses_local[i].end() )
								{
									// find previous somatic spikes that are not too old
									if ( !previous_somatic_spike_times_local[i].empty() )
									{
										std::vector<double>::iterator it_relevant_spikes = std::lower_bound(previous_somatic_spike_times_local[i].begin(),
																											previous_somatic_spike_times_local[i].end(),
																											internal_time - STDP_WINDOW);
																							
										for (auto it = it_relevant_spikes; it != previous_somatic_spike_times_local[i].end(); it++)
										{
											double dt = internal_time - *it;
									
										
											//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
											
											if (dt <= synaptic_params.T_0)
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												LTD(weights_RA_RA_local[i][postsyn_ID], dt);
												
												//~ std::cout   << "LTD from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
											}
											else
											{
												double w_before = weights_RA_RA_local[i][postsyn_ID];
												
												
												LTP(weights_RA_RA_local[i][postsyn_ID], dt);
												
												//~ std::cout   << "LTP from saturated " << presyn_ID << " -> " << postsyn_ID
															//~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
															//~ << std::endl;
												
											}	
											//double w = weights_local[i][fired_ID];
											
											//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
											 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
											 //           dt, weights_local[i][fired_ID] - w);
											
											update_synapse(i, postsyn_ID);
											
										}
									}
								}
							}
                        }

                    }
                    // if not saturated apply LTP for all dendritic spikes
                    else
                    {
                        for (size_t j = 0; j < RA_neurons_bursted_global.size(); j++)
                        {
                            int postsyn_ID = RA_neurons_bursted_global[j];
                            // don't allow self-to-self connections
                            if ( postsyn_ID != presyn_ID )
                            {
                                // find previous somatic spikes that are not too old
								if ( !previous_somatic_spike_times_local[i].empty() )
								{
									std::vector<double>::iterator it_relevant_spikes = std::lower_bound(previous_somatic_spike_times_local[i].begin(),
																										previous_somatic_spike_times_local[i].end(),
																										internal_time - STDP_WINDOW);
																						
									for (auto it = it_relevant_spikes; it != previous_somatic_spike_times_local[i].end(); it++)
									{
										double dt = internal_time - *it;
								

										//std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
										
									    if (dt <= synaptic_params.T_0)
									    {
											double w_before = weights_RA_RA_local[i][postsyn_ID];
												
											LTD(weights_RA_RA_local[i][postsyn_ID], dt);
											
											// show LTD on pool neurons only
											auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											
											if (p.first == p.second)
											{
												std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  << std::endl;
											}
											
											//~ std::cout << "LTD from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
											
										}
										else
										{
											double w_before = weights_RA_RA_local[i][postsyn_ID];
											
											LTP(weights_RA_RA_local[i][postsyn_ID], dt);
											
											// show LTP on pool neurons only
											auto p = std::equal_range(training_neurons.begin(), training_neurons.end(), postsyn_ID);
											
											if (p.first == p.second)
											{
												std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  << std::endl;
											}
											
											
											//~ std::cout << "LTP from  " << presyn_ID << " -> " << postsyn_ID
													  //~ << " dt = " << dt << " w_before = " << w_before << " dw = " << weights_RA_RA_local[i][postsyn_ID] - w_before
													  //~ << std::endl;
										}	
										//double w = weights_local[i][fired_ID];
										
										//printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
										 //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
										 //           dt, weights_local[i][fired_ID] - w);
										
										update_synapse(i, postsyn_ID);
									
									}
									
									
								}
                                
                            }
                        }
                   }

                } // end for i -> N_RA_local

				       

               // check if we need axon remodeling

                for (int i = 0; i < N_RA_local; i++)
                {
                    if ( (static_cast<int>(supersynapses_local[i].size()) == synaptic_params.Nss) && (remodeled_local[i] == 0) )
                    {
					    this->axon_remodeling(i);

				    }
                }
            
				// clear bursts
				RA_neurons_bursted_local.clear();
				RA_neurons_bursted_global.clear();
				
            
            } // end if some HVC-RA neuron bursted

           
            network_time += NETWORK_UPDATE_FREQUENCY;
            
            //printf("network_time = %f\n", network_time);
        }


        //MPI_Barrier(MPI_COMM_WORLD);
    }
    this->potentiation_decay();
    //printf("After potentiation decay")
    this->update_all_synapses();
	
	// calculate new gaba reverse potential
	// update rate
	for (int i = 0; i < N_RA_local; i++)
	{
		num_spikes_in_recent_trials_local[i].push_front(static_cast<int>(spikes_in_trial_soma_local[i].size()));
        
        //firing_rate_short_local[i] = std::accumulate(num_spikes_in_recent_trials[i].begin(), num_spikes_in_recent_trials[i].begin() + RATE_WINDOW_SHORT, 0.0)
          //                          / static_cast<double>(RATE_WINDOW_SHORT);

		// calculate firing rate in large window:
        firing_rate_long_local[i] = std::accumulate(num_spikes_in_recent_trials_local[i].begin(), num_spikes_in_recent_trials_local[i].end(), 0.0)
                                                                                            / static_cast<double>(maturation_params.RATE_WINDOW_LONG);
                                                                                            
		//if (std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].end(), 0.0) >= 1)
		//{
		//	std::cout << "Recent bursts of neuron " << Id_RA_local[i] << " :" << std::endl;
			
		//	for (size_t j = 0; j < num_bursts_in_recent_trials[i].size(); j++)
		//		std::cout << num_bursts_in_recent_trials[i][j] << "\t";
			
		//	std::cout << std::endl;
		//}
	}
    // update GABA potential based on firing rates
	//this->update_Ei();
    
    // gather all neurons that are to be replaced
    //this->gather_neurons_2replace();
    
    // if some neurons are to be replaced, replace them
    //if (replace_real_id_global.size() > 0)
    //    this->replace_neurons();
    
    
    
    
}

//~ 
//~ void NetworkGrowthSimulator::trial_burst_stdp(int training)
//~ {
	//~ int some_RA_neuron_fired_soma_local;
	//~ int some_RA_neuron_fired_dend_local;
//~ 
	//~ int some_I_neuron_fired_local;
	//~ int some_RA_neuron_fired_soma_global;
	//~ int some_RA_neuron_fired_dend_global;
//~ 
	//~ int some_I_neuron_fired_global;
//~ 
    //~ std::vector<double> spike_times_fired_dend_local;
    //~ 
    //~ std::vector<int> RA_neurons_fired_dend_global;
    //~ std::vector<int> RA_neurons_fired_dend_realID;
//~ 
    //~ 
    //~ internal_time = trial_number * trial_duration;
    //~ network_time = internal_time + network_update_frequency;
//~ 
    //~ // advance internal time of each neuron
    //~ for (int i = 0; i < N_RA_local; i++)
        //~ num_trials_after_replacement_local[i] += 1;
//~ 
    //~ //printf("network time = %f\n", network_time);
//~ 
    //~ // if training trial
    //~ if (training > 0)
        //~ this->set_training_current();
//~ 
	//~ // initialize update arrays and fired indicators
	//~ std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
	//~ std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
	//~ std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
	//~ std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
	//~ std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
	//~ std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
//~ 
    //~ some_RA_neuron_fired_soma_local = 0;
	//~ some_RA_neuron_fired_soma_global = 0;
//~ 
	//~ some_RA_neuron_fired_dend_local = 0;
    //~ some_RA_neuron_fired_dend_global = 0;
//~ 
    //~ some_I_neuron_fired_local = 0;
    //~ some_I_neuron_fired_global = 0;
	//~ 
    //~ // evolve dynamics
    //~ for (int t = 1; t < size; t++)
	//~ {
		//~ internal_time += timeStep;
		//~ 
		//~ for (int i = 0; i < N_RA_local; i++)
		//~ {
            //~ // set GABA potential
            //~ HVCRA_local[i].set_Ei(gaba_potential_local[i]);
            //~ 
            //~ // Debraband step
            //~ HVCRA_local[i].Debraband_step_no_target_update();
            //~ 
            //~ // if some neuron produced somatic spike, do LTD for all previous dendritic spikes
            //~ if (HVCRA_local[i].get_fired_soma())
            //~ {
                //~ spikes_in_trial_soma_local[i].push_back(internal_time - trial_number * trial_duration);
               //~ 
				//~ some_RA_neuron_fired_soma_local = 1;
				//~ 
				//~ // loop over all inhibitory targets of fired neurons
				//~ size_t num_I_targets = syn_ID_RA_I_local[i].size();
				//~ for (size_t j = 0; j < num_I_targets; j++)
				//~ {
					//~ int syn_ID = syn_ID_RA_I_local[i][j];
					//~ update_Ge_I_local[syn_ID] += weights_RA_I_local[i][j];
//~ 
				//~ }
				//~ 
				//~ // loop over all excitatory targets
				//~ size_t num_RA_targets = active_synapses_local[i].size();
				//~ for (size_t j = 0; j < num_RA_targets; j++)
				//~ {
					//~ int syn_ID = active_synapses_local[i][j];
					//~ update_Ge_RA_local[syn_ID] += weights_local[i][syn_ID];
				//~ }
//~ 
			//~ 
                //~ //for (int j = 0; j < last_soma_spikes_local[i].size(); j++)
                //~ //{
                  //~ //  printf("My rank = %d; Soma spikes of neuron %d:  spike_time = %f\n", MPI_rank, Id_RA_local[i],
                    //~ //    last_soma_spikes_local[i][j]);
//~ 
                //~ //}
//~ 
                //~ //printf("My rank = %d; Soma neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], spike_times_soma_local[i]);
//~ 
            //~ } // end if get fired soma
            //~ 
            //~ // if some neuron produced dendritic spike, store this neuron in array
            //~ if (HVCRA_local[i].get_fired_dend())
            //~ {
                //~ some_RA_neuron_fired_dend_local = 1;
                //~ spikes_in_trial_dend_local[i].push_back(internal_time - trial_number * trial_duration); // dend spike time relative to trial onset
                //~ 
                //~ last_dend_spike_time_local[i] = internal_time; // last burst time of HVC(RA) neuron
                //~ //printf("My rank = %d; RA neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], internal_time);
                //~ RA_neurons_fired_dend_realID.push_back(Id_RA_local[i]);
                //~ spike_times_fired_dend_local.push_back(internal_time);
                //~ //printf("My rank = %d; Dend neuron %d fired; spike_time = %f\n", MPI_rank, Id_RA_local[i], spike_times_dend_local[i]);
//~ 
                //~ // if neuron is saturated apply LTD only to supersynapses
                //~ if (static_cast<int>(supersynapses_local[i].size()) == Nss)
                //~ {
                    //~ for (size_t k = 0; k < supersynapses_local[i].size(); k++)
                    //~ {
                        //~ int supersynapse_id = supersynapses_local[i][k];
                        //~ 
                        //~ if ( Id_RA_local[i] != supersynapse_id )
                        //~ {
                            //~ double dt = last_dend_spike_time_global[supersynapse_id] - internal_time;
//~ 
                            //~ //std::cout << "dt in saturated LTD = " << dt << std::endl;
//~ 
                            //~ if (fabs(dt) < STDP_WINDOW)
                            //~ {
                                //~ //double w = weights_local[i][supersynapse_id];     
                       //~ 
                                //~ LTD_burst(weights_local[i][supersynapse_id], dt);
                                //~ 
                                //~ //printf("LTD from saturated neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n",
                                //~ //            Id_RA_local[i], supersynapse_id, spikes_in_trial_soma_local[i].back(), spike_times_dend_global[supersynapse_id],
                                //~ //            dt, weights_local[i][supersynapse_id] - w);
                                //~ 
                                //~ update_synapse(i, supersynapse_id);
                                //~ 
                            //~ }
                        //~ }
                    //~ }
//~ 
					//~ // if some supersynapse desaturated, update all synapses
                	//~ if (static_cast<int>(supersynapses_local[i].size()) < Nss)
						//~ for (int j = 0; j < N_RA; j++)
							//~ this->update_synapse(i, j);
                //~ }
                //~ // if not saturated apply LTD rule with the last dendritic spike of all neurons and add glutamate to all neuron except the fired one
                //~ else
                //~ {
                    //~ for (int j = 0; j < N_RA; j++)
                    //~ {
                        //~ if ( Id_RA_local[i] != j )
                        //~ {
                            //~ double dt = last_dend_spike_time_global[j] - internal_time;
                            //~ 
                            //~ //std::cout << "dt in LTD = " << dt << std::endl;
//~ 
                            //~ if (fabs(dt) < STDP_WINDOW)
                            //~ {
//~ 
                                //~ //double w = weights_local[i][j];
//~ 
                                //~ LTD_burst(weights_local[i][j], dt);
                     //~ 
                                //~ //printf("LTD from neuron %d onto %d; somatic spike: %f; dendritic spike: %f; dt = %f; dw = %f\n", Id_RA_local[i], j,
                                //~ //	spikes_in_trial_soma_local[i].back(), spike_times_dend_global[j], dt, weights_local[i][j] - w);
                                //~ 
                                //~ update_synapse(i, j);
                                    //~ 
                            //~ }
                        //~ }
                    //~ }
                //~ }
//~ 
            //~ } // end if get_fired_dend
//~ 
		//~ } // end for i -> N_RA_local
//~ 
		//~ for (int i = 0; i < N_I_local; i++)
		//~ {
            //~ HVCI_local[i].DP8_step_no_target_update();
            //~ 
            //~ //  if some I neuron spikes, change conductance update array
            //~ if (HVCI_local[i].get_fired())
            //~ {
                //~ //printf("My rank = %d; I neuron %d fired; spike_time = %f\n", MPI_rank, Id_I_local[i], internal_time);
                //~ some_I_neuron_fired_local = 1;
                //~ spikes_in_trial_interneuron_local[i].push_back(internal_time - trial_number * trial_duration);
//~ 
                //~ size_t num_RA_targets = syn_ID_I_RA_local[i].size();
                //~ // loop over all targets of fired neurons
                //~ for (size_t j = 0; j < num_RA_targets; j++)
                //~ {
                    //~ int syn_ID = syn_ID_I_RA_local[i][j];
                    //~ update_Gi_RA_local[syn_ID] += weights_I_RA_local[i][j];
                    //~ //printf("Rank = %d; i = %d; update_Gi_RA_local[%d] = %f; weights_I_RA_local[%d][%d] = %f\n", MPI_rank, i, syn_ID,
                     //~ //   update_Gi_RA_local[syn_ID], weights_I_RA_local[fired_ID][j], fired_ID, j);
                //~ }
            //~ }
		//~ } // end if i -> N_I_local
//~ 
        //~ // if we need to update network state or if we reached the end of the trial
        //~ // get if any neurons fired in some process
//~ 
        //~ if ( (internal_time > network_time) || (t == size-1) )
        //~ {
            //~ //if (MPI_rank == 0)
            //~ //{
            //~ //    std::cout << "internal time = " << internal_time << " > network_time = " << network_time << std::endl; 
            //~ //}
//~ 
            //~ MPI_Allreduce(&some_RA_neuron_fired_soma_local, &some_RA_neuron_fired_soma_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            //~ MPI_Allreduce(&some_RA_neuron_fired_dend_local, &some_RA_neuron_fired_dend_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            //~ MPI_Allreduce(&some_I_neuron_fired_local, &some_I_neuron_fired_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        //~ 
            //~ if (some_I_neuron_fired_global > 0)
            //~ {
            //~ // sum update array and send to all processes
//~ 
                //~ MPI_Allreduce(&update_Gi_RA_local[0], &update_Gi_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//~ 
                //~ for (int i = 0; i < N_RA_local; i++)
                //~ {
                    //~ HVCRA_local[i].raiseI(update_Gi_RA_global[Id_RA_local[i]]);
                //~ }
            	//~ 
				//~ // update conductance arrays and fired indicators
				//~ some_I_neuron_fired_local = 0;
            	//~ some_I_neuron_fired_global = 0;
				//~ 
				//~ std::fill(update_Gi_RA_local.begin(), update_Gi_RA_local.end(), 0.0);
				//~ std::fill(update_Gi_RA_global.begin(), update_Gi_RA_global.end(), 0.0);
            //~ }
//~ 
//~ 
            //~ //if (some_RA_neuron_fired_global == 1)
            //~ //    printf("Rank %d; some_RA_neuron_fired_global: %d\n", MPI_rank, some_RA_neuron_fired_global);
//~ 
            //~ // if somatic compartment of any neuron in the pool fired, update synaptic conductances
             //~ if (some_RA_neuron_fired_soma_global > 0)
             //~ {
            //~ // sum all update arrays and send to all processes
//~ 
                //~ MPI_Allreduce(&update_Ge_RA_local[0], &update_Ge_RA_global[0], N_RA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                //~ MPI_Allreduce(&update_Ge_I_local[0], &update_Ge_I_global[0], N_I, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//~ 
                //~ // now update excitatory conductances of all neurons
                //~ for (int i = 0; i < N_RA_local; i++)
                //~ {
                    //~ HVCRA_local[i].raiseE(update_Ge_RA_global[Id_RA_local[i]]); // update conductance
				//~ }
//~ 
                //~ for (int i = 0; i < N_I_local; i++)
                //~ {
                    //~ HVCI_local[i].raiseE(update_Ge_I_global[Id_I_local[i]]);
                //~ }
				//~ 
				//~ // update conductance arrays and fired indicators
            	//~ some_RA_neuron_fired_soma_local = 0;
	        	//~ some_RA_neuron_fired_soma_global = 0;
				//~ 
				//~ std::fill(update_Ge_RA_local.begin(), update_Ge_RA_local.end(), 0.0);
				//~ std::fill(update_Ge_RA_global.begin(), update_Ge_RA_global.end(), 0.0);
				//~ 
				//~ std::fill(update_Ge_I_local.begin(), update_Ge_I_local.end(), 0.0);
				//~ std::fill(update_Ge_I_global.begin(), update_Ge_I_global.end(), 0.0);
            //~ }
//~ 
            //~ // if dendritic compartment of any neuron in the pool fired, update weights
            //~ if (some_RA_neuron_fired_dend_global > 0)
            //~ {
                //~ // gather all bursted neurons
                //~ this->gather_bursts(RA_neurons_fired_dend_realID, RA_neurons_fired_dend_global, spike_times_fired_dend_local);
//~ 
              	//~ /* 
                //~ if (MPI_rank == 0)
                //~ {
					//~ for (size_t i = 0; i < RA_neurons_fired_dend_global.size(); i++)
                        //~ printf("neuron %d bursted; spike_time_dend = %f\n", RA_neurons_fired_dend_global[i], 
								//~ spikes_in_trial_dend_global[RA_neurons_fired_dend_global[i]].back());
                //~ }
                //~ */
//~ 
                //~ // apply LTP rule to the last dendritic bursts and dendritic spike of fired neurons
//~ 
                //~ for (int i = 0; i < N_RA_local; i++)
                //~ {
                    //~ // if neuron is saturated apply LTP only if dendritic spike occured in supersynapse
                    //~ if ( static_cast<int>(supersynapses_local[i].size()) == Nss)
                    //~ {
                        //~ for (size_t j = 0; j < RA_neurons_fired_dend_global.size(); j++)
                        //~ {
                            //~ int fired_ID = RA_neurons_fired_dend_global[j];
//~ 
                            //~ std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
                                        //~ supersynapses_local[i].end(), fired_ID);
//~ 
                            //~ if (pos!=supersynapses_local[i].end())
                            //~ {
						   //~ 
								//~ double dt = last_dend_spike_time_global[fired_ID] - last_dend_spike_time_local[i];
//~ 
								//~ //std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
								//~ if (dt < STDP_WINDOW)
								//~ {
									//~ if (dt <= T_0)
										//~ LTD_burst(weights_local[i][fired_ID], dt);
									//~ else
										//~ LTP_burst(weights_local[i][fired_ID], dt);
										//~ 
									//~ //double w = weights_local[i][fired_ID];
									//~ 
									//~ //printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
									 //~ //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
									 //~ //           dt, weights_local[i][fired_ID] - w);
									//~ 
									//~ update_synapse(i, fired_ID);
								//~ }
							//~ 
                            //~ }
                        //~ }
//~ 
                    //~ }
                    //~ // if not saturated apply LTP for all dendritic spikes
                    //~ else
                    //~ {
                        //~ for (size_t j = 0; j < RA_neurons_fired_dend_global.size(); j++)
                        //~ {
                            //~ int fired_ID = RA_neurons_fired_dend_global[j];
                            //~ // don't allow self-to-self connections
                            //~ if (fired_ID != Id_RA_local[i])
                            //~ {
                                //~ // apply to the last dendritic burst
                                //~ 
								//~ double dt = last_dend_spike_time_global[fired_ID] - last_dend_spike_time_local[i];
//~ 
								//~ //std::cout << "From neuron " << Id_RA_local[i] << " to neuron " << fired_ID << " dt = " << dt << std::endl;
								//~ if (dt < STDP_WINDOW)
								//~ {
								   //~ if (dt <= T_0)
										//~ LTD_burst(weights_local[i][fired_ID], dt);
									//~ else
										//~ LTP_burst(weights_local[i][fired_ID], dt);
										//~ 
									//~ //double w = weights_local[i][fired_ID];
									//~ 
									//~ //printf("LTP from saturated neuron %d onto %d; somatic spike at %f; dendritic spike at %f; dt = %f; dw = %f\n", 
									 //~ //           Id_RA_local[i], fired_ID, spikes_in_trial_soma_local[i][k], spike_times_dend_global[fired_ID], 
									 //~ //           dt, weights_local[i][fired_ID] - w);
									//~ 
									//~ update_synapse(i, fired_ID);
									//~ 
								//~ }
                                //~ 
                            //~ }
                        //~ }
                   //~ }
//~ 
                //~ } // end for i -> N_RA_local
//~ 
               //~ // check if we need axon remodeling
//~ 
                //~ for (int i = 0; i < N_RA_local; i++)
                //~ {
                    //~ if ( (static_cast<int>(supersynapses_local[i].size()) == Nss) && (remodeled_local[i] == 0) )
                    //~ {
					    //~ this->axon_remodeling(i);
//~ 
				    //~ }
                //~ }
	        	//~ 
				//~ // update fired arrays and indicators
				//~ some_RA_neuron_fired_dend_local = 0;
            	//~ some_RA_neuron_fired_dend_global = 0;
            	//~ 
				//~ RA_neurons_fired_dend_global.clear();
            	//~ RA_neurons_fired_dend_realID.clear();
//~ 
                //~ spike_times_fired_dend_local.clear();
            //~ } // end if some_neuron_dend_fired
            //~ 
            //~ network_time += network_update_frequency;
            //~ 
            //~ //printf("network_time = %f\n", network_time);
//~ 
        //~ }
//~ 
//~ 
        //~ //MPI_Barrier(MPI_COMM_WORLD);
    //~ }
    //~ this->potentiation_decay();
    //~ //printf("After potentiation decay")
    //~ this->update_all_synapses();
	//~ 
	//~ // calculate new gaba reverse potential
	//~ // update rate
	//~ for (int i = 0; i < N_RA_local; i++)
	//~ {
		//~ num_bursts_in_recent_trials[i].push_front(static_cast<int>(spikes_in_trial_dend_local[i].size()));
        //~ 
        //~ firing_rate_short_local[i] = std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].begin() + RATE_WINDOW_SHORT, 0.0)
                                    //~ / static_cast<double>(RATE_WINDOW_SHORT);
//~ 
		//~ // calculate firing rate in large window:
        //~ firing_rate_long_local[i] = std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].end(), 0.0)
                                                                                            //~ / static_cast<double>(RATE_WINDOW_LONG);
                                                                                            //~ 
		//~ //if (std::accumulate(num_bursts_in_recent_trials[i].begin(), num_bursts_in_recent_trials[i].end(), 0.0) >= 1)
		//~ //{
		//~ //	std::cout << "Recent bursts of neuron " << Id_RA_local[i] << " :" << std::endl;
			//~ 
		//~ //	for (size_t j = 0; j < num_bursts_in_recent_trials[i].size(); j++)
		//~ //		std::cout << num_bursts_in_recent_trials[i][j] << "\t";
			//~ 
		//~ //	std::cout << std::endl;
		//~ //}
	//~ }
    //~ // update GABA potential based on firing rates
	//~ this->update_Ei();
    //~ 
    //~ // gather all neurons that are to be replaced
    //~ //this->gather_neurons_2replace();
    //~ 
    //~ // if some neurons are to be replaced, replace them
    //~ //if (replace_real_id_global.size() > 0)
    //~ //    this->replace_neurons();
//~ }

void NetworkGrowthSimulator::gather_localVectorsToGlobal(const std::vector<int>& vector_local, 
													std::vector<int>& vector_global)
{
        int size_local = static_cast<int>(vector_local.size());
        int size_global;
        
        int* recvcounts = new int[MPI_size];
        int* displs = new int[MPI_size];

        // get total number of elements
        MPI_Allreduce(&size_local, &size_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
        vector_global.resize(size_global);
        
        
        // get array with local vectors size in each process
        MPI_Allgather(&size_local, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, MPI_COMM_WORLD);

        displs[0] = 0;
        for (int i = 1; i < MPI_size; i++)
            displs[i] = displs[i-1] + recvcounts[i-1];


        // get local elements
        MPI_Allgatherv(&vector_local[0], size_local, MPI_INT,
            &vector_global[0], recvcounts, displs, MPI_INT, MPI_COMM_WORLD);
        

        delete [] recvcounts;
        delete [] displs;
}

void NetworkGrowthSimulator::gather_spiked_or_bursted_neurons(const std::vector<int>& RA_neurons_local, 
													std::vector<int>& RA_neurons_global)
{
        int num_RA_spiked_local = static_cast<int>(RA_neurons_local.size());
        int num_RA_spiked_global;
        //printf("Rank %d; num_RA_fired_local: %d\n", MPI_rank, num_RA_fired_local);

        int* recvcounts = new int[MPI_size];
        int* displs = new int[MPI_size];

        // get total number of fired neurons
        MPI_Allreduce(&num_RA_spiked_local, &num_RA_spiked_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
        RA_neurons_global.resize(num_RA_spiked_global);
        
        //printf("Rank %d; fired RA num: %d\n", MPI_rank, num_RA_fired_local);

        //if (MPI_rank == 0)
        //{
        //   printf("Master; fired RA num global: %d\n", num_RA_fired_global);
        //   printf("Master; RA_neurons_fired_global.size(): %d\n", RA_neurons_fired_global.size());
        //}

        // get array with number of fired neurons in each process
        MPI_Allgather(&num_RA_spiked_local, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, MPI_COMM_WORLD);

        displs[0] = 0;
        for (int i = 1; i < MPI_size; i++)
            displs[i] = displs[i-1] + recvcounts[i-1];


        //for (int i = 0; i < RA_neurons_fired_realID.size(); i++)
        //        printf("Rank %d; fired RA neuron: %d\n", MPI_rank, RA_neurons_fired_realID[i]);

        // get fired neurons
        MPI_Allgatherv(&RA_neurons_local[0], num_RA_spiked_local, MPI_INT,
            &RA_neurons_global[0], recvcounts, displs, MPI_INT, MPI_COMM_WORLD);
        
        /*
        if (MPI_rank == 0)
         {
            printf("Master; RA_neurons_fired_dend_global.size(): %d\n", RA_neurons_fired_dend_global.size());
            for (int i = 0; i < RA_neurons_fired_dend_global.size(); i++)
                printf("Rank %d; Dend fired RA neuron: %d; spike_time = %f\n", MPI_rank, RA_neurons_fired_dend_global[i],
                        internal_time);
        }
         */
        

        delete [] recvcounts;
        delete [] displs;
}

//~ 
//~ void NetworkGrowthSimulator::gather_bursts(std::vector<int>& RA_neurons_fired_dend_realID, std::vector<int>& RA_neurons_fired_dend_global, 
                                 //~ std::vector<double>& spike_times_fired_dend_local)
//~ {
        //~ std::vector<double> spike_times_fired_dend_global;
        //~ int num_RA_fired_dend_local = RA_neurons_fired_dend_realID.size();
        //~ int num_RA_fired_dend_global;
        //~ //printf("Rank %d; num_RA_fired_local: %d\n", MPI_rank, num_RA_fired_local);
//~ 
        //~ int* recvcounts = new int[MPI_size];
        //~ int* displs = new int[MPI_size];
//~ 
        //~ // get total number of fired neurons
        //~ MPI_Allreduce(&num_RA_fired_dend_local, &num_RA_fired_dend_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        //~ RA_neurons_fired_dend_global.resize(num_RA_fired_dend_global);
        //~ spike_times_fired_dend_global.resize(num_RA_fired_dend_global);
//~ 
        //~ //printf("Rank %d; fired RA num: %d\n", MPI_rank, num_RA_fired_local);
//~ 
        //~ //if (MPI_rank == 0)
        //~ //{
        //~ //   printf("Master; fired RA num global: %d\n", num_RA_fired_global);
        //~ //   printf("Master; RA_neurons_fired_global.size(): %d\n", RA_neurons_fired_global.size());
        //~ //}
//~ 
        //~ // get array with number of fired neurons in each process
        //~ MPI_Allgather(&num_RA_fired_dend_local, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, MPI_COMM_WORLD);
//~ 
        //~ displs[0] = 0;
        //~ for (int i = 1; i < MPI_size; i++)
        //~ {
            //~ displs[i] = displs[i-1] + recvcounts[i-1];
        //~ }
//~ 
//~ 
        //~ //for (int i = 0; i < RA_neurons_fired_realID.size(); i++)
        //~ //        printf("Rank %d; fired RA neuron: %d\n", MPI_rank, RA_neurons_fired_realID[i]);
//~ 
        //~ // get fired neurons
        //~ MPI_Allgatherv(&RA_neurons_fired_dend_realID[0], num_RA_fired_dend_local, MPI_INT,
            //~ &RA_neurons_fired_dend_global[0], recvcounts, displs, MPI_INT, MPI_COMM_WORLD);
        //~ 
        //~ MPI_Allgatherv(&spike_times_fired_dend_local[0], num_RA_fired_dend_local, MPI_DOUBLE,
            //~ &spike_times_fired_dend_global[0], recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        //~ /*
        //~ if (MPI_rank == 0)
         //~ {
            //~ printf("Master; RA_neurons_fired_dend_global.size(): %d\n", RA_neurons_fired_dend_global.size());
            //~ for (int i = 0; i < RA_neurons_fired_dend_global.size(); i++)
                //~ printf("Rank %d; Dend fired RA neuron: %d; spike_time = %f\n", MPI_rank, RA_neurons_fired_dend_global[i],
                        //~ internal_time);
        //~ }
         //~ */
        //~ // change spike times
        //~ for (size_t i = 0; i < RA_neurons_fired_dend_global.size(); i++)
        //~ {
			//~ //if (MPI_rank == 0)
			//~ //	std::cout << "neuron " << RA_neurons_fired_dend_global[i] << " bursted at " << spike_times_fired_dend_global[i] << std::endl;
            //~ 
			//~ spikes_in_trial_dend_global[RA_neurons_fired_dend_global[i]].push_back(spike_times_fired_dend_global[i]); 
			//~ last_dend_spike_time_global[RA_neurons_fired_dend_global[i]] = spike_times_fired_dend_global[i]; 
			  //~ 
        //~ }
//~ 
        //~ delete [] recvcounts;
        //~ delete [] displs;
//~ }
//~ 
//~ void NetworkGrowthSimulator::gather_neurons_2replace()
//~ {
	//~ // gather all neurons that are to be replaced
	//~ // get number of neurons to replace
	//~ int num_2replace_local = static_cast<int>(replace_local_id_local.size()); // number of neurons to replace on the process
//~ 
	//~ int num_2replace_total; // total number of neurons to replace
//~ 
	//~ MPI_Allreduce(&num_2replace_local, &num_2replace_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//~ 
	//~ replace_local_id_global.resize(num_2replace_total);
	//~ replace_real_id_global.resize(num_2replace_total);
	//~ replace_process_rank_global.resize(num_2replace_total);
//~ 
	//~ // gather ids and process ranks of neurons to replace
	//~ int* recvcounts = new int[MPI_size];
	//~ int* displs = new int[MPI_size];
//~ 
	//~ // get array with number of neurons to replace in each process
	//~ MPI_Allgather(&num_2replace_local, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, MPI_COMM_WORLD);
//~ 
	//~ displs[0] = 0;
	//~ for (int i = 1; i < MPI_size; i++)
	//~ {
	//~ displs[i] = displs[i-1] + recvcounts[i-1];
	//~ }
//~ 
	//~ // get neurons to replace
	//~ MPI_Allgatherv(&replace_local_id_local[0], num_2replace_local, MPI_INT,
	//~ &replace_local_id_global[0], recvcounts, displs, MPI_INT, MPI_COMM_WORLD);
//~ 
	//~ MPI_Allgatherv(&replace_real_id_local[0], num_2replace_local, MPI_INT,
	//~ &replace_real_id_global[0], recvcounts, displs, MPI_INT, MPI_COMM_WORLD);
//~ 
	//~ MPI_Allgatherv(&replace_process_rank_local[0], num_2replace_local, MPI_INT,
	//~ &replace_process_rank_global[0], recvcounts, displs, MPI_INT, MPI_COMM_WORLD);
//~ 
	//~ delete[] recvcounts;
	//~ delete[] displs;
//~ }
//~ 
//~ void NetworkGrowthSimulator::kill_neuron(int local_id, int global_id, int process_rank)
//~ {
    //~ // if process contains neuron to replace
    //~ if (MPI_rank == process_rank)
    //~ {
        //~ // zero all local weight matrices and delete all existing synapses from and out of the neuron
        //~ // zero weights from neuron
        //~ for (int i = 0; i < N_RA; i++)
            //~ std::fill(weights_local[local_id].begin(), weights_local[local_id].end(), 0.0);
//~ 
        //~ // delete active and super synapses
        //~ supersynapses_local[local_id].clear();
        //~ active_synapses_local[local_id].clear();
        //~ 
        //~ // set active and supersynapse indicators to zero
        //~ active_indicators_local[local_id].reset();
        //~ supersynapses_indicators_local[local_id].reset();
//~ 
        //~ // reset gaba_potential, remodeled and mature states
        //~ gaba_potential_local[local_id] = E_GABA_IMMATURE;
        //~ remodeled_local[local_id] = 0;
    //~ 
        //~ gaba_reached_mature_local[local_id] = 0;
	//~ 
        //~ firing_rate_short_local[local_id] = 0.0;
		//~ firing_rate_long_local[local_id] = 0.0;
    //~ 
        //~ num_trials_after_replacement_local[local_id] = 0;
        //~ 
        //~ last_dend_spike_time_local[local_id] = -200.0;
//~ 
        //~ // clear all recent bursts
        //~ for (int j = 0; j < RATE_WINDOW_LONG; j++)
            //~ num_bursts_in_recent_trials[local_id].push_back(0);
            //~ 
        //~ // set all variables to steady state
        //~ HVCRA_local[local_id].renew_neuron();
    //~ }
//~ 
    //~ // do for all processes
    //~ // erase all super and active synapses to a neuron
    //~ for (int i = 0; i < N_RA_local; i++)
    //~ {
        //~ std::vector<int>::iterator pos = std::find(active_synapses_local[i].begin(), active_synapses_local[i].end(), global_id);
//~ 
        //~ if (pos != active_synapses_local[i].end())
        //~ {
            //~ active_synapses_local[i].erase(pos);
            //~ 
            //~ // set active indicator to false
            //~ active_indicators_local[i][global_id] = false;
        //~ }
//~ 
        //~ pos = std::find(supersynapses_local[i].begin(), supersynapses_local[i].end(), global_id);
//~ 
        //~ if (pos != supersynapses_local[i].end())
        //~ {
            //~ // if neuron was saturated, set remodeled to false
            //~ if (static_cast<int>(supersynapses_local[i].size()) == Nss)
            //~ {
                //~ if (remodeled_local[i] != 1)
                    //~ std::cerr << "In replace_neuron. Neuron " << Id_RA_local[i] << " has maximum supersynapses but is not remodeled!" << std::endl;
                //~ 
            //~ }
            //~ else
            //~ {
                //~ if (remodeled_local[i] == 1)
                    //~ std::cerr << "In replace_neuron. Neuron " << Id_RA_local[i] << " is remodeled but its number of output supersynapses is "
                              //~ <<  supersynapses_local[i].size() << " which is smaller than maximum " << Nss << std::endl;
                //~ 
            //~ }
//~ 
            //~ supersynapses_local[i].erase(pos);
            //~ remodeled_local[i] = 0;   
            //~ 
            //~ // set supersynapses indicator to false
            //~ supersynapses_indicators_local[i][global_id] = false;
        //~ }
//~ 
    //~ }
//~ 
    //~ // zero weights to neuron
    //~ for (int i = 0; i < N_RA_local; i++)
        //~ weights_local[i][global_id] = 0.0;
//~ }
//~ 
//~ void NetworkGrowthSimulator::find_chain_neurons(std::vector<int>& chain_neurons)
//~ {
	//~ this->gather_data();
	//~ 
	//~ if (MPI_rank == 0)
	//~ {
		//~ // start checking with training neurons
		//~ std::vector<int> current_neurons_to_check = training_neurons; // current group of neurons to check outgoing connections for supersynapses
		//~ std::vector<int> next_neurons_to_check; // next group of neurons to check outgoing connections for supersynapses
	//~ 
		//~ do
		//~ {
			//~ for (size_t i = 0; i < current_neurons_to_check.size(); i++)
			//~ {
				//~ size_t num_supersynapses = supersynapses_global[current_neurons_to_check[i]].size();
				//~ 
				//~ for (size_t j = 0; j < num_supersynapses; j++)
				//~ {
					//~ // check if neuron is already in the chain
					//~ std::vector<int>::iterator pos = std::find(chain_neurons.begin(), chain_neurons.end(), supersynapses_global[current_neurons_to_check[i]][j]);
				//~ 
					//~ // if neuron is not already in the chain, add it to the chain and to neurons that will be checked next iteration 
					//~ if (pos == chain_neurons.end())
					//~ {
						//~ chain_neurons.push_back(supersynapses_global[current_neurons_to_check[i]][j]);
						//~ next_neurons_to_check.push_back(supersynapses_global[current_neurons_to_check[i]][j]);
					//~ }
				//~ }
				//~ 
				//~ 
			//~ }
			//~ current_neurons_to_check = next_neurons_to_check;
			//~ next_neurons_to_check.clear();
		//~ }
		//~ while (current_neurons_to_check.size() > 0);
	//~ }
	//~ 
//~ }
//~ 
//~ void NetworkGrowthSimulator::make_lesion(double fraction)
//~ {
	//~ // first estimate how many neurons are in the chain
	//~ std::vector<int> chain_neurons; // putative chain neurons
	//~ 
	//~ this->find_chain_neurons(chain_neurons);
	//~ 
	//~ if (MPI_rank == 0)
	//~ {
		//~ std::cout << "Number of neurons in the chain: " << chain_neurons.size() << std::endl 
				//~ << "Chain neurons: " << std::endl;
		//~ 
		//~ for (size_t i = 0; i < chain_neurons.size(); i++)
			//~ std::cout << chain_neurons[i] << "\t";
		//~ 
		//~ std::cout << std::endl;
	//~ }
	//~ 
	//~ int num_to_kill; // num neurons to kill
	//~ std::vector<int> neurons_to_kill; // neurons to kill
	//~ 
	//~ if (MPI_rank == 0)
	//~ {	
		//~ // estimate number of neurons to kill
		//~ num_to_kill = static_cast<int>(std::round(fraction * chain_neurons.size()));
	//~ 
		//~ // sample random indices of neurons to kill
		//~ std::vector<int>::iterator pos;
		//~ 
		//~ for (int i = 0; i < num_to_kill; i++)
		//~ {
			//~ int random_ind; // random index of neuron in chain_neurons array 
			//~ 
			//~ do
			//~ {
				//~ random_ind = generator.sample_integer(0, chain_neurons.size()-1);
				//~ 
				//~ // check if neuron is already selected
				//~ pos = std::find(neurons_to_kill.begin(), neurons_to_kill.end(), chain_neurons[random_ind]);
			//~ }
			//~ while (pos != neurons_to_kill.end());
			//~ 
			//~ neurons_to_kill.push_back(chain_neurons[random_ind]);
		//~ }
		//~ 
		//~ std::cout << "Neurons to kill: " << std::endl;
			//~ 
		//~ for (size_t i = 0; i < neurons_to_kill.size(); i++)
			//~ std::cout << neurons_to_kill[i] << "\t";
			//~ 
		//~ std::cout << std::endl;
	//~ }
	//~ 
	//~ // send number of neurons to kill to all processes
	//~ MPI_Bcast(&num_to_kill, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//~ 
	//~ neurons_to_kill.resize(num_to_kill);
	//~ 
	//~ // send neurons to kill to all processes
	//~ MPI_Bcast(&neurons_to_kill[0], num_to_kill, MPI_INT, 0, MPI_COMM_WORLD);
	//~ 
	//~ // calculate location of neurons to kill
	//~ int rank; // rank of the process that contains neuron to kill
	//~ int shift; // location within the process of the neuron to kill
	//~ 
	//~ for (int i = 0; i < num_to_kill; i++)
	//~ {
		//~ this->get_neuronRA_location(neurons_to_kill[i], &rank, &shift);
		//~ 
		//~ replace_local_id_global.push_back(shift);
		//~ replace_real_id_global.push_back(neurons_to_kill[i]);
		//~ replace_process_rank_global.push_back(rank);
		//~ 
	//~ }
	//~ 
	//~ if (MPI_rank == 1)
	//~ {
		//~ std::cout << "Rank 1; Neurons to kill: " << std::endl;
			//~ 
		//~ for (size_t i = 0; i < neurons_to_kill.size(); i++)
			//~ std::cout << "real id = " << replace_real_id_global[i] << " process = " << replace_process_rank_global[i]
					  //~ << " local id = " << replace_local_id_global[i] << std::endl;
		//~ 
			//~ 
		//~ std::cout << std::endl;
	//~ }
	//~ 
	//~ // replace neurons
	//~ this->replace_neurons();
//~ }
//~ 

//~ void NetworkGrowthSimulator::remove_neurons_from_network_local(std::vector<int>& neurons_to_remove)
//~ {
	//~ // remove input HVC-RA -> HVC-RA connections
	//~ for (int i = 0; i < N_RA_local; i++)
	//~ {
		//~ for (int j = 0; j < neurons_to_remove.size(); j++)
		//~ {
			//~ // remove input active connections
			//~ if ( active_indicators_local[i][neurons_to_remove[j]] == 1 )
			//~ {
				//~ auto it_active = std::find(active_synapses_local[i].begin(), active_synapses_local[i].end(), neurons_to_remove[j]);
			//~ 
				//~ if ( it_active != active_synapses_local[i].end() )
					//~ active_synapses_local[i].erase(it_active);
				//~ else
					//~ std::cerr << "Synapse with active indicator was not found in active synapse array" << std::endl;
				//~ 
				//~ active_indicators_local[i][neurons_to_remove[j]] = 0;
				//~ 
			//~ }
			//~ 
			//~ // remove input super connections
			//~ if ( supersynapses_indicators_local[i][neurons_to_remove[j]] == 1 )
			//~ {
				//~ std::cout << "Neuron " << neurons_to_remove[j] << " that is replaced had a supersynaptic connection "
						  //~ << Id_RA_local[i] << " -> " << neurons_to_remove[j] << std::endl;
						  //~ 
				//~ auto it_super = std::find(supersynapses_local[i].begin(), supersynapses_local[i].end(), neurons_to_remove[j]);
			//~ 
				//~ if ( it_super != supersynapses_local[i].end() )
					//~ supersynapses_local[i].erase(it_super);
				//~ else
					//~ std::cerr << "Synapse with super indicator was not found in super synapse array" << std::endl;
				//~ 
				//~ supersynapses_indicators_local[i][neurons_to_remove[j]] = 0;
				//~ 
			//~ }
			//~ 
			//~ // set weights to zero
			//~ weights_RA_RA_local[i][neurons_to_remove[j]] = 0.0;
		//~ }
	//~ }
	//~ 
	//~ // remove input HVC-I -> HVC-RA connections
	//~ for (int i = 0; i < N_I_local; i++)
	//~ {
		//~ for (int j = 0; j < neurons_to_remove.size(); j++)
		//~ {
			//~ auto it = std::find(syn_ID_I_RA_local[i].begin(), syn_ID_I_RA_local[i].end(), neurons_to_remove[j]);
			//~ 
			//~ if ( it != syn_ID_I_RA_local[i].end() )
			//~ {
				//~ int ind = std::distance(syn_ID_I_RA_local[i].begin(), it);
				//~ 
				//~ syn_ID_I_RA_local[i].erase(it);
				//~ weights_I_RA_local[i].erase(weights_I_RA_local[i].begin() + ind);
			//~ }
		//~ }
		//~ 
	//~ }
	//~ 
	//~ // remove output connections
	//~ for (int i = 0; i < neurons_to_remove.size(); i++)
	//~ {
		//~ // remove output connections to HVC-RA neurons
		//~ active_synapses_local[neurons_to_remove[i]].clear();
		//~ 
		//~ if ( !supersynapses_local[neurons_to_remove[i]].empty() )
			//~ std::cout << "Neuron " << neurons_to_remove[i] << " being replaced has supersynaptic outputs" << std::endl;
		//~ 
		//~ supersynapses_local[neurons_to_remove[i]].clear();
		//~ 
		//~ for (int j = 0; j < N_RA; j++)
		//~ {
			//~ supersynapses_indicators_local[neurons_to_remove[i]][j] = 0;
			//~ active_synapses_indicators_local[neurons_to_remove[i]][j] = 0;
			//~ 
			//~ weights_RA_RA_local[neurons_to_remove[i]][j] = 0.0;
		//~ }
		//~ 
		//~ // remove output connections to HVC-I neurons
		//~ syn_ID_RA_I_local[neurons_to_remove[i]].clear();
		//~ weights_RA_I_local[neurons_to_remove[i]].clear();
	//~ }
	//~ 
//~ }

void NetworkGrowthSimulator::remove_neurons_from_network(std::vector<int>& neurons_to_remove)
{
	// remove input HVC-RA -> HVC-RA connections
	for (int i = 0; i < N_RA_local; i++)
	{
		for (size_t j = 0; j < neurons_to_remove.size(); j++)
		{
			// remove input active connections
			if ( active_indicators_local[i][neurons_to_remove[j]] == 1 )
			{
				auto it_active = std::find(active_synapses_local[i].begin(), active_synapses_local[i].end(), neurons_to_remove[j]);
			
				if ( it_active != active_synapses_local[i].end() )
					active_synapses_local[i].erase(it_active);
				else
					std::cerr << "Synapse with active indicator was not found in active synapse array" << std::endl;
				
				active_indicators_local[i][neurons_to_remove[j]] = 0;
				
			}
			
			// remove input super connections
			if ( supersynapses_indicators_local[i][neurons_to_remove[j]] == 1 )
			{
				std::cout << "Neuron " << neurons_to_remove[j] << " that is replaced had an input supersynaptic connection "
						  << Id_RA_local[i] << " -> " << neurons_to_remove[j] << std::endl;
						  
				auto it_super = std::find(supersynapses_local[i].begin(), supersynapses_local[i].end(), neurons_to_remove[j]);
			
				if ( it_super != supersynapses_local[i].end() )
					supersynapses_local[i].erase(it_super);
				else
					std::cerr << "Synapse with super indicator was not found in super synapse array" << std::endl;
				
				supersynapses_indicators_local[i][neurons_to_remove[j]] = 0;
				remodeled_local[i] = 0;
			}
			
			// set weights to zero
			weights_RA_RA_local[i][neurons_to_remove[j]] = 0.0;
		}
	}
	
	// remove output connections
	for (size_t i = 0; i < neurons_to_remove.size(); i++)
	{
		// find location of HVC-RA neuron
		int rank, shift;
		
		this->get_neuronRA_location(neurons_to_remove[i], &rank, &shift);
		
		if (MPI_rank == rank)
		{	
			// remove output connections to HVC-RA neurons
			active_synapses_local[shift].clear();
			
			if ( !supersynapses_local[shift].empty() )
				std::cout << "Neuron " << neurons_to_remove[i] << " being replaced has supersynaptic outputs" << std::endl;
			
			supersynapses_local[shift].clear();
			
			for (int j = 0; j < N_RA; j++)
			{
				supersynapses_indicators_local[shift][j] = 0;
				active_indicators_local[shift][j] = 0;
				
				weights_RA_RA_local[shift][j] = 0.0;
			}
		}
	}
	
	if (MPI_rank == 0)
	{
		// remove input HVC-I -> HVC-RA connections
		for (int i = 0; i < N_I; i++)
		{
			for (size_t j = 0; j < neurons_to_remove.size(); j++)
			{
				auto it = std::find(syn_ID_I_RA_global[i].begin(), syn_ID_I_RA_global[i].end(), neurons_to_remove[j]);
				
				if ( it != syn_ID_I_RA_global[i].end() )
				{
					int ind = std::distance(syn_ID_I_RA_global[i].begin(), it);
					
					syn_ID_I_RA_global[i].erase(it);
					weights_I_RA_global[i].erase(weights_I_RA_global[i].begin() + ind);
					syn_lengths_I_RA_global[i].erase(syn_lengths_I_RA_global[i].begin() + ind);
					axonal_delays_I_RA_global[i].erase(axonal_delays_I_RA_global[i].begin() + ind);
				}
			}
			
		}
		
		for (size_t i = 0; i < neurons_to_remove.size(); i++)
		{
			// remove output connections to HVC-I neurons
			syn_ID_RA_I_global[neurons_to_remove[i]].clear();
			syn_lengths_RA_I_global[neurons_to_remove[i]].clear();
			axonal_delays_RA_I_global[neurons_to_remove[i]].clear();
			weights_RA_I_global[neurons_to_remove[i]].clear();
		}
	}
}

void NetworkGrowthSimulator::resample_neurons(std::vector<int>& neurons_to_replace)
{
	this->sample_coordinates_for_replaced(neurons_to_replace);
	this->sample_connectionsAndDelays_for_replaced(neurons_to_replace);
	
	
}

void NetworkGrowthSimulator::replace_neurons(std::vector<int>& neurons_to_replace)
{
	this->remove_neurons_from_network(neurons_to_replace);

	if (MPI_rank == 0)
		this->resample_neurons(neurons_to_replace);
		
	this->send_connections_RAandI();
	this->send_axonal_delays_RA2RA();
	
	this->update_replaced_neurons(neurons_to_replace);
	this->set_neuron_properties();
}

//~ void NetworkGrowthSimulator::replace_neurons()
//~ {
    //~ // write active, supersynapses and weights before replacing neurons
    //~ this->gather_data();
//~ 
    //~ this->write_weights((outputDirectory + "weights_before_replacement_trial_" + std::to_string(trial_number) + "_.bin").c_str());
    //~ this->write_active_synapses((outputDirectory + "RA_RA_active_connections_before_replacement_trial_" + std::to_string(trial_number) + "_.bin").c_str());
    //~ this->write_supersynapses((outputDirectory + "RA_RA_super_connections_before_replacement_trial_" + std::to_string(trial_number) + "_.bin").c_str());
    //~ this->write_maturation_info((outputDirectory + "mature_before_replacement_trial_" + std::to_string(trial_number) + "_.bin").c_str());
	//~ this->write_last_dend_spike_times((outputDirectory + "last_dendritic_spike_times_before_replacement_trial_" + std::to_string(trial_number) + "_.bin").c_str());
//~ 
    //~ // loop through all neurons that are to be replaced
    //~ for (size_t i = 0; i < replace_local_id_global.size(); i++)
    //~ {
        //~ // zero all weight matrices and delete all existing synapses from and out of the neurons
        //~ this->kill_neuron(replace_local_id_global[i], replace_real_id_global[i], replace_process_rank_global[i]);
        //~ 
    //~ }
//~ 
	//~ // initialize new coordinates and connections for replaced neurons
	//~ if (MPI_rank == 0)
	//~ {
		//~ std::string extension = "_after_replacement_trial_" + std::to_string(trial_number) + "_";
		//~ std::string fileReplaced = outputDirectory + "replaced_neurons.bin";
		//~ 
		//~ networkGen.replace_neurons(replace_real_id_global, extension, outputDirectory);
		//~ 
		//~ // get new network
		//~ networkGen.get_network(&syn_ID_RA_I_global, &weights_RA_I_global, &syn_ID_I_RA_global, &weights_I_RA_global, &training_neurons);
		//~ this->write_replaced_neurons(replace_real_id_global, fileReplaced.c_str());
	//~ }
//~ 
	//~ // distribute new network among all processes
    //~ this->send_connections();
//~ 
    //~ // clear all arrays with neurons to be replaced
    //~ replace_local_id_local.clear();
    //~ replace_real_id_local.clear();
    //~ replace_process_rank_local.clear();
    //~ replace_local_id_global.clear();
    //~ replace_real_id_global.clear();
    //~ replace_process_rank_global.clear();
    //~ 
    //~ // update all synapses to change their states if necessary
    //~ this->update_all_synapses();
    //~ 
    //~ // write active, supersynapses and weights after replacing neurons
    //~ 
    //~ this->gather_data();
    //~ 
    //~ this->write_weights((outputDirectory + "weights_after_replacement_trial_" + std::to_string(trial_number) + "_.bin").c_str());
    //~ this->write_active_synapses((outputDirectory + "RA_RA_active_connections_after_replacement_trial_" + std::to_string(trial_number) + "_.bin").c_str());
    //~ this->write_supersynapses((outputDirectory + "RA_RA_super_connections_after_replacement_trial_" + std::to_string(trial_number) + "_.bin").c_str());
    //~ this->write_maturation_info((outputDirectory + "mature_after_replacement_trial_" + std::to_string(trial_number) + "_.bin").c_str());
	//~ this->write_last_dend_spike_times((outputDirectory + "last_dendritic_spike_times_after_replacement_trial_" + std::to_string(trial_number) + "_.bin").c_str());
//~ 
//~ }

//~ void NetworkGrowthSimulator::add_new_neurons(int N)
//~ {
    //~ this->gather_data();
//~ 
    //~ this->write_weights((outputDirectory + "weights_before_neuron_addition.bin").c_str());
    //~ this->write_active_synapses((outputDirectory + "RA_RA_active_connections_before_neuron_addition.bin").c_str());
    //~ this->write_supersynapses((outputDirectory + "RA_RA_super_connections_before_neuron_addition.bin").c_str());
    //~ this->write_maturation_info((outputDirectory + "mature_before_neuron_addition.bin").c_str());
	//~ 
    //~ int new_remain = N; // new neurons remaining to be distributed
    //~ // distribute new neurons among processes
//~ 
    //~ std::vector<int> num_new_on_each_process; // number of new neurons added to each process
//~ 
    //~ num_new_on_each_process.resize(MPI_size); // array with number of RA neurons per process
//~ 
    //~ // first make number of HVC(RA) neurons equal on each process if it differs
    //~ // find process where number of neurons is smaller
    //~ int first_process_with_smaller_neuron_num = -1; // id of the process starting from which all processes have number of neurons less by one
//~ 
    //~ for (int i = 1; i < MPI_size; i++)
    //~ {
        //~ if (N_RA_sizes[i] < N_RA_sizes[i-1])
            //~ first_process_with_smaller_neuron_num = i;       
    //~ }
    //~ 
    //~ // if such process is found, add 1 neuron to all processes with id >= first_process_with_smaller_neuron_num
    //~ if (first_process_with_smaller_neuron_num > 0)
    //~ {
        //~ for (int i = first_process_with_smaller_neuron_num; i < MPI_size; i++)
        //~ {
            //~ N_RA_sizes[i]++;
            //~ num_new_on_each_process[i]++; 
        //~ }
        //~ new_remain -= MPI_size - first_process_with_smaller_neuron_num; // fix number of neurons to be distributed
    //~ }
    //~ 
    //~ // now number of neurons is already equal. distribute all remaining new neurons
	//~ for (int i = 0; i < MPI_size; i++)
	//~ {
		//~ num_new_on_each_process[i] += new_remain / MPI_size;
        //~ N_RA_sizes[i] += new_remain / MPI_size;
	//~ }
	//~ new_remain = new_remain % MPI_size;
	//~ 
    //~ int j = 0;
//~ 
	//~ // distribute RA neurons
	//~ while (new_remain > 0)
	//~ {
		//~ N_RA_sizes[j] += 1;
        //~ num_new_on_each_process[j] += 1;
		//~ new_remain -= 1;
		//~ j += 1;
//~ 
		//~ if (j >= MPI_size)
			//~ j -= MPI_size;
	//~ }
   //~ 
    //~ int N_RA_local_old = N_RA_local; // number of HVC(RA) neurons on each process before new neurons were added
    //~ int N_RA_old = N_RA; // total number of HVC(RA) neurons before new neurons were added
    //~ N_RA_local = N_RA_sizes[MPI_rank]; // new number of neurons on each process after new neurons were added
    //~ N_RA = N_RA + N; // new number of HVC(RA) neurons in the network after new neurons were added
//~ 
	//~ printf("My rank = %d; neurons added = %d; N_RA_local_old = %d; N_RA_local = %d; N_RA_old = %d; N_RA = %d;\n", MPI_rank, num_new_on_each_process[MPI_rank], N_RA_local_old, N_RA_local, N_RA_old, N_RA);
//~ 
    //~ // resize all arrays
    //~ this->resize_arrays_for_RA(N_RA_local, N_RA);
//~ 
    //~ 
    //~ // assign id to added neurons
    //~ Id_RA_local.resize(N_RA_local);
//~ 
    //~ for (int i = 0; i < num_new_on_each_process[MPI_rank]; i++)
    //~ {
        //~ int start_id_on_each_process = N_RA_old; // real neuron id which is assigned to the first new neuron on each process   
//~ 
        //~ for (int j = 0; j < MPI_rank; j++)
            //~ start_id_on_each_process += num_new_on_each_process[j];
//~ 
        //~ Id_RA_local[N_RA_local_old + i] = start_id_on_each_process + i;
    //~ }
//~ 
    //~ // resize global array containing all ids of HVC(RA) neurons
    //~ Id_RA_global.resize(N_RA);
//~ 
	//~ int* recvcounts_id_RA = new int[MPI_size];
	//~ int* displs_id_RA = new int[MPI_size];
//~ 
	//~ if (MPI_rank == 0)
	//~ {
		//~ // send connections to all prclusters/ocesses
//~ 
		//~ recvcounts_id_RA[0] = N_RA_sizes[0];
		//~ displs_id_RA[0] = 0;
//~ 
		//~ for (int i = 1; i < MPI_size; i++)
		//~ {
			//~ recvcounts_id_RA[i] = N_RA_sizes[i];
			//~ displs_id_RA[i] = displs_id_RA[i-1] + recvcounts_id_RA[i-1];
		//~ }
    //~ }
//~ 
	//~ // get global array with neuronal ids
	//~ MPI_Gatherv(&Id_RA_local[0], N_RA_local, MPI_INT, &Id_RA_global[0], recvcounts_id_RA, displs_id_RA, MPI_INT, 0, MPI_COMM_WORLD);
//~ 
	//~ delete[] recvcounts_id_RA;
	//~ delete[] displs_id_RA;
//~ 
    //~ this->write_global_index_array((outputDirectory + "global_index_array_after_neuron_addition.bin").c_str());
    //~ /*
    //~ std::cout << "My rank = " << MPI_rank << "\nId_RA_local = ";
    //~ for (size_t i = 0; i < Id_RA_local.size(); i++)
        //~ std::cout << Id_RA_local[i] << "\t";
    //~ std::cout << std::endl;
    //~ */
    //~ if (MPI_rank == 0)
    //~ {
        //~ std::cout << "\nId_RA_global = ";
        //~ for (size_t i = 0; i < Id_RA_global.size(); i++)
            //~ std::cout << Id_RA_global[i] << "\t";
        //~ std::cout << std::endl;
	//~ }
    //~ 
    //~ std::fill(num_bursts_in_recent_trials.begin() + N_RA_local_old, num_bursts_in_recent_trials.end(), intBuffer(RATE_WINDOW_LONG));
	//~ 
    //~ for (size_t i = N_RA_local_old; i < num_bursts_in_recent_trials.size(); i++)
		//~ for (int j = 0; j < RATE_WINDOW_LONG; j++)
			//~ num_bursts_in_recent_trials[i].push_back(0);
//~ 
    //~ // initialize coordinates for new neurons
    //~ this->initialize_coordinates_for_added_neurons(N_RA_old);
//~ 
    //~ // initialize_connections for new neurons
    //~ this->initialize_connections_for_added_neurons(N_RA_old);
//~ 
    //~ // initialize new neurons
//~ 
    //~ for (int i = N_RA_local_old; i < N_RA_local; i++)
	//~ {	
		//~ HVCRA_local[i].set_noise_generator(&generator);
		//~ HVCRA_local[i].set_white_noise(white_noise_mean_soma, white_noise_std_soma, white_noise_mean_dend, white_noise_std_dend);
		//~ HVCRA_local[i].set_dynamics(trial_duration, timeStep);
	//~ 
        //~ gaba_potential_local[i] = E_GABA_IMMATURE;      
    //~ }
//~ 
    //~ this->gather_data();
//~ 
    //~ this->write_weights((outputDirectory + "weights_after_neuron_addition.bin").c_str());
    //~ this->write_active_synapses((outputDirectory + "RA_RA_active_connections_after_neuron_addition.bin").c_str());
    //~ this->write_supersynapses((outputDirectory + "RA_RA_super_connections_after_neuron_addition.bin").c_str());
    //~ this->write_maturation_info((outputDirectory + "mature_after_neuron_addition.bin").c_str());
//~ }

void NetworkGrowthSimulator::resize_arrays_for_master_process()
{
	// resize connections
    
	weights_RA_RA_global.resize(N_RA);
	axonal_delays_RA_RA_global.resize(N_RA);
	syn_lengths_RA_RA_global.resize(N_RA);
	
	weights_RA_I_global.resize(N_RA);
	syn_ID_RA_I_global.resize(N_RA);
	syn_lengths_RA_I_global.resize(N_RA);
	axonal_delays_RA_I_global.resize(N_RA);
	
	
	weights_I_RA_global.resize(N_I);
	syn_ID_I_RA_global.resize(N_I);
	syn_lengths_I_RA_global.resize(N_I);
	axonal_delays_I_RA_global.resize(N_I);
	
	for (int i = 0; i < N_RA; i++)
	{
		weights_RA_RA_global[i].resize(N_RA);
		axonal_delays_RA_RA_global[i].resize(N_RA);
		syn_lengths_RA_RA_global[i].resize(N_RA);
		
	}
    
	// active and super synapses
	supersynapses_global.resize(N_RA);
	active_synapses_global.resize(N_RA);

    // spikes in trials
    //previous_somatic_spike_times_local.resize(n_local);
    
    spikes_in_trial_soma_global.resize(N_RA);
    spikes_in_trial_dend_global.resize(N_RA);
    spikes_in_trial_interneuron_global.resize(N_I);
    
    remodeled_global.resize(N_RA);
    mature_global.resize(N_RA);
	maturation_scale_global.resize(N_RA);
	
	std::fill(maturation_scale_global.begin(), maturation_scale_global.end(), MATURATION_SCALE_SPONTANEOUS);
	
	gaba_potential_global.resize(N_RA);
	rest_potential_global.resize(N_RA);
	GCa_global.resize(N_RA);
	Gk_global.resize(N_RA);
	GNa_global.resize(N_RA);
	Ad_global.resize(N_RA);
	Rc_global.resize(N_RA);
	
	std::fill(GCa_global.begin(), GCa_global.end(), maturation_params.GCA_IMMATURE);
	std::fill(rest_potential_global.begin(), rest_potential_global.end(), maturation_params.E_REST_IMMATURE);
	
	firing_rate_long_global.resize(N_RA);
	
	num_trials_after_replacement_global.resize(N_RA);
	num_spikes_in_recent_trials_global.resize(N_RA);
	
	for (size_t i = 0; i < num_spikes_in_recent_trials_global.size(); i++)
		num_spikes_in_recent_trials_global[i].resize(maturation_params.RATE_WINDOW_LONG);
}

void NetworkGrowthSimulator::initialize_global_inhibitory_conductances(double time_resolution_conductance)
{
	Ginh_global.resize(N_RA);
	
	for (int i = 0; i < N_RA; i++)
		Ginh_global[i].resize(int(TRIAL_DURATION / time_resolution_conductance));
		
	std::cout << "trial duration = " << TRIAL_DURATION << " "
			  << "time resolution for conductance = " << time_resolution_conductance << " "
			  << "size of inhibitory array = " << Ginh_global[0].size() << "\n";

}

void NetworkGrowthSimulator::initialize_local_inhibitory_conductances(double time_resolution_conductance)
{
	Ginh_local.resize(N_RA_local);
	
	for (int i = 0; i < N_RA_local; i++)
		Ginh_local[i].resize(int(TRIAL_DURATION / time_resolution_conductance));
}

void NetworkGrowthSimulator::resize_arrays_for_all_processes()
{
	HVCRA_local.resize(N_RA_local);
	HVCI_local.resize(N_I_local);

    Id_RA_local.resize(N_RA_local);
    Id_I_local.resize(N_I_local);

    // connections and their ids
    weights_RA_RA_local.resize(N_RA_local);
    axonal_delays_RA_RA_local.resize(N_RA_local);
    
    weights_RA_I_local.resize(N_RA_local);
    syn_ID_RA_I_local.resize(N_RA_local);
    axonal_delays_RA_I_local.resize(N_RA_local);

	weights_I_RA_local.resize(N_I_local);
	syn_ID_I_RA_local.resize(N_I_local);
    axonal_delays_I_RA_local.resize(N_I_local);

    // super and active indicators
    supersynapses_local.resize(N_RA_local);
	active_synapses_local.resize(N_RA_local);
    
    active_indicators_local.resize(N_RA_local);
    supersynapses_indicators_local.resize(N_RA_local);
	
	rescaled_indicators_local.resize(N_RA_local);
    
	// spikes in trials
    spikes_in_trial_soma_local.resize(N_RA_local);
    spikes_in_trial_dend_local.resize(N_RA_local);
    spikes_in_trial_interneuron_local.resize(N_I_local);
	
	//previous_somatic_spike_times_global.resize(N_RA);
	previous_dendritic_spike_times_global.resize(N_RA);
	
	previous_somatic_spike_times_local.resize(N_RA_local);
	
	// delivery queues
	delivery_queue_RA_RA_soma.resize(N_RA_local);
	delivery_queue_RA_I.resize(N_RA_local);
	delivery_queue_I_RA.resize(N_I_local);

	delivered_spike_times.resize(N_RA_local);
	
	mature_global.resize(N_RA);
	
	for (int i = 0; i < N_RA_local; i++)
		delivered_spike_times[i].resize(N_RA);

    // axon remodeling
	remodeled_local.resize(N_RA_local);
	
	// maturation parameters
	gaba_potential_local.resize(N_RA_local);
	rest_potential_local.resize(N_RA_local);
	
	GCa_local.resize(N_RA_local);
	Gk_local.resize(N_RA_local);
	GNa_local.resize(N_RA_local);
	Ad_local.resize(N_RA_local);
	Rc_local.resize(N_RA_local);
	
	//gaba_reached_mature_local.resize(n_local);
	
	mature_local.resize(N_RA_local);
	maturation_scale_local.resize(N_RA_local);
	
	std::fill(maturation_scale_local.begin(), maturation_scale_local.end(), MATURATION_SCALE_SPONTANEOUS);
	

	// keeping track of neuron replacement
	firing_rate_long_local.resize(N_RA_local);
	
	//firing_rate_short_local.resize(n_local);
	//firing_rate_short_global.resize(n_total);

    num_trials_after_replacement_local.resize(N_RA_local);
   
    // initialize circular buffers for somatic spikes
	num_spikes_in_recent_trials_local.resize(N_RA_local);
	
	std::fill(num_spikes_in_recent_trials_local.begin(), num_spikes_in_recent_trials_local.end(), intBuffer(maturation_params.RATE_WINDOW_LONG));

	for (size_t i = 0; i < num_spikes_in_recent_trials_local.size(); i++)
		for (int j = 0; j < maturation_params.RATE_WINDOW_LONG; j++)
			num_spikes_in_recent_trials_local[i].push_back(0);
	
		
    // initialize last dendritic spike times
    //std::fill(last_dend_spike_time_global.begin(), last_dend_spike_time_global.end(), -200.0);
    //std::fill(last_dend_spike_time_local.begin(), last_dend_spike_time_local.end(), -200.0);
    
    // initialize arrays
	for (int i = 0; i < N_RA_local; i++)
	{
        weights_RA_RA_local[i].resize(N_RA);
        axonal_delays_RA_RA_local[i].resize(N_RA);
        active_indicators_local[i].resize(N_RA);
        supersynapses_indicators_local[i].resize(N_RA);
        rescaled_indicators_local[i].resize(N_RA);
	}
}
//~ 
//~ void NetworkGrowthSimulator::update_Ei()
//~ {
	//~ for (int i = 0; i < N_RA_local; i++)
	//~ {
        //~ // check if neuron is among training neurons
        //~ std::vector<int>::iterator pos = std::find(training_neurons.begin(), training_neurons.end(), Id_RA_local[i]);
//~ 
		//~ 
	//~ // if not a training neuron, GABA didnt hit the bottom and firing rate exceeds gaba rate threshold, decrease gaba
		//~ if (pos == training_neurons.end())
		//~ {
		//~ 
			//~ if ( (!gaba_reached_mature_local[i]) && (firing_rate_short_local[i] >= GABA_RATE_THRESHOLD))
			//~ {
				//~ gaba_potential_local[i] = gaba_potential_local[i] - GABA_DOWN * static_cast<int>(spikes_in_trial_dend_local[i].size()); // decrease gaba potential
						//~ 
				//~ if (gaba_potential_local[i] <= E_GABA_MATURE) // if gaba potential hit the bottom, make neuron mature
				//~ {
					//~ gaba_potential_local[i] = E_GABA_MATURE;
					//~ gaba_reached_mature_local[i] = true;
				//~ }
//~ 
			//~ }
//~ 
			//~ // if not a training neuron and if firing rate in wide window is smaller than death threshold, replace neuron with a new one
			//~ else if ( (firing_rate_long_local[i] <= DEATH_RATE_THRESHOLD) && (num_trials_after_replacement_local[i] >= RATE_WINDOW_LONG) )
			//~ {
				//~ replace_local_id_local.push_back(i);
				//~ replace_real_id_local.push_back(Id_RA_local[i]);
				//~ replace_process_rank_local.push_back(MPI_rank);
			//~ }
		//~ }
	//~ }
//~ }

void NetworkGrowthSimulator::potentiation_decay_sudden_maturation()
{
    for (int i = 0; i < N_RA_local; i++)
    {
        for (int j = 0; j < N_RA; j++)
        {
			if ( rescaled_indicators_local[i][j] != 1 )
			{
				if (weights_RA_RA_local[i][j] >= synaptic_params.SUPERSYNAPSE_THRESHOLD)
					weights_RA_RA_local[i][j] -= synaptic_params.BETA_SUPERSYNAPSE;
				else
					weights_RA_RA_local[i][j] -= synaptic_params.BETA;
					
				
			}
			else if (active_indicators_local[i][j] == 0)
				weights_RA_RA_local[i][j] -= synaptic_params.BETA;
				
			if (weights_RA_RA_local[i][j] < 0)
					weights_RA_RA_local[i][j] = 0;
        }
    }
}


void NetworkGrowthSimulator::potentiation_decay()
{
    for (int i = 0; i < N_RA_local; i++)
    {
        for (int j = 0; j < N_RA; j++)
        {
            if (weights_RA_RA_local[i][j] >= synaptic_params.SUPERSYNAPSE_THRESHOLD)
                weights_RA_RA_local[i][j] -= synaptic_params.BETA_SUPERSYNAPSE;
            else
                weights_RA_RA_local[i][j] -= synaptic_params.BETA;
                
            if (weights_RA_RA_local[i][j] < 0)
				weights_RA_RA_local[i][j] = 0;
        }
    }
}

void NetworkGrowthSimulator::axon_remodeling(int i)
{
	// erase all active synapses if they are not among supersynapses
	
	// first erase all active synapses
	for (size_t j = 0; j < active_synapses_local[i].size(); j++)
	{
		int target_ID = active_synapses_local[i][j];
		active_indicators_local[i][target_ID] = 0;
	}
	
	// now add all supersynapses as active synapses
	for ( size_t j = 0; j < supersynapses_local[i].size(); j++ )
	{
		int target_ID = supersynapses_local[i][j];
		active_indicators_local[i][target_ID] = 1;
	}
	
	active_synapses_local[i] = supersynapses_local[i];
	
	
	//~ std::vector<int> active_to_erase; // active synapses to erase
//~ 
	//~ // find all active synapses needed to be erased
    //~ for (size_t j = 0; j < active_synapses_local[i].size(); j++)
    //~ {
        //~ int syn_ID = active_synapses_local[i][j];
//~ 
        //~ // check if active synapse is among supersynapses
        //~ std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
                                //~ supersynapses_local[i].end(), syn_ID);
        //~ // if not erase it
        //~ if (pos == supersynapses_local[i].end())
        //~ {
            //~ active_to_erase.push_back(active_synapses_local[i][j]);
            //~ active_indicators_local[i][syn_ID] = 0;
        //~ }
    //~ }
//~ 
	//~ // erase synapses
	//~ for (size_t j = 0; j < active_to_erase.size(); j++)
	//~ {
		//~ std::vector<int>::iterator pos = std::find(active_synapses_local[i].begin(), active_synapses_local[i].end(), active_to_erase[j]);
//~ 
		//~ if (pos != active_synapses_local[i].end())
			//~ active_synapses_local[i].erase(pos);
		//~ else
			//~ std::cout << "Active synapse " << Id_RA_local[i] << " -> " << active_to_erase[j]
					  //~ << " which should be removed in axon remodeling is not found!" << std::endl;
	//~ }

    remodeled_local[i] = 1;
}

void NetworkGrowthSimulator::update_all_synapses_sudden_maturation()
{
    //printf("Process %d; Updating all synapses after potentiation decay\n", MPI_rank);
    for (int i = 0; i < N_RA_local; i++)
    {
        for (int j = 0; j < N_RA; j++)
        {
			if ( mature_global[j] != 1 )
				this->update_synapse(i, j);
			
        }
    }

}

void NetworkGrowthSimulator::update_all_synapses()
{
    //printf("Process %d; Updating all synapses after potentiation decay\n", MPI_rank);
    for (int i = 0; i < N_RA_local; i++)
    {
        for (int j = 0; j < N_RA; j++)
        {
            this->update_synapse(i, j);
			
        }
    }
}

void NetworkGrowthSimulator::rescale_inhibition_to_mature(int neuron_id)
{
	// rescale inhibitory synapses to mature neuron
	for (int i = 0; i < N_I_local; i++)
	{
		auto it = std::find(syn_ID_I_RA_local[i].begin(), syn_ID_I_RA_local[i].end(), neuron_id);
		
		if ( it != syn_ID_I_RA_local[i].end() )
		{
			int ind = std::distance(syn_ID_I_RA_local[i].begin(), it);
			weights_I_RA_local[i][ind] /= MATURE_SYNAPSE_SCALING;
		}
	}
}

void NetworkGrowthSimulator::rescale_synapses_to_mature(int neuron_id)
{
	
	for (int i = 0; i < N_RA_local; i++)
	{
		if ( active_indicators_local[i][neuron_id] == 1 ){
			rescaled_indicators_local[i][neuron_id] = 1;
			weights_RA_RA_local[i][neuron_id] *= MATURE_SYNAPSE_SCALING;
		}
	}
	
	//~ // rescale inhibitory synapses to mature neuron
	//~ for (int i = 0; i < N_I_local; i++)
	//~ {
		//~ auto it = std::find(syn_ID_I_RA_local[i].begin(), syn_ID_I_RA_local[i].end(), neuron_id);
		//~ 
		//~ if ( it != syn_ID_I_RA_local[i].end() )
		//~ {
			//~ int ind = std::distance(syn_ID_I_RA_local[i].begin(), it);
			//~ weights_I_RA_local[i][ind] /= MATURE_SYNAPSE_SCALING;
		//~ }
	//~ }
	//~ 
	//~ // rescale synapses in global array
	//~ if ( MPI_rank == 0 ){
		//~ for (int i = 0; i < N_I; i++){
			//~ auto it = std::find(syn_ID_I_RA_global[i].begin(), syn_ID_I_RA_global[i].end(), neuron_id);
		//~ 
			//~ if ( it != syn_ID_I_RA_global[i].end() ){
				//~ int ind = std::distance(syn_ID_I_RA_global[i].begin(), it);
				//~ weights_I_RA_global[i][ind] /= MATURE_SYNAPSE_SCALING;
			//~ }
		//~ }
	//~ }
}

void NetworkGrowthSimulator::update_synapse_sudden_maturation(int i, int j)
{
	double w = weights_RA_RA_local[i][j];

    if ( ( w >= synaptic_params.ACTIVATION_THRESHOLD ) && ( active_indicators_local[i][j] == 0 ) && ( remodeled_local[i] == 0) 
						&& ( rescaled_indicators_local[i][j] == 0) )
    {
       // printf("Activated synapse from %d onto %d; weight = %f\n", Id_RA_local[i], j, w);
        active_indicators_local[i][j] = 1;
        active_synapses_local[i].push_back(j);
    }

    if ( (w >= synaptic_params.SUPERSYNAPSE_THRESHOLD) && (supersynapses_indicators_local[i][j] == 0) && (static_cast<int>(supersynapses_local[i].size()) < synaptic_params.Nss) 
									&& ( rescaled_indicators_local[i][j] == 0) )
    {
       // printf("Activated supersynapse from %d onto %d; weight = %f\n", Id_RA_local[i], j, w);
        supersynapses_indicators_local[i][j] = 1;
        supersynapses_local[i].push_back(j);
    }

    if ( (w < synaptic_params.SUPERSYNAPSE_THRESHOLD) && (supersynapses_indicators_local[i][j] == 1) )
    {
       // printf("Deactivated supersynapse from %d onto %d; weight = %f\n", Id_RA_local[i], j, w);
        supersynapses_indicators_local[i][j] = 0;
        remodeled_local[i] = 0;
        std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
                                                    supersynapses_local[i].end(), j);

        if (pos!= supersynapses_local[i].end())
            supersynapses_local[i].erase(pos);
        else
            std::cout << "Supersynapse " << Id_RA_local[i] << " -> " << j << " to be erased is not found!" << std::endl;
    }

    if ( (w < synaptic_params.ACTIVATION_THRESHOLD) && (active_indicators_local[i][j] == 1) )
    {
      //  printf("Deactivated synapse from %d onto %d; w = %f\n", Id_RA_local[i], j, w);
        active_indicators_local[i][j] = 0;

        std::vector<int>::iterator pos = std::find(active_synapses_local[i].begin(),
                                                    active_synapses_local[i].end(), j);

        if (pos!= active_synapses_local[i].end())
            active_synapses_local[i].erase(pos);
        else
            std::cout << "Active synapse " << Id_RA_local[i] << " -> " << j << " to be erased is not found!" << std::endl;

    }
}

void NetworkGrowthSimulator::update_synapse(int i, int j)
{
	double w = weights_RA_RA_local[i][j];

    if ( ( w >= synaptic_params.ACTIVATION_THRESHOLD ) && ( active_indicators_local[i][j] == 0 ) && ( remodeled_local[i] == 0) )
    {
       // printf("Activated synapse from %d onto %d; weight = %f\n", Id_RA_local[i], j, w);
        active_indicators_local[i][j] = 1;
        active_synapses_local[i].push_back(j);
    }

    if ( (w >= synaptic_params.SUPERSYNAPSE_THRESHOLD) && (supersynapses_indicators_local[i][j] == 0) && (static_cast<int>(supersynapses_local[i].size()) < synaptic_params.Nss) )
    {
       // printf("Activated supersynapse from %d onto %d; weight = %f\n", Id_RA_local[i], j, w);
        supersynapses_indicators_local[i][j] = 1;
        supersynapses_local[i].push_back(j);
    }

    if ( (w < synaptic_params.SUPERSYNAPSE_THRESHOLD) && (supersynapses_indicators_local[i][j] == 1) )
    {
       // printf("Deactivated supersynapse from %d onto %d; weight = %f\n", Id_RA_local[i], j, w);
        supersynapses_indicators_local[i][j] = 0;
        remodeled_local[i] = 0;
        std::vector<int>::iterator pos = std::find(supersynapses_local[i].begin(),
                                                    supersynapses_local[i].end(), j);

        if (pos!= supersynapses_local[i].end())
            supersynapses_local[i].erase(pos);
        else
            std::cout << "Supersynapse " << Id_RA_local[i] << " -> " << j << " to be erased is not found!" << std::endl;
    }

    if ( (w < synaptic_params.ACTIVATION_THRESHOLD) && (active_indicators_local[i][j] == 1) )
    {
      //  printf("Deactivated synapse from %d onto %d; w = %f\n", Id_RA_local[i], j, w);
        active_indicators_local[i][j] = 0;

        std::vector<int>::iterator pos = std::find(active_synapses_local[i].begin(),
                                                    active_synapses_local[i].end(), j);

        if (pos!= active_synapses_local[i].end())
            active_synapses_local[i].erase(pos);
        else
            std::cout << "Active synapse " << Id_RA_local[i] << " -> " << j << " to be erased is not found!" << std::endl;

    }
}

void NetworkGrowthSimulator::LTP_toRescaled(double &w, double dt)
{
	// dt in LTP is postsynaptic spike time - presynaptic spike time  
	if (dt <= synaptic_params.T_0)
	{
		std::cerr << "Time t = " << dt << " is smaller than T_0 in LTP!" << std::endl;
		return;
	}
	
	if (dt <= synaptic_params.T_0 + synaptic_params.T_P)
		w = w + MATURE_SYNAPSE_SCALING * synaptic_params.A_P * (dt - synaptic_params.T_0) / synaptic_params.T_P;
	else
		w = w + MATURE_SYNAPSE_SCALING * synaptic_params.A_P * exp(-(dt - synaptic_params.T_0 - synaptic_params.T_P) / synaptic_params.TAU_P);
		
	if (w > MATURE_SYNAPSE_SCALING * synaptic_params.WEIGHT_MAX)
        w = MATURE_SYNAPSE_SCALING * synaptic_params.WEIGHT_MAX;

    if (w < 0)
        w = 0;
}

void NetworkGrowthSimulator::LTP(double &w, double dt)
{
	// dt in LTP is postsynaptic spike time - presynaptic spike time  
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

void NetworkGrowthSimulator::LTD(double &w, double dt)
{
	// dt in LTD is postsynaptic spike time - presynaptic spike time  
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

//~ void NetworkGrowthSimulator::LTP_burst(double &w, double dt)
//~ {
	//~ // dt in LTP_burst is postsynaptic burst time - presynaptic burst time  
	//~ if (dt <= T_0)
	//~ {
		//~ std::cerr << "Time t = " << dt << " is smaller than T_0 in LTP!" << std::endl;
		//~ return;
	//~ }
	//~ 
	//~ if (dt <= T_0 + T_P)
		//~ w = w + A_P * (dt - T_0) / T_P;
	//~ else
		//~ w = w + A_P * exp(-(dt - T_0 - T_P) / TAU_P);
		//~ 
	//~ if (w > WEIGHT_MAX)
        //~ w = WEIGHT_MAX;
//~ 
    //~ if (w < 0)
        //~ w = 0;
//~ }
//~ 
//~ 
//~ void NetworkGrowthSimulator::LTD_burst(double &w, double dt)
//~ {
	//~ // dt in LTD_burst is postsynaptic burst time - presynaptic burst time  
    //~ if (dt > T_0)
    //~ {
		//~ std::cerr << "Time in LTD dt = " << dt << " is bigger than T0 = " << T_0 << std::endl;
		//~ return;
	//~ }
	//~ 
	//~ if (dt >= T_0 - T_D)
		//~ w = w + w * A_D * (dt - T_0) / T_D;
	//~ else
		//~ w = w - w * A_D * exp((dt - T_0 + T_D) / TAU_D);
		//~ 
	//~ if (w < 0)
		//~ w = 0;
//~ }

//~ void NetworkGrowthSimulator::gather_mature_data(std::vector<std::vector<double>>& average_dendritic_spike_time) 
//~ {
//~ 
    //~ MPI_Status status;
    //~ // Gather all data to master process
    //~ int *spike_num_soma_local = new int[N_RA_local];
    //~ int *spike_num_soma_global = new int[N_RA];
    //~ int *spike_num_dend_local = new int[N_RA_local];
    //~ int *spike_num_dend_global = new int[N_RA];
    //~ int *spike_num_interneuron_local = new int[N_I_local];
    //~ int *spike_num_interneuron_global = new int[N_I];
//~ 
    //~ int *recvcounts_RA = new int[MPI_size];
    //~ int *displs_RA = new int[MPI_size];
	//~ int *recvcounts_I = new int[MPI_size];
    //~ int *displs_I = new int[MPI_size];
//~ 
//~ 
    //~ if (MPI_rank == 0)
    //~ {
        //~ recvcounts_RA[0] = N_RA_sizes[0];
        //~ displs_RA[0] = 0;
//~ 
        //~ for (int i = 1; i < MPI_size; i++)
        //~ {
            //~ recvcounts_RA[i] = N_RA_sizes[i];
            //~ displs_RA[i] = displs_RA[i-1] + recvcounts_RA[i-1];
        //~ }
		//~ 
		//~ recvcounts_I[0] = N_I_sizes[0];
        //~ displs_I[0] = 0;
//~ 
        //~ for (int i = 1; i < MPI_size; i++)
        //~ {
            //~ recvcounts_I[i] = N_I_sizes[i];
            //~ displs_I[i] = displs_I[i-1] + recvcounts_I[i-1];
        //~ }
//~ 
    //~ }
//~ 
    //~ for (int i = 0; i < N_RA_local; i++)
    //~ {
        //~ spike_num_soma_local[i] = spikes_in_trial_soma_local[i].size();
        //~ spike_num_dend_local[i] = spikes_in_trial_dend_local[i].size();
//~ 
        //~ //printf("Rank = %d, supersyn_sizes_local[%d] = %d\n", MPI_rank, Id_RA_local[i], supersyn_sizes_local[i]);
        //~ //printf("Rank = %d, syn_sizes_local[%d] = %d\n", MPI_rank, Id_RA_local[i], syn_sizes_local[i]);
        //~ //printf("Rank = %d, spike_num_local[%d] = %d\n", MPI_rank, Id_RA_local[i], spike_num_local[i]);
    //~ }
//~ 
	//~ for (int i = 0; i < N_I_local; i++)
		//~ spike_num_interneuron_local[i] = (int) spikes_in_trial_interneuron_local[i].size();
//~ 
	//~ MPI_Gatherv(&spike_num_interneuron_local[0], N_I_local, MPI_INT,
        //~ &spike_num_interneuron_global[0], recvcounts_I, displs_I, MPI_INT, 0, MPI_COMM_WORLD);
//~ 
    //~ MPI_Gatherv(&spike_num_soma_local[0], N_RA_local, MPI_INT,
        //~ &spike_num_soma_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);
//~ 
    //~ MPI_Gatherv(&spike_num_dend_local[0], N_RA_local, MPI_INT,
        //~ &spike_num_dend_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);
//~ 
    //~ 
	//~ if (MPI_rank == 0)
    //~ {
        //~ for (int i = 0; i < N_RA; i++)
        //~ {
            //~ spikes_in_trial_soma_global[Id_RA_global[i]].resize(spike_num_soma_global[i]);
            //~ spikes_in_trial_dend_global[Id_RA_global[i]].resize(spike_num_dend_global[i]);
//~ 
            //~ //printf("Master, spike_num_global[%d] = %d\n", i, spike_num_global[i]);
        //~ }
//~ 
//~ 
        //~ // Copy from master's local arrays
        //~ for (int i = 0; i < N_RA_local; i++)
        //~ {
            //~ spikes_in_trial_soma_global[Id_RA_global[i]] = spikes_in_trial_soma_local[i];
            //~ spikes_in_trial_dend_global[Id_RA_global[i]] = spikes_in_trial_dend_local[i];
        //~ }
//~ 
			//~ 
    //~ // Gather from others
		//~ int N = N_RA_sizes[0]; // number of RA neurons in the processes with lower rank
//~ 
        //~ for (int i = 1; i < MPI_size; i++)
        //~ {
//~ 
            //~ for (int j = 0; j < N_RA_sizes[i]; j++)
            //~ {
                //~ int count;
				//~ int index_in_global_array = N + j;
                //~ int receive_index = Id_RA_global[index_in_global_array];
                //~ 
                //~ if (spike_num_soma_global[index_in_global_array] != 0)
                //~ {
                    //~ MPI_Recv(&spikes_in_trial_soma_global[receive_index][0],
                        //~ spike_num_soma_global[index_in_global_array], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
//~ 
                    //~ MPI_Get_count(&status, MPI_INT, &count);
                    //~ //printf("Recv spikes in trial; from i = %d  count = %d\n", i, count);
                //~ }
//~ 
                //~ if (spike_num_dend_global[index_in_global_array] != 0)
                //~ {
                    //~ MPI_Recv(&spikes_in_trial_dend_global[receive_index][0],
                        //~ spike_num_dend_global[index_in_global_array], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
//~ 
                    //~ MPI_Get_count(&status, MPI_INT, &count);
                    //~ //printf("Recv spikes in trial; from i = %d  count = %d\n", i, count);
                //~ }
//~ 
               //~ // for (int k = 0; k < spikes_in_trial_global[N_RA_local + (i-1)*offset + j].size(); k++)
                   //~ // printf("Master; spikes_in_trial_global[%d][%d] = %f\n", N_RA_local + (i-1)*offset + j, k,
                      //~ //  spikes_in_trial_global[N_RA_local + (i-1)*offset + j][k]);
            //~ }
//~ 
			//~ N += N_RA_sizes[i];
//~ 
        //~ }
//~ 
    //~ }
//~ 
    //~ else
    //~ {
        //~ for (int i = 0; i < N_RA_local; i++)
        //~ {
//~ 
            //~ if (spike_num_soma_local[i] != 0)
                //~ MPI_Send(&spikes_in_trial_soma_local[i][0],
                        //~ spike_num_soma_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
//~ 
            //~ if (spike_num_dend_local[i] != 0)
                //~ MPI_Send(&spikes_in_trial_dend_local[i][0],
                        //~ spike_num_dend_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        //~ }
    //~ }
//~ 
	//~ // gather spikes of interneurons
	//~ if (MPI_rank == 0)
    //~ {
    //~ 
        //~ for (int i = 0; i < N_I; i++)
			//~ spikes_in_trial_interneuron_global[i].resize(spike_num_interneuron_global[i]);
        //~ 
	    //~ for (int i = 0; i < N_I_local; i++)
			//~ spikes_in_trial_interneuron_global[i] = spikes_in_trial_interneuron_local[i];
	//~ 
        //~ int N = N_I_sizes[0]; // number of I neurons in the processes with lower rank
//~ 
        //~ for (int i = 1; i < MPI_size; i++)
        //~ {
            //~ for (int j = 0; j < N_I_sizes[i]; j++)
            //~ {
                //~ int count;
				//~ int receive_index = N + j;
//~ 
                //~ if (spike_num_interneuron_global[receive_index] != 0)
                //~ {
                    //~ MPI_Recv(&spikes_in_trial_interneuron_global[receive_index][0],
                        //~ spike_num_interneuron_global[receive_index], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
//~ 
                    //~ MPI_Get_count(&status, MPI_INT, &count);
                    //~ //printf("Recv spikes in trial; from i = %d  count = %d\n", i, count);
                //~ }
//~ 
               //~ // for (int k = 0; k < spikes_in_trial_global[N_RA_local + (i-1)*offset + j].size(); k++)
                   //~ // printf("Master; spikes_in_trial_global[%d][%d] = %f\n", N_RA_local + (i-1)*offset + j, k,
                      //~ //  spikes_in_trial_global[N_RA_local + (i-1)*offset + j][k]);
            //~ }
//~ 
			//~ N += N_I_sizes[i];
        //~ }
    //~ }
//~ 
    //~ 
    //~ else
    //~ {
        //~ for (int i = 0; i < N_I_local; i++)
        //~ {
            //~ if (spike_num_interneuron_local[i] != 0)
                //~ MPI_Send(&spikes_in_trial_interneuron_local[i][0],
                    //~ spike_num_interneuron_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
//~ 
//~ 
        //~ }
    //~ }
    //~ 
	//~ // process dendritic spikes
	//~ if (MPI_rank == 0)
	//~ {
		//~ for (int i = 0; i < N_RA; i++)
		//~ {
			//~ if (spike_num_dend_global[i] > 0)
            //~ {
                //~ double average_spike_time = std::accumulate(spikes_in_trial_dend_global[i].begin(), spikes_in_trial_dend_global[i].end(), 0.0) / static_cast<double>(spike_num_dend_global[i]);
//~ 
                //~ //printf("Average dendritic spike time = %f\n", average_spike_time);
				//~ average_dendritic_spike_time[i].push_back(average_spike_time);
            //~ }
		//~ }
//~ 
//~ 
	//~ }
//~ 
//~ 
    //~ delete [] recvcounts_RA;
    //~ delete [] displs_RA;
    //~ delete [] recvcounts_I;
    //~ delete [] displs_I;
    //~ delete [] spike_num_soma_local;
    //~ delete [] spike_num_soma_global;
    //~ delete [] spike_num_dend_local;
    //~ delete [] spike_num_dend_global;
	//~ delete [] spike_num_interneuron_local;
	//~ delete [] spike_num_interneuron_global;
//~ }

void NetworkGrowthSimulator::gather_graph_state_data()
{
    MPI_Status status;
    // Gather all data to master process
    int *supersyn_sizes_global = new int[N_RA];
    int *supersyn_sizes_local = new int[N_RA_local];
    int *syn_sizes_global = new int[N_RA];
    int *syn_sizes_local = new int[N_RA_local];
    int *spike_num_soma_local = new int[N_RA_local];
    int *spike_num_soma_global = new int[N_RA];
    int *spike_num_dend_local = new int[N_RA_local];
    int *spike_num_dend_global = new int[N_RA];
    int *spike_num_interneuron_local = new int[N_I_local];
    int *spike_num_interneuron_global = new int[N_I];
    
   
    int *recvcounts_RA = new int[MPI_size];
    int *displs_RA = new int[MPI_size];
	int *recvcounts_I = new int[MPI_size];
    int *displs_I = new int[MPI_size];


    if (MPI_rank == 0)
    {
        recvcounts_RA[0] = N_RA_sizes[0];
        displs_RA[0] = 0;

        for (int i = 1; i < MPI_size; i++)
        {
            recvcounts_RA[i] = N_RA_sizes[i];
            displs_RA[i] = displs_RA[i-1] + recvcounts_RA[i-1];
        }
		
		recvcounts_I[0] = N_I_sizes[0];
        displs_I[0] = 0;

        for (int i = 1; i < MPI_size; i++)
        {
            recvcounts_I[i] = N_I_sizes[i];
            displs_I[i] = displs_I[i-1] + recvcounts_I[i-1];
        }

    }

    for (int i = 0; i < N_RA_local; i++)
    {
        supersyn_sizes_local[i] = supersynapses_local[i].size();
        syn_sizes_local[i] = active_synapses_local[i].size();
        spike_num_soma_local[i] = spikes_in_trial_soma_local[i].size();
        spike_num_dend_local[i] = spikes_in_trial_dend_local[i].size();

        //printf("Rank = %d, supersyn_sizes_local[%d] = %d\n", MPI_rank, Id_RA_local[i], supersyn_sizes_local[i]);
        //printf("Rank = %d, syn_sizes_local[%d] = %d\n", MPI_rank, Id_RA_local[i], syn_sizes_local[i]);
        //printf("Rank = %d, spike_num_soma_local[%d] = %d\n", MPI_rank, Id_RA_local[i], spike_num_soma_local[i]);
        //printf("Rank = %d, spike_num_dend_local[%d] = %d\n", MPI_rank, Id_RA_local[i], spike_num_dend_local[i]);
    }

    // all gatherv functions collect data from local processes. Data is a concatenated version of local 
    // data on the processes. Therefore neuronal ids may not be monotonic!!! Data is rearranged only 
    // writing to files

	for (int i = 0; i < N_I_local; i++)
		spike_num_interneuron_local[i] = (int) spikes_in_trial_interneuron_local[i].size();

    MPI_Gatherv(&supersyn_sizes_local[0], N_RA_local, MPI_INT,
        &supersyn_sizes_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gatherv(&syn_sizes_local[0], N_RA_local, MPI_INT,
        &syn_sizes_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);
	
	MPI_Gatherv(&spike_num_interneuron_local[0], N_I_local, MPI_INT,
        &spike_num_interneuron_global[0], recvcounts_I, displs_I, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gatherv(&spike_num_soma_local[0], N_RA_local, MPI_INT,
        &spike_num_soma_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gatherv(&spike_num_dend_local[0], N_RA_local, MPI_INT,
        &spike_num_dend_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);

	
	
	
    // Receive functions on the other hand rearrange data immediately. Thus, there is no need
    // to take special care while writing to files

    if (MPI_rank == 0)
    {
        for (int i = 0; i < N_RA; i++)
        {
            //printf("Master; supersyn_sizes_global[%d] = %d\n", i, supersyn_sizes_global[i]);
            supersynapses_global[Id_RA_global[i]].resize(supersyn_sizes_global[i]);
            active_synapses_global[Id_RA_global[i]].resize(syn_sizes_global[i]);
            spikes_in_trial_soma_global[Id_RA_global[i]].resize(spike_num_soma_global[i]);
            spikes_in_trial_dend_global[Id_RA_global[i]].resize(spike_num_dend_global[i]);

            //printf("Master, supersyn_sizes_global[%d] = %d\n", i, supersyn_sizes_global[i]);
            //printf("Master, syn_sizes_global[%d] = %d\n", i, syn_sizes_global[i]);
            //printf("Master, spike_num_global[%d] = %d\n", i, spike_num_global[i]);
        }


        // Copy from master's local arrays to master's global arrays
        // Note that data is copied according to neuronal ids stored in Id_RA_global
        // thus it is already ordered according to neuronal ids.
        for (int i = 0; i < N_RA_local; i++)
        {
			
			
            supersynapses_global[Id_RA_global[i]] = supersynapses_local[i];
            active_synapses_global[Id_RA_global[i]] = active_synapses_local[i];
            spikes_in_trial_soma_global[Id_RA_global[i]] = spikes_in_trial_soma_local[i];
            spikes_in_trial_dend_global[Id_RA_global[i]] = spikes_in_trial_dend_local[i];

            //for (int k = 0; k < active_supersynapses_global[i].size(); k++)
              //      printf("Master; active_supersynapses_global[%d][%d] = %u\n", i, k,
                //       active_supersynapses_global[i][k]);

            //for (int k = 0; k < active_synapses_global[i].size(); k++)
              //      printf("Master; active_synapses_global[%d][%d] = %u\n", i, k,
                //        active_synapses_global[i][k]);
        }

    // Receive from others
		int N = N_RA_sizes[0]; // number of RA neurons in the processes with lower rank

        for (int i = 1; i < MPI_size; i++)
        {

            for (int j = 0; j < N_RA_sizes[i]; j++)
            {
                int count;
				int index_in_global_array = N + j;
                int receive_index = Id_RA_global[index_in_global_array];

                //printf("Recv weights; from i = %d  count = %d\n", i, count);
                            
                if (supersyn_sizes_global[index_in_global_array] != 0)
                {
                    MPI_Recv(&supersynapses_global[receive_index][0],
                        supersyn_sizes_global[index_in_global_array], MPI_INT, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv supersynapses; from i = %d  count = %d\n", i, count);
                }

                if (syn_sizes_global[index_in_global_array] != 0)
                {
                    MPI_Recv(&active_synapses_global[receive_index][0],
                        syn_sizes_global[index_in_global_array], MPI_INT, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv synapses; from i = %d  count = %d\n", i, count);
                }

                if (spike_num_soma_global[index_in_global_array] != 0)
                {
                    MPI_Recv(&spikes_in_trial_soma_global[receive_index][0],
                        spike_num_soma_global[index_in_global_array], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv spikes in trial; from i = %d  count = %d\n", i, count);
                }

                if (spike_num_dend_global[index_in_global_array] != 0)
                {
                    MPI_Recv(&spikes_in_trial_dend_global[receive_index][0],
                        spike_num_dend_global[index_in_global_array], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv spikes in trial; from i = %d  count = %d\n", i, count);
                }
                //printf("Master;  N_RA_local + (i-1)*offset + j = %d\n", N_RA_local + (i-1)*offset + j);
                //printf("Master;  supersyn_sizes_global[N_RA_local + (i-1)*offset + j] = %d\n",
                 //   supersyn_sizes_global[N_RA_local + (i-1)*offset + j]);

                //for (int k = 0; k < N_RA; k++)
                  // printf("Master; weights_global[%d][%d] = %f\n", N_RA_local + (i-1)*offset + j, k,
                    //    weights_global[N_RA_local + (i-1)*offset + j][k]);

                //for (int k = 0; k < active_supersynapses_global[N_RA_local + (i-1)*offset + j].size(); k++)
                  //  printf("Master; active_supersynapses_global[%d][%d] = %u\n", N_RA_local + (i-1)*offset + j, k,
                    //   active_supersynapses_global[N_RA_local + (i-1)*offset + j][k]);

                //for (int k = 0; k < active_synapses_global[N_RA_local + (i-1)*offset + j].size(); k++)
                  //  printf("Master; active_synapses_global[%d][%d] = %u\n", N_RA_local + (i-1)*offset + j, k,
                    //    active_synapses_global[N_RA_local + (i-1)*offset + j][k]);

               // for (int k = 0; k < spikes_in_trial_global[N_RA_local + (i-1)*offset + j].size(); k++)
                   // printf("Master; spikes_in_trial_global[%d][%d] = %f\n", N_RA_local + (i-1)*offset + j, k,
                      //  spikes_in_trial_global[N_RA_local + (i-1)*offset + j][k]);
            }

			N += N_RA_sizes[i];

        }

    }

    else
    {
        for (int i = 0; i < N_RA_local; i++)
        {
            //for (int j = 0; j < N_RA; j++)
              //  printf("Rank %d; weights_local[%d][%d] = %f\n", MPI_rank, i, j, weights_local[i][j]);
            if (supersyn_sizes_local[i] != 0)
                MPI_Send(&supersynapses_local[i][0],
                        supersyn_sizes_local[i], MPI_INT, 0, 0, MPI_COMM_WORLD);

            if (syn_sizes_local[i] != 0)
                MPI_Send(&active_synapses_local[i][0],
                        syn_sizes_local[i], MPI_INT, 0, 0, MPI_COMM_WORLD);

            if (spike_num_soma_local[i] != 0)
                MPI_Send(&spikes_in_trial_soma_local[i][0],
                        spike_num_soma_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

            if (spike_num_dend_local[i] != 0)
                MPI_Send(&spikes_in_trial_dend_local[i][0],
                        spike_num_dend_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            //printf("Rank %d;  supersyn_sizes_local[i] = %d\n", MPI_rank, supersyn_sizes_local[i]);
            //for (int k = 0; k < active_synapses_local[i].size(); k++)
             //   printf("Rank %d; active_synapses_local[%d][%d] = %u\n", MPI_rank,
              //      Id_RA_local[i], k, active_synapses_local[i][k]);

           // for (int k = 0; k < active_supersynapses_local[i].size(); k++)
               // printf("Rank %d; active_supersynapses_local[%d][%d] = %u\n", MPI_rank,
                 //   Id_RA_local[i], k, active_supersynapses_local[i][k]);


        }
    }

	// gather spikes of interneurons
	if (MPI_rank == 0)
    {
    
        for (int i = 0; i < N_I; i++)
			spikes_in_trial_interneuron_global[i].resize(spike_num_interneuron_global[i]);
        
	    for (int i = 0; i < N_I_local; i++)
			spikes_in_trial_interneuron_global[i] = spikes_in_trial_interneuron_local[i];
	
        int N = N_I_sizes[0]; // number of I neurons in the processes with lower rank

        for (int i = 1; i < MPI_size; i++)
        {
            for (int j = 0; j < N_I_sizes[i]; j++)
            {
                int count;
				int receive_index = N + j;

                if (spike_num_interneuron_global[receive_index] != 0)
                {
                    MPI_Recv(&spikes_in_trial_interneuron_global[receive_index][0],
                        spike_num_interneuron_global[receive_index], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

                    MPI_Get_count(&status, MPI_INT, &count);
                    //printf("Recv spikes in trial; from i = %d  count = %d\n", i, count);
                }

               // for (int k = 0; k < spikes_in_trial_global[N_RA_local + (i-1)*offset + j].size(); k++)
                   // printf("Master; spikes_in_trial_global[%d][%d] = %f\n", N_RA_local + (i-1)*offset + j, k,
                      //  spikes_in_trial_global[N_RA_local + (i-1)*offset + j][k]);
            }

			N += N_I_sizes[i];
        }
    }

    
    else
    {
        for (int i = 0; i < N_I_local; i++)
        {
            if (spike_num_interneuron_local[i] != 0)
                MPI_Send(&spikes_in_trial_interneuron_local[i][0],
                    spike_num_interneuron_local[i], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

           // for (int k = 0; k < active_supersynapses_local[i].size(); k++)
               // printf("Rank %d; active_supersynapses_local[%d][%d] = %u\n", MPI_rank,
                 //   Id_RA_local[i], k, active_supersynapses_local[i][k]);


        }
    }
    
    delete [] recvcounts_RA;
    delete [] displs_RA;
    delete [] recvcounts_I;
    delete [] displs_I;
    delete [] supersyn_sizes_global;
    delete [] supersyn_sizes_local;
    delete [] syn_sizes_global;
    delete [] syn_sizes_local;
    delete [] spike_num_soma_local;
    delete [] spike_num_soma_global;
    delete [] spike_num_dend_local;
    delete [] spike_num_dend_global;
	delete [] spike_num_interneuron_local;
	delete [] spike_num_interneuron_global;
}

void NetworkGrowthSimulator::gather_inhibition()
{
	MPI_Status status;
	
	if (MPI_rank == 0)
	{
		for (int i = 0; i < N_RA_local; i++)
			Ginh_global[i] = Ginh_local[i];

		// Receive from others
		int N = N_RA_sizes[0]; // number of RA neurons in the processes with lower rank

		for (int i = 1; i < MPI_size; i++)
		{
			for (int j = 0; j < N_RA_sizes[i]; j++)
			{
				int count;
				int receive_index = N + j;
				
				//std::cout << "receive_index = " << receive_index << std::endl;
				//std::cout << "Size of global array with num spikes: " << num_spikes_in_recent_trials_global[receive_index].size() << std::endl;
			

				MPI_Recv(&Ginh_global[receive_index][0],
										static_cast<int>(Ginh_global[receive_index].size()), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
				MPI_Get_count(&status, MPI_DOUBLE, &count);

			}

			N += N_RA_sizes[i];
		}

    }
    else
    {
        for (int i = 0; i < N_RA_local; i++)
            MPI_Send(&Ginh_local[i][0], static_cast<int>(Ginh_local[i].size()), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    
}


void NetworkGrowthSimulator::gather_full_state_data()
{
    MPI_Status status;
    
    // Gather all data to master process
    int *recvcounts_RA = new int[MPI_size];
    int *displs_RA = new int[MPI_size];
	int *recvcounts_I = new int[MPI_size];
    int *displs_I = new int[MPI_size];


    if (MPI_rank == 0)
    {
        recvcounts_RA[0] = N_RA_sizes[0];
        displs_RA[0] = 0;

        for (int i = 1; i < MPI_size; i++)
        {
            recvcounts_RA[i] = N_RA_sizes[i];
            displs_RA[i] = displs_RA[i-1] + recvcounts_RA[i-1];
        }
		
		recvcounts_I[0] = N_I_sizes[0];
        displs_I[0] = 0;

        for (int i = 1; i < MPI_size; i++)
        {
            recvcounts_I[i] = N_I_sizes[i];
            displs_I[i] = displs_I[i-1] + recvcounts_I[i-1];
        }

    }
	
	// gather replacement history   
    MPI_Gatherv(&num_trials_after_replacement_local[0], N_RA_local, MPI_INT,
        &num_trials_after_replacement_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);
    
    MPI_Gatherv(&remodeled_local[0], N_RA_local, MPI_INT,
        &remodeled_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);
    
    
    MPI_Gatherv(&rest_potential_local[0], N_RA_local, MPI_DOUBLE,
        &rest_potential_global[0], recvcounts_RA, displs_RA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	MPI_Gatherv(&GCa_local[0], N_RA_local, MPI_DOUBLE,
        &GCa_global[0], recvcounts_RA, displs_RA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	
	//MPI_Gatherv(&gaba_potential_local[0], N_RA_local, MPI_DOUBLE,
     //   &gaba_potential_global[0], recvcounts_RA, displs_RA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	//MPI_Gatherv(&firing_rate_short_local[0], N_RA_local, MPI_DOUBLE,
      //  &firing_rate_short_global[0], recvcounts_RA, displs_RA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
    MPI_Gatherv(&firing_rate_long_local[0], N_RA_local, MPI_DOUBLE,
        &firing_rate_long_global[0], recvcounts_RA, displs_RA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
	MPI_Gatherv(&mature_local[0], N_RA_local, MPI_INT,
        &mature_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);
    
    MPI_Gatherv(&maturation_scale_local[0], N_RA_local, MPI_INT,
        &maturation_scale_global[0], recvcounts_RA, displs_RA, MPI_INT, 0, MPI_COMM_WORLD);
     
    // Receive functions on the other hand rearrange data immediately. Thus, there is no need
    // to take special care while writing to files

    if (MPI_rank == 0)
    {
		for (int i = 0; i < N_RA_local; i++)
        {
            for (int j = 0; j < N_RA; j++)
            {
               weights_RA_RA_global[Id_RA_global[i]][j] = weights_RA_RA_local[i][j];

               //printf("Master; weights_global[%d][%d] = %f\n", i, j,
                 //   weights_global[i][j]);
            }
            
            for (int j = 0; j < maturation_params.RATE_WINDOW_LONG; j++)
				num_spikes_in_recent_trials_global[Id_RA_global[i]][j] = num_spikes_in_recent_trials_local[i][j];

        }
    // Receive from others
		int N = N_RA_sizes[0]; // number of RA neurons in the processes with lower rank

        for (int i = 1; i < MPI_size; i++)
        {

            for (int j = 0; j < N_RA_sizes[i]; j++)
            {
                int count;
				int index_in_global_array = N + j;
                int receive_index = Id_RA_global[index_in_global_array];
                
                //std::cout << "receive_index = " << receive_index << std::endl;
                //std::cout << "Size of global array with num spikes: " << num_spikes_in_recent_trials_global[receive_index].size() << std::endl;
			

                MPI_Recv(&weights_RA_RA_global[receive_index][0],
                                        N_RA, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_DOUBLE, &count);
                //printf("Recv weights; from i = %d  count = %d\n", i, count);
                
                MPI_Recv(&num_spikes_in_recent_trials_global[receive_index][0],
                                        maturation_params.RATE_WINDOW_LONG, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
                                        
                MPI_Get_count(&status, MPI_INT, &count);
                
               
                //printf("Master;  N_RA_local + (i-1)*offset + j = %d\n", N_RA_local + (i-1)*offset + j);
                //printf("Master;  supersyn_sizes_global[N_RA_local + (i-1)*offset + j] = %d\n",
                 //   supersyn_sizes_global[N_RA_local + (i-1)*offset + j]);

                //for (int k = 0; k < N_RA; k++)
                  // printf("Master; weights_global[%d][%d] = %f\n", N_RA_local + (i-1)*offset + j, k,
                    //    weights_global[N_RA_local + (i-1)*offset + j][k]);

                //for (int k = 0; k < active_supersynapses_global[N_RA_local + (i-1)*offset + j].size(); k++)
                  //  printf("Master; active_supersynapses_global[%d][%d] = %u\n", N_RA_local + (i-1)*offset + j, k,
                    //   active_supersynapses_global[N_RA_local + (i-1)*offset + j][k]);

                //for (int k = 0; k < active_synapses_global[N_RA_local + (i-1)*offset + j].size(); k++)
                  //  printf("Master; active_synapses_global[%d][%d] = %u\n", N_RA_local + (i-1)*offset + j, k,
                    //    active_synapses_global[N_RA_local + (i-1)*offset + j][k]);

               // for (int k = 0; k < spikes_in_trial_global[N_RA_local + (i-1)*offset + j].size(); k++)
                   // printf("Master; spikes_in_trial_global[%d][%d] = %f\n", N_RA_local + (i-1)*offset + j, k,
                      //  spikes_in_trial_global[N_RA_local + (i-1)*offset + j][k]);
            }

			N += N_RA_sizes[i];

        }

    }

    else
    {
        for (int i = 0; i < N_RA_local; i++)
        {
			
			//~ if (Id_RA_local[i] == 100)
			//~ {
				//~ std::cout << "Num recent spikes for neuron " << Id_RA_local[i] << " in local array:\n";
				//~ for (int j = 0; j < maturation_params.RATE_WINDOW_LONG; j++)
					//~ std::cout << num_spikes_in_recent_trials_local[i][j] << " ";
				//~ std::cout << std::endl;
				//~ 
				//~ std::cout << "\nWeights for neuron " << Id_RA_local[i] << " in local array:\n";
				//~ for (int j = 0; j < N_RA; j++)
					//~ std::cout << weights_RA_RA_local[i][j] << " ";
				//~ std::cout << std::endl;
				//~ 
			//~ }
			
			//std::cout << "Rank " << MPI_rank << " rate_window = " << maturation_params.RATE_WINDOW_LONG << std::endl;
			
            //for (int j = 0; j < N_RA; j++)
              //  printf("Rank %d; weights_local[%d][%d] = %f\n", MPI_rank, i, j, weights_local[i][j]);
            MPI_Send(&weights_RA_RA_local[i][0], N_RA, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			
			//std::cout << "Size of local array with num spikes: " << num_spikes_in_recent_trials_local[i].size() << std::endl;
			
			std::vector<int> num_spikes(num_spikes_in_recent_trials_local[i].size());
			std::copy(num_spikes_in_recent_trials_local[i].begin(), num_spikes_in_recent_trials_local[i].end(), num_spikes.begin());
			
			//MPI_Send(&num_spikes_in_recent_trials_local[i][0],
              //                          maturation_params.RATE_WINDOW_LONG, MPI_INT, 0, 1, MPI_COMM_WORLD);
			
			MPI_Send(&num_spikes[0],
                                        maturation_params.RATE_WINDOW_LONG, MPI_INT, 0, 1, MPI_COMM_WORLD);
			
			
            //printf("Rank %d;  supersyn_sizes_local[i] = %d\n", MPI_rank, supersyn_sizes_local[i]);
            //for (int k = 0; k < active_synapses_local[i].size(); k++)
             //   printf("Rank %d; active_synapses_local[%d][%d] = %u\n", MPI_rank,
              //      Id_RA_local[i], k, active_synapses_local[i][k]);

           // for (int k = 0; k < active_supersynapses_local[i].size(); k++)
               // printf("Rank %d; active_supersynapses_local[%d][%d] = %u\n", MPI_rank,
                 //   Id_RA_local[i], k, active_supersynapses_local[i][k]);


        }
    }
    
    //~ if (MPI_rank == 0)
    //~ {
		//~ int neuron_id = 100;
		//~ std::cout << "\nNum recent spikes for neuron " << neuron_id << " in global array:\n";
		//~ for (int j = 0; j < maturation_params.RATE_WINDOW_LONG; j++)
			//~ std::cout << num_spikes_in_recent_trials_global[neuron_id][j] << " ";
		//~ std::cout << std::endl;
		//~ 
		//~ std::cout << "\nWeights for neuron " << neuron_id << " in global array:\n";
		//~ for (int j = 0; j < N_RA; j++)
			//~ std::cout << weights_RA_RA_global[neuron_id][j] << " ";
		//~ std::cout << std::endl;
		//~ 
	//~ }

	
    delete [] recvcounts_RA;
    delete [] displs_RA;
    delete [] recvcounts_I;
    delete [] displs_I;
}

void NetworkGrowthSimulator::get_neuronRA_location(int n, int* rank, int* shift)
{
	int i = 0;
	int N = 0;

	while (n >= N)
	{
		if (n < N + N_RA_sizes[i])
		{
			*rank = i;
			*shift = n - N;
		}
		
		N += N_RA_sizes[i];
		i++;
	}
}

void NetworkGrowthSimulator::get_neuronI_location(int n, int* rank, int* shift)
{
	int i = 0;
	int N = 0;

	while (n >= N)
	{
		if (n < N + N_I_sizes[i])
		{
			*rank = i;
			*shift = n - N;
		}
		
		N += N_I_sizes[i];
		i++;
	}
}



//~ void NetworkGrowthSimulator::write_chain_test(int num_trials, std::vector<int>& num_trials_with_dend_spikes, std::vector<double>& average_num_dend_spikes_in_trials, 
                                    //~ std::vector<double>& average_num_soma_spikes_in_trials, std::vector<double>& mean_burst_time, 
                                    //~ std::vector<double>& std_burst_time, const char* filename)
//~ {
    //~ if (MPI_rank == 0)
    //~ {
        //~ std::ofstream out;
//~ 
        //~ out.open(filename, std::ios::out | std::ios::binary);
//~ 
        //~ // write number of neurons to each file
//~ 
        //~ out.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));
        //~ out.write(reinterpret_cast<char *>(&num_trials), sizeof(num_trials));
//~ 
        //~ for (int i = 0; i < N_RA; i++)
        //~ {
            //~ double firing_robustness = num_trials_with_dend_spikes[i] / static_cast<double>(num_trials);
//~ 
            //~ out.write(reinterpret_cast<char *>(&firing_robustness), sizeof(firing_robustness));
            //~ out.write(reinterpret_cast<char *>(&average_num_dend_spikes_in_trials[i]), sizeof(average_num_dend_spikes_in_trials[i]));
            //~ out.write(reinterpret_cast<char *>(&average_num_soma_spikes_in_trials[i]), sizeof(average_num_soma_spikes_in_trials[i]));
            //~ out.write(reinterpret_cast<char *>(&mean_burst_time[i]), sizeof(mean_burst_time[i]));
            //~ out.write(reinterpret_cast<char *>(&std_burst_time[i]), sizeof(std_burst_time[i]));
        //~ }
        //~ out.close();
    //~ }
//~ }

void NetworkGrowthSimulator::write_active_synapses(const char* filename)
{
	std::ofstream out;

	out.open(filename, std::ios::out | std::ios::binary);

	// write number of neurons to each file

	out.write(reinterpret_cast<char *>(&N_RA), sizeof(int));
	out.write(reinterpret_cast<char *>(&trial_number), sizeof(int));

	for (int i = 0; i < N_RA; i++)
	{
		// write number of connections to RA neurons
		int num_synapses = static_cast<int>(active_synapses_global[i].size());
		
		out.write(reinterpret_cast<char *>(&num_synapses), sizeof(int));

		for (int j = 0; j < num_synapses; j++)
			out.write(reinterpret_cast<char *>(&active_synapses_global[i][j]), sizeof(int));
	}
	out.close();
    
}

void NetworkGrowthSimulator::write_super_synapses(const char* filename)
{
	std::ofstream out;

	out.open(filename, std::ios::out | std::ios::binary);

	// write number of neurons to each file

	out.write(reinterpret_cast<char *>(&N_RA), sizeof(int));
	out.write(reinterpret_cast<char *>(&trial_number), sizeof(int));

	for (int i = 0; i < N_RA; i++)
	{
		// write number of connections to RA neurons
		int num_synapses = static_cast<int>(supersynapses_global[i].size());
		
		out.write(reinterpret_cast<char *>(&num_synapses), sizeof(int));

		for (int j = 0; j < num_synapses; j++)
			out.write(reinterpret_cast<char *>(&supersynapses_global[i][j]), sizeof(int));
	}
	out.close();
}

void NetworkGrowthSimulator::write_weight_statistics(const char* filename)
{

	std::ofstream out;

	if (trial_number == 0)
		out.open(filename, std::ios::binary | std::ios::out);
	else
		out.open(filename, std::ios::binary | std::ios::out | std::ios::app);

	out.write(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));

	double mean_synaptic_weight = 0.0;
	
	// calculate sum of all synaptic weights
	for (size_t i = 0; i < weights_RA_RA_global.size(); i++)
		mean_synaptic_weight += std::accumulate(weights_RA_RA_global[i].begin(), weights_RA_RA_global[i].end(), 0.0);

	// normalize sum to get mean synaptic weight
	if (static_cast<int>(weights_RA_RA_global.size() > 0))
	{
		if (static_cast<int>(weights_RA_RA_global[0].size() > 0))
			mean_synaptic_weight /= (static_cast<double>(weights_RA_RA_global.size()) * static_cast<double>(weights_RA_RA_global[0].size()));
		else
			std::cerr << "Matrix weights_global has no columns!" << std::endl;
	}
	else
		std::cerr << "Matrix weights_global has no rows!" << std::endl;
	
	double std_synaptic_weight = 0; // standard deviation of synaptic weights

	for (size_t i = 0; i < weights_RA_RA_global.size(); i++)
	{
		std::for_each(weights_RA_RA_global[i].begin(), weights_RA_RA_global[i].end(), [&std_synaptic_weight, mean_synaptic_weight](const double w)
		{
			std_synaptic_weight += (w - mean_synaptic_weight) * (w - mean_synaptic_weight);
		});
	}
	
	// if there are more than two weights in the weight matrix
	if (static_cast<int>(weights_RA_RA_global.size()) * static_cast<int>(weights_RA_RA_global[0].size()) > 1)
		std_synaptic_weight = sqrt(std_synaptic_weight / (static_cast<double>(weights_RA_RA_global.size()) * static_cast<double>(weights_RA_RA_global[0].size()) - 1));
	else
		std::cerr << "Matrix weights_global has less than two elements!" << std::endl;

	out.write(reinterpret_cast<char *>(&mean_synaptic_weight), sizeof(mean_synaptic_weight));
	out.write(reinterpret_cast<char *>(&std_synaptic_weight), sizeof(std_synaptic_weight));

	out.close();

}

void NetworkGrowthSimulator::write_remodeled_indicators(const char* filename)
{
	std::ofstream output;

	output.open(filename, std::ios::out | std::ios::binary);

	output.write(reinterpret_cast<char *>(&N_RA), sizeof(int));
	output.write(reinterpret_cast<char *>(&trial_number), sizeof(int));

	for (int i = 0; i < N_RA; i++)
	{
		//printf("weigths[%d][%d] = %1.10f\n", i, j, weights[i][j]);
		output.write(reinterpret_cast<char *>(&remodeled_global[i]), sizeof(int));

	}
	output.close();
}

void NetworkGrowthSimulator::write_maturation_properties(const char* filename)
{
	std::ofstream output;

	output.open(filename, std::ios::out | std::ios::binary);

	output.write(reinterpret_cast<char *>(&N_RA), sizeof(int));
	output.write(reinterpret_cast<char *>(&trial_number), sizeof(int));

	for (int i = 0; i < N_RA; i++)
	{
		//printf("weigths[%d][%d] = %1.10f\n", i, j, weights[i][j]);
		output.write(reinterpret_cast<char *>(&mature_global[i]), sizeof(int));
		output.write(reinterpret_cast<char *>(&maturation_scale_global[i]), sizeof(int));
		output.write(reinterpret_cast<char *>(&rest_potential_global[i]), sizeof(double));
		output.write(reinterpret_cast<char *>(&GCa_global[i]), sizeof(double));
	}
	output.close();
}


void NetworkGrowthSimulator::write_mature_indicators(const char* filename)
{
	std::ofstream output;

	output.open(filename, std::ios::out | std::ios::binary);

	output.write(reinterpret_cast<char *>(&N_RA), sizeof(int));
	output.write(reinterpret_cast<char *>(&trial_number), sizeof(int));

	for (int i = 0; i < N_RA; i++)
	{
		//printf("weigths[%d][%d] = %1.10f\n", i, j, weights[i][j]);
		output.write(reinterpret_cast<char *>(&mature_global[i]), sizeof(int));

	}
	output.close();
}

void NetworkGrowthSimulator::write_inhibitory_conductance(int trial, double time_resolution_conductance, std::string outputDir)
{
	std::ofstream output;

	std::string filename = outputDir + "Ginh_trial_" + std::to_string(trial) + ".bin";	
		
	output.open(filename, std::ios::out | std::ios::binary);

	int num_points = static_cast<int>(Ginh_global[0].size());

	output.write(reinterpret_cast<char *>(&N_RA), sizeof(int));
	output.write(reinterpret_cast<char *>(&num_points), sizeof(int));
	output.write(reinterpret_cast<char *>(&time_resolution_conductance), sizeof(double));
		
	for (int i = 0; i < N_RA; i++)
	{
		for (size_t j = 0; j < Ginh_global[i].size(); j++)
			output.write(reinterpret_cast<char *>(&Ginh_global[i][j]), sizeof(double));
	}
	output.close();
}

void NetworkGrowthSimulator::write_activity_history(const char* filename)
{
	std::ofstream output;

	output.open(filename, std::ios::out | std::ios::binary);

	output.write(reinterpret_cast<char *>(&N_RA), sizeof(int));
	output.write(reinterpret_cast<char *>(&maturation_params.RATE_WINDOW_LONG), sizeof(int));
	output.write(reinterpret_cast<char *>(&trial_number), sizeof(int));

	for (int i = 0; i < N_RA; i++)
	{
		for (int j = 0; j < maturation_params.RATE_WINDOW_LONG; j++)
		{
			//printf("weigths[%d][%d] = %1.10f\n", i, j, weights[i][j]);
			output.write(reinterpret_cast<char *>(&num_spikes_in_recent_trials_global[i][j]), sizeof(int));

		}
	}
	output.close();
}

void NetworkGrowthSimulator::write_axonal_delays(const std::vector<std::vector<double>> &axonal_delays, const char* filename)
{
	std::ofstream output;

	output.open(filename, std::ios::out | std::ios::binary);

	int num_neurons = static_cast<int>(axonal_delays.size());

	output.write(reinterpret_cast<char *>(&num_neurons), sizeof(int));
	output.write(reinterpret_cast<char *>(&trial_number), sizeof(int));

	for (int i = 0; i < num_neurons; i++)
	{
		int num_targets = static_cast<int>(axonal_delays[i].size());

		output.write(reinterpret_cast<char *>(&num_targets), sizeof(int));

		for (int j = 0; j < num_targets; j++)
		{
			//printf("weigths[%d][%d] = %1.10f\n", i, j, weights[i][j]);
			output.write(reinterpret_cast<const char *>(&axonal_delays[i][j]), sizeof(double));

		}
	}
	output.close();
}

void NetworkGrowthSimulator::write_weights(const char* filename)
{
	std::ofstream output;

	output.open(filename, std::ios::out | std::ios::binary);

	output.write(reinterpret_cast<char *>(&N_RA), sizeof(int));
	output.write(reinterpret_cast<char *>(&trial_number), sizeof(int));

	for (int i = 0; i < N_RA; i++)
	{
		for (int j = 0; j < N_RA; j++)
		{
			//printf("weigths[%d][%d] = %1.10f\n", i, j, weights[i][j]);
			output.write(reinterpret_cast<char *>(&weights_RA_RA_global[i][j]), sizeof(double));

		}
	}
	output.close();
}

void NetworkGrowthSimulator::write_inhAndExc(int num_trials, 
									const std::vector<int> &num_trials_with_relevant_spikes, 
									const std::vector<double> &mean_relevant_spike_time, 
									const std::vector<double> &std_relevant_spike_time, 
									const std::vector<int> &num_trials_with_nonrelevant_spikes, 
									const std::vector<double> &mean_spike_time,
									const std::vector<double> &std_spike_time,
									const char *filename)
{
	std::ofstream output;

	output.open(filename, std::ios::out | std::ios::binary);

	output.write(reinterpret_cast<char *>(&N_RA), sizeof(int));
	output.write(reinterpret_cast<char *>(&num_trials), sizeof(int));

	for (int i = 0; i < N_RA; i++)
	{
		double probability_for_relevant_spike = static_cast<double>(num_trials_with_relevant_spikes[i]) / 
														static_cast<double>(num_trials);
			
		output.write(reinterpret_cast<char *>(&probability_for_relevant_spike), sizeof(double));
		output.write(reinterpret_cast<const char *>(&mean_relevant_spike_time[i]), sizeof(double));
		output.write(reinterpret_cast<const char *>(&std_relevant_spike_time[i]), sizeof(double));
		
		double probability_for_nonrelevant_spike = static_cast<double>(num_trials_with_nonrelevant_spikes[i]) /
														static_cast<double>(num_trials);
		
		output.write(reinterpret_cast<char *>(&probability_for_nonrelevant_spike), sizeof(double));
		output.write(reinterpret_cast<const char *>(&mean_spike_time[i]), sizeof(double));
		output.write(reinterpret_cast<const char *>(&std_spike_time[i]), sizeof(double));
		
	}
	output.close();	
}


void NetworkGrowthSimulator::write_jitter(int num_trials,
										const std::vector<int> &num_trials_with_soma_spikes, 
										const std::vector<double> &average_num_soma_spikes_in_trial,
										const std::vector<double> &mean_first_soma_spike_time,
										const std::vector<double> &std_first_soma_spike_time,
										const std::vector<int> &num_trials_with_dend_spikes, 
										const std::vector<double> &average_num_dend_spikes_in_trial,
										const std::vector<double> &mean_first_dend_spike_time,
										const std::vector<double> &std_first_dend_spike_time,
										const char *filename)
{
	std::ofstream output;

	output.open(filename, std::ios::out | std::ios::binary);

	output.write(reinterpret_cast<char *>(&N_RA), sizeof(int));
	output.write(reinterpret_cast<char *>(&num_trials), sizeof(int));

	for (int i = 0; i < N_RA; i++)
	{
		double probability_to_soma_spike = static_cast<double>(num_trials_with_soma_spikes[i]) / static_cast<double>(num_trials);
			
		output.write(reinterpret_cast<char *>(&probability_to_soma_spike), sizeof(double));
		output.write(reinterpret_cast<const char *>(&average_num_soma_spikes_in_trial[i]), sizeof(double));
		output.write(reinterpret_cast<const char *>(&mean_first_soma_spike_time[i]), sizeof(double));
		output.write(reinterpret_cast<const char *>(&std_first_soma_spike_time[i]), sizeof(double));
		
		double probability_to_dend_spike = static_cast<double>(num_trials_with_dend_spikes[i]) / static_cast<double>(num_trials);
		
		output.write(reinterpret_cast<char *>(&probability_to_dend_spike), sizeof(double));
		output.write(reinterpret_cast<const char *>(&average_num_dend_spikes_in_trial[i]), sizeof(double));
		output.write(reinterpret_cast<const char *>(&mean_first_dend_spike_time[i]), sizeof(double));
		output.write(reinterpret_cast<const char *>(&std_first_dend_spike_time[i]), sizeof(double));
		
	}
	output.close();
	
	
	
}

void NetworkGrowthSimulator::write_replaced_neurons(const std::vector<int>& replaced_neurons, 
												const char* filename)
{
	if (MPI_rank == 0)
	{
		std::ofstream out;
   
        if (trial_number == 0)
	    	out.open(filename, std::ios::binary | std::ios::out);
        else
	    	out.open(filename, std::ios::binary | std::ios::out | std::ios::app);

		int num_replaced = static_cast<int>(replaced_neurons.size());

        out.write(reinterpret_cast<char*>(&trial_number), sizeof(int)); // write current trial number
		out.write(reinterpret_cast<char*>(&num_replaced), sizeof(int)); // number of replaced neurons
		
		for (int i = 0; i < num_replaced; i++)
			out.write(reinterpret_cast<const char*>(&replaced_neurons[i]), sizeof(int)); // write ids of replaced neurons
	
        out.close();
	}
}



void NetworkGrowthSimulator::write_num_synapses(const char* fileSynapses)
{
	if (MPI_rank == 0)
	{
		std::ofstream out_synapses;
   
        if (trial_number == 0)
	    	out_synapses.open(fileSynapses, std::ios::binary | std::ios::out);
        else
	    	out_synapses.open(fileSynapses, std::ios::binary | std::ios::out | std::ios::app);

		// count active synapses and supersynapses
		int active = 0;
		int super = 0;

		for (int i = 0; i < N_RA; i++)
		{
			active += (int) active_synapses_global[i].size(); 
			super += (int) supersynapses_global[i].size(); 

		}
		
        out_synapses.write(reinterpret_cast<char*>(&trial_number), sizeof(trial_number)); // write current trial numbe
		out_synapses.write(reinterpret_cast<char*>(&active), sizeof(active)); // write number of active synapses
		out_synapses.write(reinterpret_cast<char*>(&super), sizeof(super)); // write number of supersynapses
	
        out_synapses.close();
	}
}


void NetworkGrowthSimulator::write_soma_spike_times(const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;

        out.open(filename, std::ios::out | std::ios::binary );
        out.write(reinterpret_cast<char *>(&trial_number), sizeof(int));
        out.write(reinterpret_cast<const char *>(&TRIAL_DURATION), sizeof(double));
        out.write(reinterpret_cast<char *>(&N_RA), sizeof(int));

        //for (int i = 0; i < N_RA; i++)
        //    printf("spike_times[%d] = %f\n", i, spike_times[i]);
        // write spike times
        for (int i = 0; i < N_RA; i++)
        {
            int spike_array_size = spikes_in_trial_soma_global[i].size();
            //printf("Neuron %d; number of somatic spikes in trial: %d\n", i, spike_array_size);
            out.write(reinterpret_cast<char *>(&spike_array_size), sizeof(int));

            for (int j = 0; j < spike_array_size; j++)
	        {
       		//printf("relative_spike_time_soma = %f\n", relative_spike_time);
                out.write(reinterpret_cast<char *>(&spikes_in_trial_soma_global[i][j]), sizeof(double));
	        }
        }
        //out.write(reinterpret_cast<char *>(spike_times), N_RA*sizeof(double));

        out.close();
    }
}

//~ void NetworkGrowthSimulator::write_last_dend_spike_times(const char* filename)
//~ {
	//~ if (MPI_rank == 0)
    //~ {
        //~ std::ofstream out;
//~ 
        //~ out.open(filename, std::ios::out | std::ios::binary );
        //~ out.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));
        //~ out.write(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));
        //~ 
//~ 
       //~ 
        //~ for (int i = 0; i < N_RA; i++)
        	//~ out.write(reinterpret_cast<char *>(&last_dend_spike_time_global[i]), sizeof(last_dend_spike_time_global[i]));
			//~ 
	//~ 
        //~ out.close();
    //~ }
	//~ 
//~ }

void NetworkGrowthSimulator::write_dend_spike_times(const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;

        out.open(filename, std::ios::out | std::ios::binary );
        out.write(reinterpret_cast<char *>(&trial_number), sizeof(int));
        out.write(reinterpret_cast<const char *>(&TRIAL_DURATION), sizeof(double));
        out.write(reinterpret_cast<char *>(&N_RA), sizeof(int));

        //for (int i = 0; i < N_RA; i++)
        //    printf("spike_times[%d] = %f\n", i, spike_times[i]);
        // write spike times
        for (int i = 0; i < N_RA; i++)
        {
            int spike_array_size = spikes_in_trial_dend_global[i].size();
            //printf("Neuron %d; number of dendritic spikes in trial: %d\n", i, spike_array_size);
            out.write(reinterpret_cast<char *>(&spike_array_size), sizeof(int));
	    
            for (int j = 0; j < spike_array_size; j++)
			{
                //out.write(reinterpret_cast<char *>(&spikes_in_trial_dend_global[i][j]), sizeof(double));
        	out.write(reinterpret_cast<char *>(&spikes_in_trial_dend_global[i][j]), sizeof(double));
			//printf("Neuron %d; relative spike time = %f\n", i, relative_spike_time);
			}
		}
        //out.write(reinterpret_cast<char *>(spike_times), N_RA*sizeof(double));

        out.close();
    }
}

void NetworkGrowthSimulator::write_interneuron_spike_times(const char* filename)
{
    if (MPI_rank == 0)
    {
        std::ofstream out;

        out.open(filename, std::ios::out | std::ios::binary );
        out.write(reinterpret_cast<char *>(&trial_number), sizeof(int));
        out.write(reinterpret_cast<const char *>(&TRIAL_DURATION), sizeof(double));
        out.write(reinterpret_cast<char *>(&N_I), sizeof(int));

        //for (int i = 0; i < N_RA; i++)
        //    printf("spike_times[%d] = %f\n", i, spike_times[i]);
        // write spike times
        for (int i = 0; i < N_I; i++)
        {
            int spike_array_size = spikes_in_trial_interneuron_global[i].size();
            //printf("Neuron %d; number of dendritic spikes in trial: %d\n", i, spike_array_size);
            out.write(reinterpret_cast<char *>(&spike_array_size), sizeof(int));
	    
            for (int j = 0; j < spike_array_size; j++)
	    {
                //out.write(reinterpret_cast<char *>(&spikes_in_trial_dend_global[i][j]), sizeof(double));
        	out.write(reinterpret_cast<char *>(&spikes_in_trial_interneuron_global[i][j]), sizeof(double));
			//printf("Neuron %d; relative spike time = %f\n", i, relative_spike_time);
	    }
	}
        //out.write(reinterpret_cast<char *>(spike_times), N_RA*sizeof(double));

        out.close();
    }
}

void NetworkGrowthSimulator::write_replacement_history(const char* filename)
{
	std::ofstream out;

  
	out.open(filename, std::ios::out | std::ios::binary );
	out.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));
	out.write(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));
	
	
	for (int i = 0; i < N_RA; i++)
		out.write(reinterpret_cast<char *>(&num_trials_after_replacement_global[i]), sizeof(num_trials_after_replacement_global[i]));
		   
	out.close();
}
//~ 
//~ void NetworkGrowthSimulator::write_maturation_info(const char* filename)
//~ {
    //~ if (MPI_rank == 0)
    //~ {
        //~ std::ofstream out;
//~ 
      //~ 
        //~ out.open(filename, std::ios::out | std::ios::binary );
        //~ out.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));
//~ 
//~ 
      //~ 
        //~ out.write(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));
//~ 
        //~ for (int i = 0; i < N_RA; i++)
        //~ {
          //~ 
//~ 
			//~ out.write(reinterpret_cast<char *>(&gaba_potential_global[i]), sizeof(gaba_potential_global[i]));
            //~ out.write(reinterpret_cast<char *>(&firing_rate_short_global[i]), sizeof(firing_rate_short_global[i]));
            //~ out.write(reinterpret_cast<char *>(&firing_rate_long_global[i]), sizeof(firing_rate_long_global[i]));
            //~ out.write(reinterpret_cast<char *>(&remodeled_global[i]), sizeof(remodeled_global[i]));
	        //~ 
	    //~ }
		//~ out.close();
	//~ }
//~ 
//~ }

//~ void NetworkGrowthSimulator::write_maturation_time_sequence(const std::vector<int>& neurons, const char* filename)
//~ {
    //~ if (MPI_rank == 0)
    //~ {
        //~ std::ofstream out;
//~ 
      //~ 
        //~ if (trial_number == 0)
        //~ {
            //~ out.open(filename, std::ios::out | std::ios::binary );
            //~ 
            //~ int size = static_cast<int>(neurons.size());
            //~ 
            //~ out.write(reinterpret_cast<char *>(&size), sizeof(size));
//~ 
            //~ for (size_t i = 0; i < neurons.size(); i++)
                //~ out.write(reinterpret_cast<const char *>(&neurons[i]), sizeof(neurons[i]));
//~ 
        //~ }
        //~ else
            //~ out.open(filename, std::ios::out | std::ios::binary | std::ios::app);
//~ 
        //~ out.write(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));
//~ 
        //~ for (size_t i = 0; i < neurons.size(); i++)
        //~ {
            //~ // find position in global id array
            //~ std::vector<int>::iterator pos = std::find(Id_RA_global.begin(), Id_RA_global.end(), neurons[i]);
//~ 
            //~ if (pos == Id_RA_global.end())
                //~ std::cerr << "Neuronal id " << neurons[i] << " is not found in Id_RA_global in write_maturation_time_sequence" << std::endl; 
            //~ else
            //~ {
                //~ int ind = std::distance(Id_RA_global.begin(), pos);
 //~ 
                //~ out.write(reinterpret_cast<char *>(&gaba_potential_global[ind]), sizeof(gaba_potential_global[ind]));
                //~ out.write(reinterpret_cast<char *>(&firing_rate_global[ind]), sizeof(firing_rate_global[ind]));
                //~ out.write(reinterpret_cast<char *>(&remodeled_global[ind]), sizeof(remodeled_global[ind]));
                //~ out.write(reinterpret_cast<char *>(&mature_global[ind]), sizeof(mature_global[ind]));
	        //~ }
	    //~ 
	    //~ }
		//~ out.close();
	//~ }
//~ }
//~ 
//~ void NetworkGrowthSimulator::write_weights_time_sequence_from_source_to_target(const std::vector<int>& source, const std::vector<int>& target, const char* filename)
//~ {
    //~ if (MPI_rank == 0)
    //~ {
        //~ std::ofstream out;
//~ 
      //~ 
        //~ if (trial_number == 0)
        //~ {
            //~ out.open(filename, std::ios::out | std::ios::binary );
            //~ 
            //~ int size_source = static_cast<int>(source.size());
            //~ int size_target = static_cast<int>(target.size());
            //~ 
            //~ out.write(reinterpret_cast<char *>(&size_source), sizeof(size_source));
            //~ out.write(reinterpret_cast<char *>(&size_target), sizeof(size_target));
//~ 
            //~ for (size_t i = 0; i < source.size(); i++)
                //~ out.write(reinterpret_cast<const char *>(&source[i]), sizeof(source[i]));
//~ 
            //~ for (size_t i = 0; i < target.size(); i++)
                //~ out.write(reinterpret_cast<const char *>(&target[i]), sizeof(target[i]));
        //~ }
        //~ else
            //~ out.open(filename, std::ios::out | std::ios::binary | std::ios::app);
//~ 
        //~ out.write(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));
//~ 
        //~ for (size_t i = 0; i < source.size(); i++)
        //~ {
            //~ for (size_t j = 0; j < target.size(); j++)
            //~ {
                //~ out.write(reinterpret_cast<char *>(&weights_global[source[i]][target[j]]), sizeof(double));
	        //~ }
	    //~ }
		//~ out.close();
	//~ }
//~ 
//~ }
//~ 
//~ void NetworkGrowthSimulator::write_global_index_array(const char* filename)
//~ {
    //~ if (MPI_rank == 0)
    //~ {
        //~ std::ofstream out;
//~ 
      //~ 
        //~ out.open(filename, std::ios::out | std::ios::binary );
//~ 
        //~ out.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));
//~ 
        //~ for (size_t i = 0; i < Id_RA_global.size(); i++)
            //~ out.write(reinterpret_cast<char *>(&Id_RA_global[i]), sizeof(Id_RA_global[i]));
		//~ 
        //~ out.close();
	//~ }
//~ }
//~ 



//~ void NetworkGrowthSimulator::write_replaced_neurons(std::vector<int>& global_id_2replace, const char* filename)
//~ {
    //~ if (MPI_rank == 0)
    //~ {
        //~ std::ofstream out;
//~ 
		//~ if (trial_number == 0)
	    	//~ out.open(filename, std::ios::binary | std::ios::out);
        //~ else
	    	//~ out.open(filename, std::ios::binary | std::ios::out | std::ios::app);
//~ 
		//~ // write trial number
		//~ out.write(reinterpret_cast<char *>(&trial_number), sizeof(trial_number));
//~ 
        //~ int num_replaced = static_cast<int>(global_id_2replace.size()); // number of replaced neurons
//~ 
        //~ // write number of replaced neurons
        //~ out.write(reinterpret_cast<char *>(&num_replaced), sizeof(num_replaced));
//~ 
        //~ for (int i = 0; i < num_replaced; i++)
            //~ out.write(reinterpret_cast<char *>(&global_id_2replace[i]), sizeof(global_id_2replace[i]));
		//~ 
        //~ out.close();
	//~ }
//~ }

//~ void NetworkGrowthSimulator::write_RA(const char* filename, int n)
//~ {
	//~ if (n >= N_RA)
    //~ {
//~ 
        //~ if (MPI_rank == 0)
        //~ {
			//~ printf("Selected neuron ID doesn't exist in the pool! Instead wrote neuron 0 to the file.\n");
            //~ HVCRA_local[0].writeToFile(filename);
    	//~ }
	//~ }
    //~ else
    //~ {
		//~ int N = 0; // number of neurons in all previous processes
        //~ int i = 0;
//~ 
		//~ while (n >= N)
		//~ {
			//~ if (n < N + N_RA_sizes[i])
			//~ {
				//~ if (MPI_rank == i)
				//~ {
					//~ HVCRA_local[n-N].writeToFile(filename);
					//~ //printf("My rank is %d; Writing neuron %d for n-N = %d\n", MPI_rank, n, n-N);
				//~ }
			//~ }
//~ 
			//~ N += N_RA_sizes[i];
			//~ i++;
		//~ }
		//~ 
    //~ }
//~ }

//~ 
//~ void NetworkGrowthSimulator::write_I(const char* filename, int n)
//~ {
    //~ if (n >= N_I)
    //~ {
        //~ if (MPI_rank == 0)
        //~ {
			//~ printf("Selected neuron %d doesn't exist in the pool! Instead wrote neuron 0 to the file.\n",n);
            //~ HVCI_local[0].writeToFile(filename);
    	//~ }
	//~ }
    //~ else
    //~ {
		//~ int N = 0; // number of neurons in all previous processes
        //~ int i = 0;
//~ 
		//~ while (n >= N)
		//~ {
			//~ if (n < N + N_I_sizes[i])
			//~ {
				//~ //printf("My rank is %d; n = %d; N = %d\n", MPI_rank, n, N);
				//~ if (MPI_rank == i)
				//~ {
					//~ //printf("My rank is %d; Writing neuron %d\n", MPI_rank, n);
					//~ HVCI_local[n-N].writeToFile(filename);
//~ 
				//~ }
			//~ }
//~ 
			//~ N += N_I_sizes[i];
			//~ i++;
		//~ }
		//~ 
    //~ }
//~ }

void NetworkGrowthSimulator::write_alterable_network(std::string fileEnding, std::string outputDirectory)
{
	std::string filename_RA_coordinates = outputDirectory + "RA_xy" + fileEnding + ".bin";
	std::string filename_RA_I = outputDirectory + "RA_I_connections" + fileEnding + ".bin";
	std::string filename_I_RA = outputDirectory + "I_RA_connections" + fileEnding + ".bin";
	
	this->write_coordinates(xx_RA, yy_RA, zz_RA, filename_RA_coordinates.c_str());
    
    this->write_connections_RAandI(filename_RA_I.c_str(), filename_I_RA.c_str());
}

void NetworkGrowthSimulator::write_network_topology(std::string outputDirectory)
{
	std::string filename_I_coordinates = outputDirectory + "I_xy.bin";
	std::string filename_RA_coordinates = outputDirectory + "RA_xy.bin";
	
	std::string filename_RA_I = outputDirectory + "RA_I_connections.bin";
	std::string filename_I_RA = outputDirectory + "I_RA_connections.bin";
	
	std::string filename_num_neurons = outputDirectory + "num_neurons.bin";
	
	this->write_coordinates(xx_I,  yy_I,  zz_I,  filename_I_coordinates.c_str());
	this->write_coordinates(xx_RA, yy_RA, zz_RA, filename_RA_coordinates.c_str());
	
	this->write_connections_RAandI(filename_RA_I.c_str(), filename_I_RA.c_str());
	
	this->write_number_of_neurons(filename_num_neurons.c_str());
}

void NetworkGrowthSimulator::write_connections_RAandI(const char* RA_I, const char* I_RA)
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
	for (int i = 0; i < N_RA; i++)
	{
		size = static_cast<int>(syn_ID_RA_I_global[i].size());

		out_RA_I.write(reinterpret_cast<char *>(&i), sizeof(i));
		out_RA_I.write(reinterpret_cast<char *>(&size), sizeof(size)); // write neuron's ID

		for (int j = 0; j < size; j++)
		{
			out_RA_I.write(reinterpret_cast<char *>(&syn_ID_RA_I_global[i][j]), sizeof(int));
			out_RA_I.write(reinterpret_cast<char *>(&weights_RA_I_global[i][j]), sizeof(double));
			out_RA_I.write(reinterpret_cast<char *>(&syn_lengths_RA_I_global[i][j]), sizeof(double)); // write length of connections
			out_RA_I.write(reinterpret_cast<char *>(&axonal_delays_RA_I_global[i][j]), sizeof(double)); // write axonal delays
	
		}

	}

	// write connections from I to RA
	for (int i = 0; i < N_I; i++)
	{
		out_I_RA.write(reinterpret_cast<char *>(&i), sizeof(i)); // write neuron's ID number

		size = static_cast<int>(syn_ID_I_RA_global[i].size());

		out_I_RA.write(reinterpret_cast<char *>(&size), sizeof(size)); // write number of targets a neuron has

		for (int j = 0; j < size; j++)
		{		
			out_I_RA.write(reinterpret_cast<char *>(&syn_ID_I_RA_global[i][j]), sizeof(int)); // write targets ID
			out_I_RA.write(reinterpret_cast<char *>(&weights_I_RA_global[i][j]), sizeof(double)); // write targets conductance
			out_I_RA.write(reinterpret_cast<char *>(&syn_lengths_I_RA_global[i][j]), sizeof(double)); // write length of connections
			out_I_RA.write(reinterpret_cast<char *>(&axonal_delays_I_RA_global[i][j]), sizeof(double)); // write axonal delays
	
		}
	}
	// close files
	out_I_RA.close();
	out_RA_I.close();
}

void NetworkGrowthSimulator::write_pajek_topology(const char* filename)
{
	std::ofstream out;
	out.open(filename, std::ios::out);

	out << "*Vertices " << N_RA + N_I << "\n";
	
	// check dimensionality of a network
	if ( topology_params.arrangement == "square" )
	{
	
		for (int i = 0; i < N_RA; i++)
		{
			// if training neuron, paint with Yellow, otherwise with Green color
			std::vector<int>::iterator iter = std::find(training_neurons.begin(), training_neurons.end(), i);
			
			if (iter != training_neurons.end())
				out << i + 1 << " \"" << i << "\" " << xx_RA[i] << " " << yy_RA[i] << " ic Yellow\n";
			else
				out << i + 1 << " \"" << i << "\" " << xx_RA[i] << " " << yy_RA[i] << " ic Green\n";
		}

		for (int i = 0; i < N_I; i++)
			out << i + N_RA + 1 << " \"" << i + N_RA << "\" " << xx_I[i] << " " << yy_I[i] << " ic Red\n";
		
	}
	if ( topology_params.arrangement == "sphere" )
	{
		
		for (int i = 0; i < N_RA; i++)
		{
			// if training neuron, paint with Yellow, otherwise with Green color
			std::vector<int>::iterator iter = std::find(training_neurons.begin(), training_neurons.end(), i);
			
			if (iter != training_neurons.end())
				out << i + 1 << " \"" << i << "\" " << xx_RA[i] << " " << yy_RA[i] << " " << zz_RA[i] << " ic Yellow\n";
			else
				out << i + 1 << " \"" << i << "\" " << xx_RA[i] << " " << yy_RA[i] << " " << zz_RA[i] << " ic Green\n";
		}

		for (int i = 0; i < N_I; i++)
			out << i + N_RA + 1 << " \"" << i + N_RA << "\" " << xx_I[i] << " " << yy_I[i] << " " << zz_I[i] <<" ic Red\n";
	}
	

	out << "*Arcs\n";
	for (int i = 0; i < N_RA; i++)
	{
		for (size_t j = 0; j < syn_ID_RA_I_global[i].size(); j++)
		{
			int syn_ID = syn_ID_RA_I_global[i][j] + N_RA;
			out << i + 1 << " " << syn_ID + 1 << " " << weights_RA_I_global[i][j] << " c Green\n";
		}
	}
	for (int i = 0; i < N_I; i++)
	{
		for (size_t j = 0; j < syn_ID_I_RA_global[i].size(); j++)
		{
			int syn_ID = syn_ID_I_RA_global[i][j];
			out << i + N_RA + 1 << " " << syn_ID + 1 << " " << weights_I_RA_global[i][j] << " c Red\n";
		}
	}
}

void NetworkGrowthSimulator::write_coordinates(const std::vector<double> &xx, const std::vector<double> &yy,
										 const std::vector<double> &zz, const char* filename)
{
	std::ofstream out;

	// open files
	out.open(filename, std::ios::binary | std::ios::out);

	
	int dimensionality = 3;
	
	if ( xx.empty() )
	{
		std::cerr << "xx array is empty in write coordinates!" << std::endl;
		return;
	}
	
	if ( yy.empty() )
		dimensionality = 1;
		
	if ( zz.empty() )
		dimensionality = 2;
		
	std::cout << "Dimensionality of coordinates in write_coordinates = " << dimensionality << std::endl;

	int num_neurons = static_cast<int>(xx.size());

	// write number of neurons and dimensionality
	out.write(reinterpret_cast<char*>(&num_neurons), sizeof(int));
	out.write(reinterpret_cast<char*>(&dimensionality), sizeof(int));
	out.write(reinterpret_cast<char*>(&MODEL_INTERNEURON_DISTANCE), sizeof(double));
	


	// write coordinates;
	// check dimensionality of a network
	switch (dimensionality)
	{
		case 1:
			for (int i = 0; i < num_neurons; i++)
				out.write(reinterpret_cast<const char*>(&xx[i]), sizeof(double));
			break;
		
		
		case 2:
			for (int i = 0; i < num_neurons; i++)
			{
				out.write(reinterpret_cast<const char*>(&xx[i]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&yy[i]), sizeof(double));
			}
			break;
		
		case 3:
			for (int i = 0; i < num_neurons; i++)
			{
				out.write(reinterpret_cast<const char*>(&xx[i]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&yy[i]), sizeof(double));
				out.write(reinterpret_cast<const char*>(&zz[i]), sizeof(double));
			}
			break;
		
		default:
			std::cerr << "Dimensionality " << dimensionality << " is not supported!" << std::endl;
			break;
	}
	
	// close files
	out.close();
	
}

void NetworkGrowthSimulator::write_training_spread(const char* filename)
{
	std::ofstream out;
  
	out.open(filename, std::ios::out | std::ios::binary );
	
	if (static_cast<int>(training_spread_times.size()) != N_TR)
		std::cerr << "Size of vector with spread_times = " << training_spread_times.size() << " is different from N_TR = " << N_TR << std::endl;
	else
	{
		out.write(reinterpret_cast<char *>(&N_TR), sizeof(int));

		for (int i = 0; i < N_TR; i++)
			out.write(reinterpret_cast<char *>(&training_spread_times[i]), sizeof(double));
	}
		
	out.close();
}

void NetworkGrowthSimulator::write_training_neurons(const char* filename)
{
	std::ofstream out;
  
	out.open(filename, std::ios::out | std::ios::binary );
	
	if (static_cast<int>(training_neurons.size()) != N_TR)
		std::cerr << "Size of vector with training neurons = " << training_neurons.size() << " is different from N_TR = " << N_TR << std::endl;
	else
	{
		out.write(reinterpret_cast<char *>(&N_TR), sizeof(int));

		for (int i = 0; i < N_TR; i++)
			out.write(reinterpret_cast<char *>(&training_neurons[i]), sizeof(training_neurons[i]));
	}
		
	out.close();
}

void NetworkGrowthSimulator::write_number_of_neurons(const char* filename)
{
	std::ofstream out;
  
	out.open(filename, std::ios::out | std::ios::binary );
	
	out.write(reinterpret_cast<char *>(&N_RA), sizeof(int));
	out.write(reinterpret_cast<char *>(&N_I) , sizeof(int));

	out.close();
}

void NetworkGrowthSimulator::write_full_network_state(std::string fileEnding, std::string outputDirectory)
{
	std::string fileWeights = outputDirectory + "weights"  + fileEnding + ".bin";
	
	std::string fileAxonalDelaysRA2RA = outputDirectory + "axonal_delays_RA2RA"  + fileEnding + ".bin";
	std::string fileAxonalDelaysRA2I = outputDirectory + "axonal_delays_RA2I"  + fileEnding + ".bin";
	std::string fileAxonalDelaysI2RA = outputDirectory + "axonal_delays_I2RA"  + fileEnding + ".bin";

	std::string fileWeightStatistics = outputDirectory + "weight_statistics.bin";
	std::string fileActive = outputDirectory + "RA_RA_active_connections" + fileEnding + ".bin";
    std::string fileSuper = outputDirectory + "RA_RA_super_connections" + fileEnding + ".bin";
    std::string fileMaturation = outputDirectory + "mature" + fileEnding + ".bin";
	std::string fileMaturationProperties = outputDirectory + "maturation_properties" + fileEnding + ".bin";
	
	std::string fileActivityHistory = outputDirectory + "activity_history" + fileEnding + ".bin";
    std::string fileReplacementHistory = outputDirectory + "replacement_history" + fileEnding + ".bin";
    std::string fileRemodeledIndicators = outputDirectory + "remodeled_indicators" + fileEnding + ".bin";
    
    std::string fileTimeSoma = outputDirectory + "spike_times_soma" + fileEnding + ".bin";
    std::string fileTimeDend = outputDirectory + "spike_times_dend" + fileEnding + ".bin";
    std::string fileTimeInterneuron = outputDirectory + "spike_times_interneuron" + fileEnding + ".bin";
    
    this->write_alterable_network(fileEnding, outputDirectory);
    this->write_weights(fileWeights.c_str());
    this->write_weight_statistics(fileWeightStatistics.c_str());
    
    this->write_axonal_delays(axonal_delays_RA_RA_global, fileAxonalDelaysRA2RA.c_str());
	this->write_axonal_delays(axonal_delays_RA_I_global, fileAxonalDelaysRA2I.c_str());
	this->write_axonal_delays(axonal_delays_I_RA_global, fileAxonalDelaysI2RA.c_str());
				
    this->write_super_synapses(fileSuper.c_str());
    this->write_active_synapses(fileActive.c_str());
    this->write_mature_indicators(fileMaturation.c_str());
    this->write_maturation_properties(fileMaturationProperties.c_str());
    this->write_activity_history(fileActivityHistory.c_str());
    this->write_replacement_history(fileReplacementHistory.c_str());
    this->write_remodeled_indicators(fileRemodeledIndicators.c_str());
    
    this->write_soma_spike_times(fileTimeSoma.c_str());
	this->write_dend_spike_times(fileTimeDend.c_str());
	this->write_interneuron_spike_times(fileTimeInterneuron.c_str());
}

void NetworkGrowthSimulator::write_chain_test(std::string outputDirectory, int trial_number)
{
	std::string fileTimeSoma = outputDirectory + "test_spike_times_soma_" + std::to_string(trial_number) + ".bin";
    std::string fileTimeDend = outputDirectory + "test_spike_times_dend_" + std::to_string(trial_number) + ".bin";
    std::string fileTimeInterneuron = outputDirectory + "test_spike_times_interneuron_" + std::to_string(trial_number) + ".bin";
    
	this->write_soma_spike_times(fileTimeSoma.c_str());
	this->write_dend_spike_times(fileTimeDend.c_str());
	this->write_interneuron_spike_times(fileTimeInterneuron.c_str());
}

void NetworkGrowthSimulator::write_inhibition_tracking_state(int trial_number, double time_resolution_conductance, std::string outputDirectory)
{
	std::string fileNumSynapses = outputDirectory + "num_synapses_" + std::to_string(trial_number) + ".bin";
	std::string fileTimeSoma = outputDirectory + "spike_times_soma_" + std::to_string(trial_number) + ".bin";
    std::string fileTimeDend = outputDirectory + "spike_times_dend_" + std::to_string(trial_number) + ".bin";
    std::string fileTimeInterneuron = outputDirectory + "spike_times_interneuron_" + std::to_string(trial_number) + ".bin";
    
	this->write_num_synapses(fileNumSynapses.c_str());
	this->write_soma_spike_times(fileTimeSoma.c_str());
	this->write_dend_spike_times(fileTimeDend.c_str());
	this->write_interneuron_spike_times(fileTimeInterneuron.c_str());
	this->write_inhibitory_conductance(trial_number, time_resolution_conductance, outputDirectory);
}

void NetworkGrowthSimulator::write_graph_network_state(std::string outputDirectory)
{
	std::string fileNumSynapses = outputDirectory + "num_synapses.bin";
	std::string fileWeightStatistics = outputDirectory + "weight_statistics.bin";
	
    std::string fileActive = outputDirectory + "RA_RA_active_connections.bin";
    std::string fileSuper = outputDirectory + "RA_RA_super_connections.bin";
    std::string fileTimeSoma = outputDirectory + "spike_times_soma.bin";
    std::string fileTimeDend = outputDirectory + "spike_times_dend.bin";
    std::string fileTimeInterneuron = outputDirectory + "spike_times_interneuron.bin";
    

	this->write_num_synapses(fileNumSynapses.c_str());
	this->write_soma_spike_times(fileTimeSoma.c_str());
	this->write_dend_spike_times(fileTimeDend.c_str());
	this->write_interneuron_spike_times(fileTimeInterneuron.c_str());
   
	this->write_super_synapses(fileSuper.c_str());
	this->write_active_synapses(fileActive.c_str());
}
