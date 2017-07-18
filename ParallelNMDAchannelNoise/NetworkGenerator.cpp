#include "NetworkGenerator.h"
#include <algorithm>
#include <iostream>
#include <fstream>

NetworkGenerator::NetworkGenerator(ConfigurationNetworkGenerator cfg, Poisson_noise* generator) : config(cfg)
{
    struct SpatialParameters spatial_params = cfg.get_spatial_parameters();
    
    this->set_spatial_parameters(spatial_params);
 
    noise_generator = generator;  
}

void NetworkGenerator::set_spatial_parameters(const struct SpatialParameters& spatial_params)
{
	dimensionality = spatial_params.dimensionality;
	
    MIN_INTERNEURON_DISTANCE = spatial_params.MIN_INTERNEURON_DISTANCE;

    A_RA2I = spatial_params.A_RA2I;
    SIGMA_RA2I = spatial_params.SIGMA_RA2I;
    B_I2RA = spatial_params.B_I2RA;
    SIGMA_I2RA = spatial_params.SIGMA_I2RA;

    Gei_mean = spatial_params.Gei_mean;
    Gei_std = spatial_params.Gei_std;
    Gie_mean = -1.0;
    Gie_std = 0.0;
}

void NetworkGenerator::set_inhibitory_strength(double G_mean, double G_std)
{
	Gie_mean = G_mean;
	Gie_std = G_std;
}

void NetworkGenerator::get_neuron_numbers(int* N_ra, int* N_tr, int* N_i)
{
	*N_ra = N_RA;
	*N_tr = N_TR;
	*N_i = N_I;
}

void NetworkGenerator::get_axonal_delays(double delay_constant, std::vector<std::vector<double>>& axonal_delays_RA_I, 
			 std::vector<std::vector<double>>& axonal_delays_RA_RA)
{
	// check dimensionality
	switch (dimensionality)
	{
		case 2:		
			for (int i = 0; i < N_RA; i++)
			{
				for (int j = 0; j < N_I; j++)
					axonal_delays_RA_I[i][j] = delay_constant * distance(xx_RA[i], yy_RA[i], xx_I[j], yy_I[j]);
			
				for (int j = 0; j < N_RA; j++)
					axonal_delays_RA_RA[i][j] = delay_constant * distance(xx_RA[i], yy_RA[i], xx_RA[j], yy_RA[j]);
			}
			break;
		
		case 3:		
			for (int i = 0; i < N_RA; i++)
			{
				for (int j = 0; j < N_I; j++)
					axonal_delays_RA_I[i][j] = delay_constant * distance_on_sphere(xx_RA[i], yy_RA[i], zz_RA[i], xx_I[j], yy_I[j], zz_I[j]);
			
				for (int j = 0; j < N_RA; j++)
					axonal_delays_RA_RA[i][j] = delay_constant * distance_on_sphere(xx_RA[i], yy_RA[i], zz_RA[i], xx_RA[j], yy_RA[j], zz_RA[j]);
			}
			break;
			
		default:
			std::cerr << "Dimensionality " << dimensionality << " is not supported for setting up connections" << std::endl;
			break;		
	}
}

void NetworkGenerator::get_network(std::vector<std::vector<int>>* id_RA_I, std::vector<std::vector<double>>* w_RA_I, 
					std::vector<std::vector<int>>* id_I_RA, std::vector<std::vector<double>>* w_I_RA, std::vector<int>* training)
{	
	*id_RA_I = syn_ID_RA_I;
	*w_RA_I = weights_RA_I;
	*id_I_RA = syn_ID_I_RA;
	*w_I_RA = weights_I_RA;
	*training = training_neurons;
}

void NetworkGenerator::generate_default_network(int N_ra, int N_tr, int N_i, std::string outputDirectory)
{
	N_RA = N_ra;
	N_TR = N_tr;
	N_I = N_i;
	
	// resize arrays
    // coordinates
    switch (dimensionality)
    {
		case 1:
			xx_RA.resize(N_RA);
			xx_I.resize(N_I);
			break;
			
		case 2:
			xx_RA.resize(N_RA);
			yy_RA.resize(N_RA);
			xx_I.resize(N_I);
			yy_I.resize(N_I);
			break;
		
		case 3:
			xx_RA.resize(N_RA);
			yy_RA.resize(N_RA);
			zz_RA.resize(N_RA);
			xx_I.resize(N_I);
			yy_I.resize(N_I);
			zz_I.resize(N_I);
			break;
	}
    
    // connections
    weights_RA_I.resize(N_RA);
	syn_ID_RA_I.resize(N_RA);

	weights_I_RA.resize(N_I);
	syn_ID_I_RA.resize(N_I);
	
	this->initialize_coordinates();
	this->initialize_connections();	
	
	this->write_invariant_network_to_directory(outputDirectory);
	this->write_alterable_network_to_directory("_initial", outputDirectory);
}

void NetworkGenerator::generate_network_from_data(int N_ra, int N_tr, int N_i)
{
	N_RA = N_ra;
	N_TR = N_tr;
	N_I = N_i;
	
	// resize arrays
    // coordinates
    xx_RA.resize(N_RA);
    yy_RA.resize(N_RA);
    
    xx_I.resize(N_I);
    yy_I.resize(N_I);
    
    // connections
    weights_RA_I.resize(N_RA);
	syn_ID_RA_I.resize(N_RA);

	weights_I_RA.resize(N_I);
	syn_ID_I_RA.resize(N_I);
	
	weights_RA_RA.resize(N_RA);
	
	for (int i = 0; i < N_RA; i++)
		weights_RA_RA[i].resize(N_RA);
	
	this->initialize_coordinates();
	this->initialize_connections();
	this->initialize_connections_between_HVCRA();
	 
}

//~ void NetworkGenerator::generate_network_with_clustered_training(int N_ra, int N_tr, int N_i)
//~ {
	//~ this->generate_default_network(N_ra, N_tr, N_i);
	//~ 
	//~ // find neurons that are close to neuron 0 to from a training set cluster
	//~ // select neuron 0 as the first training neuron and find distances to all other neurons in the pool
	//~ double xx = xx_RA[0];
	//~ double yy = yy_RA[0];
//~ 
	//~ std::vector<double> distances_to_pool_neurons(xx_RA.size()-1); 
//~ 
	//~ for (size_t i = 1; i < xx_RA.size(); i++)
		//~ distances_to_pool_neurons[i-1] = distance(xx, yy, xx_RA[i], yy_RA[i]);
//~ 
	//~ // sort distances
	//~ std::vector<size_t> idx(distances_to_pool_neurons.size()); // vector with sorted indices of distances
//~ 
	//~ std::iota(idx.begin(), idx.end(), 0);
//~ 
	//~ std::sort(idx.begin(), idx.end(), [&distances_to_pool_neurons](size_t i1, size_t i2)
										//~ {return distances_to_pool_neurons[i1] < distances_to_pool_neurons[i2];});
//~ 
	//~ // select neuron 0 and N_TR - 1 closest neurons to be training set
	//~ 
	//~ training_neurons[0] = 0;
	//~ 
	//~ for (int i = 1; i < N_TR; i++) 
		//~ training_neurons[i] = idx[i-1] + 1;
//~ }
//~ 
//~ void NetworkGenerator::generate_network_with_dispersed_training(int N_ra, int N_tr, int N_i)
//~ {
	//~ this->generate_default_network(N_ra, N_tr, N_i);
	//~ 
	//~ // fix training neuron 0 and move other training neurons if necessary
	//~ training_neurons[0] = 0;
	//~ 
	//~ int neuron_id_to_check = 1; // id of the next neuron to check
	//~ bool found_next_distant_neuron;
//~ 
	//~ for (int i = 1; i < N_TR; i++)
	//~ {
		//~ 
		//~ // calculate distances from next_neuron_to_check to previously confirmed distant training neurons
		//~ do 
		//~ {
			//~ found_next_distant_neuron = true;
			//~ 
			//~ for (int j = 0; j < i; j++)
				//~ if (distance(xx_RA[j], yy_RA[j], xx_RA[neuron_id_to_check], yy_RA[neuron_id_to_check]) < 5 / (sqrt(N_I)+1))
				//~ {
					//~ 
					//~ std::cout << "Neuron " << neuron_id_to_check << std::endl;
						//~ 
					//~ std::cout << "Distance to neuron " << j << " = " << distance(xx_RA[j], yy_RA[j], xx_RA[neuron_id_to_check], yy_RA[neuron_id_to_check]) 
							//~ << " is smaller than " << 5 / (sqrt(N_I) + 1) << std::endl;
//~ 
					//~ neuron_id_to_check++;
					//~ found_next_distant_neuron = false;
					//~ break;
				//~ }
		//~ } while(!found_next_distant_neuron);
		//~ 
		//~ training_neurons[i] = neuron_id_to_check;
		//~ 
		//~ neuron_id_to_check++;	
	//~ }
//~ }

void NetworkGenerator::replace_neurons(const std::vector<int>& neurons_to_replace, std::string extension, std::string outputDirectory)
{
	this->initialize_coordinates_for_replaced_neurons(neurons_to_replace);
	this->initialize_connections_for_replaced_neurons(neurons_to_replace);
	
	this->write_alterable_network_to_directory(extension, outputDirectory);
}

void NetworkGenerator::initialize_coordinates()
{   
	// check network dimensionality
	switch (dimensionality)
	{
		case 2:
		{
			double xx; // temporary x-coordinate
			double yy; // temporary y-coordinate

			bool close; // are neurons too close or not
			// set coordinates for HVC(I) neurons
			for (int i = 0; i < (int) sqrt(N_I); i++)
			{
				for (int k = 0; k < (int) sqrt(N_I); k++)
				{
					do
					{
						close = false;
						xx = (double) (i+1) / (sqrt(N_I)+1) + noise_generator->random(0.25 / (sqrt(N_I) + 1));
						yy = (double) (k+1) / (sqrt(N_I)+1) + noise_generator->random(0.25 / (sqrt(N_I) + 1));

						// check distances to all previous I neurons

						for (size_t j = 0; j < i * ((int) sqrt(N_I)) + k; j++)
						{
							if (distance(xx,yy,xx_I[j],yy_I[j]) < MIN_INTERNEURON_DISTANCE)
							{
								close = true;
								break;
							}
						}
					} while(close);

				xx_I[i * ((int) sqrt(N_I)) + k] = xx;
				yy_I[i * ((int) sqrt(N_I)) + k] = yy;

				}
			}

			// set coordinates for HVC(RA) neurons
			for (int i = 0; i < N_RA; i++)
			{
				do
				{
					close = false;
					xx = 0.5 / (sqrt(N_I)+1) + noise_generator->random(1.0 - 1.0 / (sqrt(N_I)+1));
					yy = 0.5 / (sqrt(N_I)+1) + noise_generator->random(1.0 - 1.0 / (sqrt(N_I)+1));

					// check distances to all I neurons

					for (int j = 0; j < N_I; j++)
					{
						if (distance(xx,yy,xx_I[j],yy_I[j]) < MIN_INTERNEURON_DISTANCE)
						{
							close = true;
							break;
						}
					}

					// check distances to all previous RA neurons
					if (!close)
						for (int j = 0; j < i; j++)
						{
							if (distance(xx, yy, xx_RA[j], yy_RA[j]) < MIN_INTERNEURON_DISTANCE)
						{
							close = true;
							break;
						}

					}
				} while(close);

				xx_RA[i] = xx;
				yy_RA[i] = yy;
			}
			break;
		}
		
		case 3:
		{
			double zz; // temporary z-coordinate

			// set coordinates for HVC(I) neurons
			// HVC(I) neurons are sampled equally distantly on sphere
			const double pi = 3.14159265358979323846; // value of constant pi
			
			double dphi = pi * (3.0 - sqrt(5.0)); // angle phi increment value
			double dz = 2 / static_cast<double>(N_I); // z increment value
			
			double phi = 0; // current value of angle phi
			double z = 1 - dz/2.0; // current value of z  
			double r; // current value of disk radius

			for (int i = 0; i < N_I; i++)
			{
				r = sqrt(1.0 - z*z);
				
				xx_I[i] = cos(phi) * r;
				yy_I[i] = sin(phi) * r;
				zz_I[i] = z;
				
				z = z - dz;
				phi = phi + dphi;
			}

			// set coordinates for HVC(RA) neurons
			// HVC(RA) neurons are sample from uniform distribution on sphere
			for (int i = 0; i < N_RA; i++)
			{
				z = 2*noise_generator->random(1.0) - 1; // generate z in range(-1, 1)
				phi = noise_generator->random(2.0 * pi); // generate phi in range(0, 2*pi)
				
				xx_RA[i] = sin(acos(z)) * cos(phi);
				yy_RA[i] = sin(acos(z)) * sin(phi);
				zz_RA[i] = z;		
			}
			break;
		}
		
		default:
			std::cerr << "Dimensionality " << dimensionality << " is not supported for setting up connections" << std::endl;
			break;
	}
}

void NetworkGenerator::initialize_connections()
{
	switch (dimensionality)
	{
		case 2:	
		{	
			// connections for HVC(RA) neurons
			for (int i = 0; i < N_RA; i++)
			{
				for (int j = 0; j < N_I; j++)
				{
					double d = distance(xx_RA[i], yy_RA[i], xx_I[j], yy_I[j]); // distance between HVC(RA) and HVC(I)
				
					if (noise_generator->random(1) < p_RA2I(d))
					{
						double G = this->sample_Ge2i();

						weights_RA_I[i].push_back(G);
						syn_ID_RA_I[i].push_back(j);
					}
				}

			 }
			// connections for HVC(I) neurons

			for (int i = 0; i < N_I; i++)
			{
				for (int j = 0; j < N_RA; j++)
				{
					double d = distance(xx_I[i], yy_I[i], xx_RA[j], yy_RA[j]); // distance between HVC(I) and HVC(RA)
				
					if (noise_generator->random(1) < p_I2RA(d))
					{
						double G = this->sample_Gi2e();

						weights_I_RA[i].push_back(G);
						syn_ID_I_RA[i].push_back(j);
					}
				}
			}
			break;
		}
		
		case 3:
		{
			// connections for HVC(RA) neurons
			for (int i = 0; i < N_RA; i++)
			{
				for (int j = 0; j < N_I; j++)
				{
					double d = distance_on_sphere(xx_RA[i], yy_RA[i], zz_RA[i], xx_I[j], yy_I[j], zz_I[j]); // distance between HVC(RA) and HVC(I)
				
					if (noise_generator->random(1) < p_RA2I(d))
					{
						double G = this->sample_Ge2i();

						weights_RA_I[i].push_back(G);
						syn_ID_RA_I[i].push_back(j);
					}
				}

			 }
			// connections for HVC(I) neurons

			for (int i = 0; i < N_I; i++)
			{
				for (int j = 0; j < N_RA; j++)
				{
					double d = distance_on_sphere(xx_I[i], yy_I[i], zz_I[i], xx_RA[j], yy_RA[j], zz_RA[j]); // distance between HVC(I) and HVC(RA)
				
					if (noise_generator->random(1) < p_I2RA(d))
					{
						double G = this->sample_Gi2e();

						weights_I_RA[i].push_back(G);
						syn_ID_I_RA[i].push_back(j);
					}
				}
			}
			break;
		}
		
		default:
			std::cerr << "Dimensionality " << dimensionality << " is not supported for setting up connections" << std::endl;
			break;
			
	}
}

void NetworkGenerator::initialize_connections_between_HVCRA()
{
	for (int i = 0; i < N_RA; i++)
		for (int j = 0; j < N_RA; j++)
		{
			if ( (i != 0) && (j % i == 0) )
				weights_RA_RA[i][j] = 0.5;
			else
				weights_RA_RA[i][j] = 0.0;
		}
}

void NetworkGenerator::write_network_from_data_to_directory(double active_threshold, double super_threshold, std::string outputDirectory)
{
	
	this->write_invariant_network_to_directory(outputDirectory);
	this->write_alterable_network_to_directory("_0_", outputDirectory);
	
	// write active synapses
	std::string file_active = outputDirectory + "RA_RA_active_connections_0_.bin";
	std::string file_super = outputDirectory + "RA_RA_super_connections_0_.bin";
	
	std::ofstream out_active, out_super;
	
	out_active.open(file_active.c_str(), std::ios::out | std::ios::binary);
	out_super.open(file_super.c_str(), std::ios::out | std::ios::binary);

	// write number of neurons to file
	out_active.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));
	out_super.write(reinterpret_cast<char *>(&N_RA), sizeof(N_RA));

	int num_active; // number of active outputs for each neuron
	int num_super; // number of super outputs for each neuron
	
	std::vector<int> id_active; // id of target neurons with active synapses
	std::vector<int> id_super; // id of target neurons with super synapses
	
	// calculate number of active and super output connections for each neuron
	// and write them to files
	for (int i = 0; i < N_RA; i++)
	{
		num_active = 0;
		num_super = 0;
		
		id_active.clear();
		id_super.clear();
		
		for (int j = 0; j < N_RA; j++)
		{
			if (weights_RA_RA[i][j] >= active_threshold)
			{
				num_active++;
				id_active.push_back(j);
			}
			
			if (weights_RA_RA[i][j] >= super_threshold)
			{
				num_super++;
				id_super.push_back(j);
			}
		}
		
		out_active.write(reinterpret_cast<char *>(&i), sizeof(i));
		out_active.write(reinterpret_cast<char *>(&num_active), sizeof(num_active));

		out_super.write(reinterpret_cast<char *>(&i), sizeof(i));
		out_super.write(reinterpret_cast<char *>(&num_super), sizeof(num_super));


		for (int k = 0; k < num_active; k++)
		{
			out_active.write(reinterpret_cast<char *>(&id_active[k]), sizeof(id_active[k]));
			out_active.write(reinterpret_cast<char *>(&weights_RA_RA[i][id_active[k]]), sizeof(weights_RA_RA[i][id_active[k]]));
		}
		
		for (int k = 0; k < num_super; k++)
		{
			out_super.write(reinterpret_cast<char *>(&id_super[k]), sizeof(id_super[k]));
			out_super.write(reinterpret_cast<char *>(&weights_RA_RA[i][id_super[k]]), sizeof(weights_RA_RA[i][id_super[k]]));
		}
	
	}
	
	out_active.close();
	out_super.close();
}


void NetworkGenerator::initialize_connections_for_replaced_neurons(const std::vector<int>& neurons_to_replace)
{
	// first clear all existing connections from and to replaced neurons
	for (size_t i = 0; i < neurons_to_replace.size(); i++)
	{
		weights_RA_I[neurons_to_replace[i]].clear();
		syn_ID_RA_I[neurons_to_replace[i]].clear();
		
		
		// erase all connections from interneurons to replaced neuron
        for (int j = 0; j < N_I; j++)
        {
            std::vector<int>::iterator pos = std::find(syn_ID_I_RA[j].begin(), syn_ID_I_RA[j].end(), neurons_to_replace[i]);

            // if interneuron targets replaced neuron, erase connection
            if (pos != syn_ID_I_RA[j].end())
            {
                int index_in_array = std::distance(syn_ID_I_RA[j].begin(), pos);

                syn_ID_I_RA[j].erase(pos);
                weights_I_RA[j].erase(weights_I_RA[j].begin() + index_in_array);
            }

        }	
		
	}

		
	// generate connections for replaced HVC(RA) neurons
	// HVC(RA) -> HVC(I) connections
	
	// first check network dimensionality
	switch (dimensionality)
	{
		case 2:				
			for (size_t i = 0; i < neurons_to_replace.size(); i++)
			{
				for (int j = 0; j < N_I; j++)
				{
					double d = distance(xx_RA[neurons_to_replace[i]], yy_RA[neurons_to_replace[i]], xx_I[j], yy_I[j]); // distance between replaced neuron and HVC(I)
					 
					 if (noise_generator->random(1) < p_RA2I(d))
					 {
						 double G = this->sample_Ge2i();

						 weights_RA_I[neurons_to_replace[i]].push_back(G);
						 syn_ID_RA_I[neurons_to_replace[i]].push_back(j);
					 }
				 }

			}
			
			// HVC(I) -> HVC(RA) connections
			for (int i = 0; i < N_I; i++)
			{
				for (size_t j = 0; j < neurons_to_replace.size(); j++)
				{
					double d = distance(xx_I[i], yy_I[i], xx_RA[neurons_to_replace[j]], yy_RA[neurons_to_replace[j]]); // distance between HVC(I) and replaced neuron
					
					 if (noise_generator->random(1) < p_I2RA(d))
					 {
						 double G = this->sample_Gi2e();

						 weights_I_RA[i].push_back(G);
						 syn_ID_I_RA[i].push_back(neurons_to_replace[j]);
					 }
				 }
			 }
			 
			 break;
			 
		case 3:
			for (size_t i = 0; i < neurons_to_replace.size(); i++)
			{
				for (int j = 0; j < N_I; j++)
				{
					double d = distance_on_sphere(xx_RA[neurons_to_replace[i]], yy_RA[neurons_to_replace[i]], zz_RA[neurons_to_replace[i]], xx_I[j], yy_I[j], zz_I[j]); // distance between replaced neuron and HVC(I)
					 
					 if (noise_generator->random(1) < p_RA2I(d))
					 {
						 double G = this->sample_Ge2i();

						 weights_RA_I[neurons_to_replace[i]].push_back(G);
						 syn_ID_RA_I[neurons_to_replace[i]].push_back(j);
					 }
				 }

			}
			
			// HVC(I) -> HVC(RA) connections
			for (int i = 0; i < N_I; i++)
			{
				for (size_t j = 0; j < neurons_to_replace.size(); j++)
				{
					double d = distance_on_sphere(xx_I[i], yy_I[i], zz_I[i], xx_RA[neurons_to_replace[j]], yy_RA[neurons_to_replace[j]], zz_RA[neurons_to_replace[j]]); // distance between HVC(I) and replaced neuron
					
					 if (noise_generator->random(1) < p_I2RA(d))
					 {
						 double G = this->sample_Gi2e();

						 weights_I_RA[i].push_back(G);
						 syn_ID_I_RA[i].push_back(neurons_to_replace[j]);
					 }
				 }
			 }
			break;
		
		default:
			std::cerr << "Dimensionality " << dimensionality << " is not supported for setting up connections for replaced neurons" << std::endl;
			break;
	 
	}
}


void NetworkGenerator::initialize_coordinates_for_replaced_neurons(const std::vector<int>& neurons_to_replace)
{
	switch (dimensionality)
	{
		case 2:
		{
			std::vector<double> new_xx(neurons_to_replace.size()); // new x-coordinates of replaced neurons
			std::vector<double> new_yy(neurons_to_replace.size()); // new y-coordinates of replaced neurons
			
			double xx; // temporary x-coordinate
			double yy; // temporary y-coordinate

			bool close; // are neurons too close or not

			
			// set coordinates for replaced HVC(RA) neurons
			for (size_t i = 0; i < neurons_to_replace.size(); i++)
			{
				do
				{
					close = false;
					xx = 0.5 / (sqrt(N_I)+1) + noise_generator->random(1.0 - 1.0 / (sqrt(N_I)+1));
					yy = 0.5 / (sqrt(N_I)+1) + noise_generator->random(1.0 - 1.0 / (sqrt(N_I)+1));

					// check distances to all HVC(I) neurons

					for (int j = 0; j < N_I; j++)
					{
						if (distance(xx,yy,xx_I[j],yy_I[j]) < MIN_INTERNEURON_DISTANCE)
						{
							close = true;
							break;
						}
					}

					// check distances to all HVC(RA) neurons that we had before
					if (!close)
						for (int j = 0; j < N_RA; j++)
						{
							if (distance(xx, yy, xx_RA[j], yy_RA[j]) < MIN_INTERNEURON_DISTANCE)
						{
							close = true;
							break;
						}

					}
				} while(close);

				new_xx[i] = xx;
				new_yy[i] = yy;
			}
			
			// assign new coordinates to replaced neurons:
			for (size_t i = 0; i < neurons_to_replace.size(); i++)
			{
				xx_RA[neurons_to_replace[i]] = new_xx[i];
				yy_RA[neurons_to_replace[i]] = new_yy[i];
			}
			
			break;
		}
		
		case 3:
		{
			std::vector<double> new_xx(neurons_to_replace.size()); // new x-coordinates of replaced neurons
			std::vector<double> new_yy(neurons_to_replace.size()); // new y-coordinates of replaced neurons
			std::vector<double> new_zz(neurons_to_replace.size()); // new z-coordinates of replaced neurons
		
			// set coordinates for HVC(RA) neurons
			// HVC(RA) neurons are sample from uniform distribution on sphere
			const double pi = 3.14159265358979323846; // value of constant pi
			
			double z; // z-coord
			double phi; // angle phi
			
			
			for (size_t i = 0; i < neurons_to_replace.size(); i++)
			{
				z = 2*noise_generator->random(1.0) - 1; // generate z in range(-1, 1)
				phi = noise_generator->random(2.0 * pi); // generate phi in range(0, 2*pi)
				
				new_xx[i] = sin(acos(z)) * cos(phi);
				new_yy[i] = sin(acos(z)) * sin(phi);
				new_zz[i] = z;		
			}
		
			// assign new coordinates to replaced neurons:
			for (size_t i = 0; i < neurons_to_replace.size(); i++)
			{
				xx_RA[neurons_to_replace[i]] = new_xx[i];
				yy_RA[neurons_to_replace[i]] = new_yy[i];
				zz_RA[neurons_to_replace[i]] = new_zz[i];
			}
			
		
			break;
		}
		
		default:
			std::cerr << "Dimensionality " << dimensionality << " is not supported for setting up coordinates for replaced neurons" << std::endl;
			break;
	}
}

void NetworkGenerator::sample_weights_based_on_connections()
{
	for (int i = 0; i < N_RA; i++)
		for (size_t j = 0; j < syn_ID_RA_I[i].size(); j++)
			weights_RA_I[i][j] = this->sample_Ge2i();
	
	for (int i = 0; i < N_I; i++)
		for (size_t j = 0; j < syn_ID_I_RA[i].size(); j++)
			weights_I_RA[i][j] = this->sample_Gi2e();
	
}

double NetworkGenerator::sample_Ge2i()
{
    return Gei_mean + noise_generator->normal_distribution() * Gei_std;
}

double NetworkGenerator::sample_Gi2e()
{
    return Gie_mean + noise_generator->normal_distribution() * Gie_std;
}

double NetworkGenerator::p_RA2I(double d)
{
	return A_RA2I * exp(-d * d / (2 * SIGMA_RA2I * SIGMA_RA2I));
}

double NetworkGenerator::p_I2RA(double d)
{
	return B_I2RA * exp(-d*d / (2 * SIGMA_I2RA * SIGMA_I2RA));
}


void NetworkGenerator::read_network_from_directory(std::string extension, std::string networkDir, std::string fileTraining)
{
	this->read_network_with_weights_from_directory(extension, networkDir, fileTraining);
	this->sample_weights_based_on_connections();
}

void NetworkGenerator::read_network_with_weights_from_directory(std::string extension, std::string networkDir, std::string fileTraining)
{
	std::string fileRA2I = networkDir + "RA_I_connections" + extension + ".bin";
    std::string fileI2RA = networkDir + "I_RA_connections" + extension + ".bin";

	std::string filename_RA_xy = networkDir + "RA_xy" + extension + ".bin";
	std::string filename_I_xy = networkDir + "I_xy.bin";
   
    this->read_coordinates(filename_RA_xy.c_str(), filename_I_xy.c_str()); 
    this->read_connections_and_weights(fileRA2I.c_str(), fileI2RA.c_str());
	this->read_training_neurons(fileTraining.c_str()); 
    
}


void NetworkGenerator::read_coordinates(const char* RA_xy, const char* I_xy)
{
	std::ifstream inp_RA, inp_I;

	// open files
	inp_RA.open(RA_xy, std::ios::binary | std::ios::in);
	inp_I.open(I_xy, std::ios::binary | std::ios::in);

	// read number of neurons
	
	inp_RA.read(reinterpret_cast<char*>(&N_RA), sizeof(N_RA));
	inp_I.read(reinterpret_cast<char*>(&N_I), sizeof(N_I));
	
	//std::cout << "Number of neurons in file with HVC(RA) coordinates N_RA = " << N_RA << std::endl;
	//std::cout << "Number of neurons in file with HVC(I) coordinates N_I = " << N_I << std::endl;
	std::cout << "Dimensionality in read_coordinates: " << dimensionality << std::endl;
	// check dimensionality of a network
	switch (dimensionality)
	{
		case 2:
			// resize arrays
			xx_RA.resize(N_RA);
			yy_RA.resize(N_RA);
			
			xx_I.resize(N_I);
			yy_I.resize(N_I);
			
			// read RA coordinates
			for (int i = 0; i < N_RA; i++)
			{
				inp_RA.read(reinterpret_cast<char*>(&xx_RA[i]), sizeof(xx_RA[i]));
				inp_RA.read(reinterpret_cast<char*>(&yy_RA[i]), sizeof(yy_RA[i]));
			}

			// read I coordinates
			for (int i = 0; i < N_I; i++)
			{
				inp_I.read(reinterpret_cast<char*>(&xx_I[i]), sizeof(xx_I[i]));
				inp_I.read(reinterpret_cast<char*>(&yy_I[i]), sizeof(yy_I[i]));
			}

			break;
			
		case 3:
			// resize arrays
			xx_RA.resize(N_RA);
			yy_RA.resize(N_RA);
			zz_RA.resize(N_RA);
			
			xx_I.resize(N_I);
			yy_I.resize(N_I);
			zz_I.resize(N_I);

			for (int i = 0; i < N_RA; i++)
			{
				inp_RA.read(reinterpret_cast<char*>(&xx_RA[i]), sizeof(xx_RA[i]));
				inp_RA.read(reinterpret_cast<char*>(&yy_RA[i]), sizeof(yy_RA[i]));
				inp_RA.read(reinterpret_cast<char*>(&zz_RA[i]), sizeof(zz_RA[i]));
			}

			// read I coordinates
			for (int i = 0; i < N_I; i++)
			{
				inp_I.read(reinterpret_cast<char*>(&xx_I[i]), sizeof(xx_I[i]));
				inp_I.read(reinterpret_cast<char*>(&yy_I[i]), sizeof(yy_I[i]));
				inp_I.read(reinterpret_cast<char*>(&zz_I[i]), sizeof(zz_I[i]));
			}
				
			break;
			
		default:
			std::cerr << "Dimensionality " << dimensionality << " is not supported in read_coordinates!" << std::endl;
			break;
		
	}
	
	
	// close files
	inp_RA.close();	
	inp_I.close();	
}

void NetworkGenerator::read_connections_and_weights(const char* RA_I, const char* I_RA)
{
	std::ifstream inp_RA_I, inp_I_RA;
        
	// open files
	inp_RA_I.open(RA_I, std::ios::binary | std::ios::in);
	inp_I_RA.open(I_RA, std::ios::binary | std::ios::in);

	int N_ra;
	int N_i;

	inp_RA_I.read(reinterpret_cast<char *>(&N_ra), sizeof(N_ra));
	inp_I_RA.read(reinterpret_cast<char *>(&N_i), sizeof(N_i));

	if (N_ra != N_RA)
		std::cerr << "Number of HVC(RA) neurons read from file with fixed synapses: N_ra = " << N_ra << "is different from N_RA = " << N_RA << std::endl;
	
	if (N_i != N_I)
		std::cerr << "Number of HVC(I) neurons read from file with fixed synapses: N_i = " << N_i << "is different from N_I = " << N_I << std::endl;
	
	// resize connections
    weights_RA_I.resize(N_RA);
	syn_ID_RA_I.resize(N_RA);

	weights_I_RA.resize(N_I);
	syn_ID_I_RA.resize(N_I);
	
	// read connections from RA to I
	for (int i = 0; i < N_RA; i++)
	{
		int n_id; // neuronal id
		int size; // number of outgoing connections

		inp_RA_I.read(reinterpret_cast<char *>(&n_id), sizeof(n_id));
		inp_RA_I.read(reinterpret_cast<char *>(&size), sizeof(size)); // write neuron's ID

		syn_ID_RA_I[i].resize(size);
		weights_RA_I[i].resize(size);

		for (int j = 0; j < size; j++)
		{
			inp_RA_I.read(reinterpret_cast<char *>(&syn_ID_RA_I[i][j]), sizeof(syn_ID_RA_I[i][j]));
			inp_RA_I.read(reinterpret_cast<char *>(&weights_RA_I[i][j]), sizeof(weights_RA_I[i][j]));

		}

	}

	// read connections from I to RA
	for (int i = 0; i < N_I; i++)
	{
		int n_id; // neuronal id
		int size; // number of outgoing connections

		inp_I_RA.read(reinterpret_cast<char *>(&n_id), sizeof(n_id));
		inp_I_RA.read(reinterpret_cast<char *>(&size), sizeof(size)); // write neuron's ID
		
		syn_ID_I_RA[i].resize(size);
		weights_I_RA[i].resize(size);
		
		for (int j = 0; j < size; j++)
		{
			inp_I_RA.read(reinterpret_cast<char *>(&syn_ID_I_RA[i][j]), sizeof(syn_ID_I_RA[i][j])); // write targets ID
			inp_I_RA.read(reinterpret_cast<char *>(&weights_I_RA[i][j]), sizeof(weights_I_RA[i][j])); // write targets conductance

		}
	}
	// close files
	inp_I_RA.close();
	inp_RA_I.close();
}

void NetworkGenerator::read_training_neurons(const char* filename)
{
	std::ifstream inp;
  
	inp.open(filename, std::ios::in | std::ios::binary);

	inp.read(reinterpret_cast<char *>(&N_TR), sizeof(N_TR));
	
	//std::cout << "Number of training neurons read from file with training neurons N_TR = " << N_TR << std::endl; 
				 
	training_neurons.resize(N_TR);
	
	for (size_t i = 0; i < N_TR; i++)
		inp.read(reinterpret_cast<char *>(&training_neurons[i]), sizeof(training_neurons[i]));
	
	inp.close();
	
}

void NetworkGenerator::write_invariant_network_to_directory(std::string outputDirectory)
{
	std::string filename_I_coordinates = outputDirectory + "I_xy.bin";
	std::string filenameTraining = outputDirectory + "training_neurons.bin";
	
	this->write_coordinates_I(filename_I_coordinates.c_str());
    //this->write_training_neurons(filenameTraining.c_str());
}

void NetworkGenerator::write_alterable_network_to_directory(std::string extension, std::string outputDirectory)
{
	std::string filename_RA_coordinates = outputDirectory + "RA_xy" + extension + ".bin";
	std::string filename_RA_I = outputDirectory + "RA_I_connections" + extension + ".bin";
	std::string filename_I_RA = outputDirectory + "I_RA_connections" + extension + ".bin";
	std::string filePajekFixed = outputDirectory + "fixed" + extension + ".net";
	
	this->write_coordinates_RA(filename_RA_coordinates.c_str());
    
    this->write_invariable_synapses(filename_RA_I.c_str(), filename_I_RA.c_str());
    this->write_pajek_fixed(filePajekFixed.c_str());
    
}

void NetworkGenerator::write_invariable_synapses(const char* RA_I, const char* I_RA)
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
		size = syn_ID_RA_I[i].size();

		out_RA_I.write(reinterpret_cast<char *>(&i), sizeof(i));
		out_RA_I.write(reinterpret_cast<char *>(&size), sizeof(size)); // write neuron's ID

		for (int j = 0; j < size; j++)
		{
			int k = syn_ID_RA_I[i][j];
			double G = weights_RA_I[i][j];

			out_RA_I.write(reinterpret_cast<char *>(&k), sizeof(k));
			out_RA_I.write(reinterpret_cast<char *>(&G), sizeof(G));

		}

	}

	// write connections from I to RA
	for (int i = 0; i < N_I; i++)
	{
		out_I_RA.write(reinterpret_cast<char *>(&i), sizeof(i)); // write neuron's ID number

		size = syn_ID_I_RA[i].size();
		out_I_RA.write(reinterpret_cast<char *>(&size), sizeof(size)); // write number of targets a neuron has
		for (int j = 0; j < size; j++)
		{
				int k = syn_ID_I_RA[i][j];
				double G = weights_I_RA[i][j];
				out_I_RA.write(reinterpret_cast<char *>(&k), sizeof(k)); // write targets ID
				out_I_RA.write(reinterpret_cast<char *>(&G), sizeof(G)); // write targets conductance

		}
	}
	// close files
	out_I_RA.close();
	out_RA_I.close();
}

void NetworkGenerator::write_pajek_fixed(const char* filename)
{
	std::ofstream out;
	out.open(filename, std::ios::out);

	out << "*Vertices " << N_RA + N_I << "\n";
	
	// check dimensionality of a network
	switch (dimensionality)
	{
		case 1:
			for (int i = 0; i < N_RA; i++)
			{
				// if training neuron, paint with Yellow, otherwise with Green color
				std::vector<int>::iterator iter = std::find(training_neurons.begin(), training_neurons.end(), i);
				
				if (iter != training_neurons.end())
					out << i + 1 << " \"" << i << "\" " << xx_RA[i] << " ic Yellow\n";
				else
					out << i + 1 << " \"" << i << "\" " << xx_RA[i] <<  " ic Green\n";
			}

			for (int i = 0; i < N_I; i++)
			{

				out << i + N_RA + 1 << " \"" << i + N_RA << "\" " << xx_I[i] << " ic Red\n";
			}

			break;
		
		
		case 2:
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
				{

					out << i + N_RA + 1 << " \"" << i + N_RA << "\" " << xx_I[i] << " " << yy_I[i] << " ic Red\n";
				}
			break;
		
		case 3:
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
			{

				out << i + N_RA + 1 << " \"" << i + N_RA << "\" " << xx_I[i] << " " << yy_I[i] << " " << zz_I[i] <<" ic Red\n";
			}
			break;
		
		default:
			std::cerr << "Dimensionality " << dimensionality << " is not supported!" << std::endl;
			break;
	}
	

	out << "*Arcs\n";
	for (int i = 0; i < N_RA; i++)
	{
		for (int j = 0; j < syn_ID_RA_I[i].size(); j++)
			{
				int syn_ID = syn_ID_RA_I[i][j] + N_RA;
				out << i + 1 << " " << syn_ID + 1 << " " << weights_RA_I[i][j] << " c Green\n";
			}
	}
	for (int i = 0; i < N_I; i++)
	{
		for (int j = 0; j < syn_ID_I_RA[i].size(); j++)
		{
			int syn_ID = syn_ID_I_RA[i][j];
			out << i + N_RA + 1 << " " << syn_ID + 1 << " " << weights_I_RA[i][j] << " c Red\n";
		}
	}
}

void NetworkGenerator::write_coordinates_RA(const char* filename)
{
	std::ofstream out;

	// open files
	out.open(filename, std::ios::binary | std::ios::out);

	// write number of neurons
	out.write(reinterpret_cast<char*>(&N_RA), sizeof(N_RA));

	// write coordinates;
	// check dimensionality of a network
	switch (dimensionality)
	{
		case 1:
			for (int i = 0; i < N_RA; i++)
				out.write(reinterpret_cast<char*>(&xx_RA[i]), sizeof(xx_RA[i]));
			break;
		
		
		case 2:
			for (int i = 0; i < N_RA; i++)
			{
				out.write(reinterpret_cast<char*>(&xx_RA[i]), sizeof(xx_RA[i]));
				out.write(reinterpret_cast<char*>(&yy_RA[i]), sizeof(yy_RA[i]));
			}
			break;
		
		case 3:
			for (int i = 0; i < N_RA; i++)
			{
				out.write(reinterpret_cast<char*>(&xx_RA[i]), sizeof(xx_RA[i]));
				out.write(reinterpret_cast<char*>(&yy_RA[i]), sizeof(yy_RA[i]));
				out.write(reinterpret_cast<char*>(&zz_RA[i]), sizeof(zz_RA[i]));
			}
			break;
		
		default:
			std::cerr << "Dimensionality " << dimensionality << " is not supported!" << std::endl;
			break;
	}
	
	// close files
	out.close();
	
}

void NetworkGenerator::write_coordinates_I(const char* filename)
{
	std::ofstream out;

	// open files
	out.open(filename, std::ios::binary | std::ios::out);

	// write number of neurons
	out.write(reinterpret_cast<char*>(&N_I), sizeof(N_I));

	// check dimensionality of a network
	switch (dimensionality)
	{
		case 1:
			for (int i = 0; i < N_I; i++)
				out.write(reinterpret_cast<char*>(&xx_I[i]), sizeof(xx_I[i]));
			break;
		
		
		case 2:
			for (int i = 0; i < N_I; i++)
			{
				out.write(reinterpret_cast<char*>(&xx_I[i]), sizeof(xx_I[i]));
				out.write(reinterpret_cast<char*>(&yy_I[i]), sizeof(yy_I[i]));
			}
			break;
		
		case 3:
			for (int i = 0; i < N_I; i++)
			{
				out.write(reinterpret_cast<char*>(&xx_I[i]), sizeof(xx_I[i]));
				out.write(reinterpret_cast<char*>(&yy_I[i]), sizeof(yy_I[i]));
				out.write(reinterpret_cast<char*>(&zz_I[i]), sizeof(zz_I[i]));
			}
			break;
		
		default:
			std::cerr << "Dimensionality " << dimensionality << " is not supported!" << std::endl;
			break;
	}

	// close files
	out.close();	
}

void NetworkGenerator::write_training_neurons(const char* filename)
{
	std::ofstream out;
  
	out.open(filename, std::ios::out | std::ios::binary );
	
	if (static_cast<int>(training_neurons.size()) != N_TR)
		std::cerr << "Size of vector with training neurons = " << training_neurons.size() << " is different from N_TR = " << N_TR << std::endl;
	else
	{
		out.write(reinterpret_cast<char *>(&N_TR), sizeof(N_TR));

		for (size_t i = 0; i < N_TR; i++)
			out.write(reinterpret_cast<char *>(&training_neurons[i]), sizeof(training_neurons[i]));
	}
		
	out.close();

}

void NetworkGenerator::write_configuration_to_directory(std::string outputDirectory)
{
	std::string filename = outputDirectory + "network_parameters.cfg";
	
	config.write_configuration(filename.c_str());
	
}
