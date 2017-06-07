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
	SIDE = spatial_params.SIDE;
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

void NetworkGenerator::get_network(std::vector<std::vector<int>>* id_RA_I, std::vector<std::vector<double>>* w_RA_I, 
					std::vector<std::vector<int>>* id_I_RA, std::vector<std::vector<double>>* w_I_RA, std::vector<int>* training)
{	
	*id_RA_I = syn_ID_RA_I;
	*w_RA_I = weights_RA_I;
	*id_I_RA = syn_ID_I_RA;
	*w_I_RA = weights_I_RA;
	*training = training_neurons;
}

void NetworkGenerator::generate_default_network(int N_ra, int N_tr, int N_i)
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
	
	training_neurons.resize(N_TR);
	
	this->initialize_coordinates();
	this->initialize_connections();	
}

void NetworkGenerator::generate_network_with_clustered_training(int N_ra, int N_tr, int N_i)
{
	this->generate_default_network(N_ra, N_tr, N_i);
	
	// find neurons that are close to neuron 0 to from a training set cluster
	// select neuron 0 as the first training neuron and find distances to all other neurons in the pool
	double xx = xx_RA[0];
	double yy = yy_RA[0];

	std::vector<double> distances_to_pool_neurons(xx_RA.size()-1); 

	for (size_t i = 1; i < xx_RA.size(); i++)
		distances_to_pool_neurons[i-1] = distance(xx, yy, xx_RA[i], yy_RA[i]);

	// sort distances
	std::vector<size_t> idx(distances_to_pool_neurons.size()); // vector with sorted indices of distances

	std::iota(idx.begin(), idx.end(), 0);

	std::sort(idx.begin(), idx.end(), [&distances_to_pool_neurons](size_t i1, size_t i2)
										{return distances_to_pool_neurons[i1] < distances_to_pool_neurons[i2];});

	// select neuron 0 and N_TR - 1 closest neurons to be training set
	
	training_neurons[0] = 0;
	
	for (int i = 1; i < N_TR; i++) 
		training_neurons[i] = idx[i-1] + 1;
}

void NetworkGenerator::generate_network_with_dispersed_training(int N_ra, int N_tr, int N_i)
{
	this->generate_default_network(N_ra, N_tr, N_i);
	
	// fix training neuron 0 and move other training neurons if necessary
	training_neurons[0] = 0;
	
	int neuron_id_to_check = 1; // id of the next neuron to check
	bool found_next_distant_neuron;

	for (int i = 1; i < N_TR; i++)
	{
		
		// calculate distances from next_neuron_to_check to previously confirmed distant training neurons
		do 
		{
			found_next_distant_neuron = true;
			
			for (int j = 0; j < i; j++)
				if (distance(xx_RA[j], yy_RA[j], xx_RA[neuron_id_to_check], yy_RA[neuron_id_to_check]) < 5 * SIDE / (sqrt(N_I)+1))
				{
					
					std::cout << "Neuron " << neuron_id_to_check << std::endl;
						
					std::cout << "Distance to neuron " << j << " = " << distance(xx_RA[j], yy_RA[j], xx_RA[neuron_id_to_check], yy_RA[neuron_id_to_check]) 
							<< " is smaller than " << 5 * SIDE / (sqrt(N_I) + 1) << std::endl;

					neuron_id_to_check++;
					found_next_distant_neuron = false;
					break;
				}
		} while(!found_next_distant_neuron);
		
		training_neurons[i] = neuron_id_to_check;
		
		neuron_id_to_check++;	
	}
}

void NetworkGenerator::replace_neurons(const std::vector<int>& neurons_to_replace, std::string extension, std::string outputDirectory)
{
	this->initialize_coordinates_for_replaced_neurons(neurons_to_replace);
	this->initialize_connections_for_replaced_neurons(neurons_to_replace);
	
	this->write_alterable_network_to_directory(extension, outputDirectory);
}

void NetworkGenerator::initialize_coordinates()
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
				xx = (double) (i+1) * SIDE / (sqrt(N_I)+1) + noise_generator->random(0.25 * SIDE / (sqrt(N_I) + 1));
				yy = (double) (k+1) * SIDE / (sqrt(N_I)+1) + noise_generator->random(0.25 * SIDE / (sqrt(N_I) + 1));

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
			xx = 0.5*SIDE / (sqrt(N_I)+1) + noise_generator->random(SIDE - SIDE / (sqrt(N_I)+1));
			yy = 0.5*SIDE / (sqrt(N_I)+1) + noise_generator->random(SIDE - SIDE / (sqrt(N_I)+1));

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
}

void NetworkGenerator::initialize_connections()
{
	// connections for HVC(RA) neurons
	for (int i = 0; i < N_RA; i++)
	{
		for (int j = 0; j < N_I; j++)
		{
			 if (noise_generator->random(1) < p_RA2I(i,j))
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
			 if (noise_generator->random(1) < p_I2RA(i,j))
			 {
				 double G = this->sample_Gi2e();

				 weights_I_RA[i].push_back(G);
				 syn_ID_I_RA[i].push_back(j);
			 }
		 }
	 }
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
	for (size_t i = 0; i < neurons_to_replace.size(); i++)
	{
		for (int j = 0; j < N_I; j++)
		{
			 if (noise_generator->random(1) < p_RA2I(neurons_to_replace[i],j))
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
			 if (noise_generator->random(1) < p_I2RA(i,neurons_to_replace[j]))
			 {
				 double G = this->sample_Gi2e();

				 weights_I_RA[i].push_back(G);
				 syn_ID_I_RA[i].push_back(neurons_to_replace[j]);
			 }
		 }
	 }
}


void NetworkGenerator::initialize_coordinates_for_replaced_neurons(const std::vector<int>& neurons_to_replace)
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
			xx = 0.5*SIDE / (sqrt(N_I)+1) + noise_generator->random(SIDE - SIDE / (sqrt(N_I)+1));
			yy = 0.5*SIDE / (sqrt(N_I)+1) + noise_generator->random(SIDE - SIDE / (sqrt(N_I)+1));

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

double NetworkGenerator::p_RA2I(int i_RA, int j_I)
{
	double prob;
    double d;
    d = distance(xx_RA[i_RA], yy_RA[i_RA], xx_I[j_I], yy_I[j_I]);

	if (d < MIN_INTERNEURON_DISTANCE)
		return 0;
	else
		prob = A_RA2I * exp(-d * d / (2 * SIGMA_RA2I * SIGMA_RA2I));

	return prob;
}

double NetworkGenerator::p_I2RA(int i_I, int j_RA)
{
	double prob;
    double d;
    d = distance(xx_RA[j_RA], yy_RA[j_RA], xx_I[i_I], yy_I[i_I]);

	if (d < MIN_INTERNEURON_DISTANCE)
		return 0;
	else
		prob = B_I2RA * exp(-d*d / (2 * SIGMA_I2RA * SIGMA_I2RA));

	return prob;
}


void NetworkGenerator::read_network_from_directory(std::string extension, std::string networkDir)
{
	this->read_network_with_weights_from_directory(extension, networkDir);
	this->sample_weights_based_on_connections();
}

void NetworkGenerator::read_network_with_weights_from_directory(std::string extension, std::string networkDir)
{
	std::string fileRA2I = networkDir + "RA_I_connections" + extension + ".bin";
    std::string fileI2RA = networkDir + "I_RA_connections" + extension + ".bin";

	std::string filename_RA_xy = networkDir + "RA_xy" + extension + ".bin";
	std::string filename_I_xy = networkDir + "I_xy.bin";
   
	std::string filename_training = networkDir + "training_neurons.bin";
   
   
    this->read_coordinates(filename_RA_xy.c_str(), filename_I_xy.c_str()); 
    this->read_connections_and_weights(fileRA2I.c_str(), fileI2RA.c_str());
	this->read_training_neurons(filename_training.c_str()); 
    
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
	
	std::cout << "Number of neurons in file with HVC(RA) coordinates N_RA = " << N_RA << std::endl;
	std::cout << "Number of neurons in file with HVC(I) coordinates N_I = " << N_I << std::endl;

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
	
	std::cout << "Number of training neurons read from file with training neurons N_TR = " << N_TR << std::endl; 
				 
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
    this->write_training_neurons(filenameTraining.c_str());
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

	for (int i = 0; i < N_RA; i++)
	{
		// if training neuron, paint with Yellow, otherwise with Green color
		std::vector<int>::iterator iter = std::find(training_neurons.begin(), training_neurons.end(), i);
		
		if (iter != training_neurons.end())
			out << i + 1 << " \"" << i << "\" " << xx_RA[i]/SIDE << " " << yy_RA[i]/SIDE << " ic Yellow\n";
		else
			out << i + 1 << " \"" << i << "\" " << xx_RA[i]/SIDE << " " << yy_RA[i]/SIDE << " ic Green\n";
	}

	for (int i = 0; i < N_I; i++)
	{

		out << i + N_RA + 1 << " \"" << i + N_RA << "\" " << xx_I[i]/SIDE << " " << yy_I[i]/SIDE << " ic Red\n";
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

	// write coordinates
	for (int i = 0; i < N_RA; i++)
	{
		out.write(reinterpret_cast<char*>(&xx_RA[i]), sizeof(xx_RA[i]));
		out.write(reinterpret_cast<char*>(&yy_RA[i]), sizeof(yy_RA[i]));
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

	// write coordinates
	for (int i = 0; i < N_I; i++)
	{
		out.write(reinterpret_cast<char*>(&xx_I[i]), sizeof(xx_I[i]));
		out.write(reinterpret_cast<char*>(&yy_I[i]), sizeof(yy_I[i]));
	}

	// close files
	out.close();
	
}


void NetworkGenerator::write_training_neurons(const char* filename)
{
	std::ofstream out;
  
	out.open(filename, std::ios::out | std::ios::binary );

	out.write(reinterpret_cast<char *>(&N_TR), sizeof(N_TR));

	for (size_t i = 0; i < N_TR; i++)
		out.write(reinterpret_cast<char *>(&training_neurons[i]), sizeof(training_neurons[i]));
	
	out.close();

}

void NetworkGenerator::write_configuration_to_directory(std::string outputDirectory)
{
	std::string filename = outputDirectory + "network_parameters.cfg";
	
	config.write_configuration(filename.c_str());
	
}
