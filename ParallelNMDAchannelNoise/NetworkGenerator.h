#ifndef NETWORK_GENERATOR_H
#define NETWORK_GENERATOR_H

#include "ConfigurationNetworkGenerator.h"
#include <vector>
#include "poisson_noise.h"
#include <cmath>
#include <string>

class NetworkGenerator
{
public:
    NetworkGenerator(ConfigurationNetworkGenerator cfg, Poisson_noise* generator); // initialize with parameters from configuration file
    
    void generate_default_network(int N_ra, int N_tr, int N_i); // generate network by initializing coordinates of HVC(RA) and HVC(I) neurons. 
									 // training neurons are chosen randomly from pool neurons
									 // connections between HVC(RA) and HVC(I) neurons are also initialized
									 
	void generate_network_with_clustered_training(int N_ra, int N_tr, int N_i); // generate network by initializing coordinates of HVC(RA) and HVC(I) neurons. 
									 // training neurons are chosen so that they form spatial cluster
									 // connections between HVC(RA) and HVC(I) neurons are also initialized
									 
									 
	void generate_network_with_dispersed_training(int N_ra, int N_tr, int N_i); // generate network by initializing coordinates of HVC(RA) and HVC(I) neurons. 
									 // training neurons are chosen so that they are far apart in space
									 // connections between HVC(RA) and HVC(I) neurons are also initialized
									 
	void replace_neurons(const std::vector<int>& neurons_to_replace, 
				std::string extension, std::string outputDirectory); // replace neurons with ids in neurons_to_replace
																					// replaced neurons get new coordinates and new
																					// input and output connections. Extension is added to the end of all filenames
																					// Files are written to directory outputDirectory
	
	void get_network(std::vector<std::vector<int>>* id_RA_I, std::vector<std::vector<double>>* w_RA_I, 
					std::vector<std::vector<int>>* id_I_RA, std::vector<std::vector<double>>* w_I_RA, std::vector<int>* training); 
																										 // get connections and weights between HVC(RA) 
																									     // and HVC(I) neurons and ids of training neurons
	
	void get_neuron_numbers(int* N_ra, int* N_tr, int* N_i); // get number of HVC(RA), HVC(I) and training neurons in the network
					 
	void read_network_with_weights_from_directory(std::string extension, std::string networkDir); // read network configuration with all weights as in files
																		   // stored in directory networkDir. Coordinates of HVC(RA) and HVC(I) neurons
																		   // should be in files RA_xy"extension".bin and I_xy.bin. Connections between HVC(RA) and
																		   // HVC(I) neurons should be in files RA_I_connections"extension".bin and I_RA_connections"extension".bin
																		   // ids of training neurons should be in training_neurons.bin file
	
	void read_network_from_directory(std::string extension, std::string networkDir); // read network configuration as in files
															  // stored in directory networkDir. Coordinates of HVC(RA) and HVC(I) neurons
															  // should be in files RA_xy"extension".bin and I_xy.bin. Connections between HVC(RA) and
															  // HVC(I) neurons should be in files RA_I_connections"extension".bin and I_RA_connections"extension".bin
															  // ids of training neurons should be in training_neurons.bin file
															  // After connections are read, weights are sampled according to distributions
															  // read from configuration file
	
	void set_inhibitory_strength(double G_mean, double G_std); // set inhibitory strength parameters. G_mean - average inhibitory strength;
																// G_std - standard deviation of inhibitory strength
	
	void write_invariant_network_to_directory(std::string outputDirectory); // write invariant part of generated network to output directory; 
																			// it includes coordinates of interneurons and set of training neurons
	
	void write_alterable_network_to_directory(std::string extension, std::string outputDirectory); // write alterable part of generated network to output directory; 
																								   // it includes coordinates of HVC(RA) neurons; connections between
																								   // HVC(RA) and interneurons; and pajek file with HVC(RA) <-> HVC(I) network.
																								   // extension is added to filenames
	
	void write_configuration_to_directory(std::string outputDirectory); // write configuration file to directory
	
private:
	ConfigurationNetworkGenerator config; // configuration containing network parameters

	std::vector <double> xx_RA; // x-coordinates of RA neurons
	std::vector <double> yy_RA; // y-coordinates of RA neurons
	std::vector <double> xx_I; // x-coordinates of I neurons
	std::vector <double> yy_I; // y-coordinates of I neurons

	Poisson_noise* noise_generator; // noise generator
	
	std::vector<std::vector<int>> syn_ID_RA_I; // array with synaptic ID numbers from RA to I neurons
	std::vector<std::vector<int>> syn_ID_I_RA; // array with synaptic ID numbers from I to RA neurons
	
	std::vector<std::vector<double>> weights_RA_I; // array with synaptic weights from RA to I neurons
	std::vector<std::vector<double>> weights_I_RA; // array with synaptic weights from I to RA neurons
		
	// network params
	int N_TR; // number of training HVC(RA) neurons
	int N_RA; // number of HVC(RA) neurons
	int N_I; // number of HVC(I) neurons
	
	std::vector<int> training_neurons; // array with indices of HVC(RA) training neurons
	
	// connectivity params
	
	double Gei_mean; // mean synaptic weight from HVC(RA) to HVC(I) neuron
    double Gei_std; // standard deviation of synaptic weight from HVC(RA) to HVC(I) neuron
    double Gie_mean; // mean synaptic weight from HVC(I) to HVC(RA) neuron
    double Gie_std; // standard deviation of synaptic weight from HVC(I) to HVC(RA) neuron
    
	double SIDE; // length of HVC side
	double MIN_INTERNEURON_DISTANCE; // minimum distance between neurons
	
	double A_RA2I; // constant for nearby HVC(RA) -> HVC(I) connections 
	double SIGMA_RA2I; // spatial scale of probability of connections decay

	double B_I2RA; // constant for nearby HVC(RA) -> HVC(I) connections
	double SIGMA_I2RA; // spatial scale of probability of connections decay
	
	void initialize_coordinates(); // initialize coordinates of neurons
    void initialize_connections(); // initialize connections for neurons based on Gaussian distributions
        
    void initialize_coordinates_for_replaced_neurons(const std::vector<int>& neurons_to_replace); // change coordinates of neurons with ids in
																								  // neurons_to_replace to new coordinates not within
																								  // MIN_INTERNEURON_DISTANCE from all previous HVC(RA)
																								  // locations and all HVC(I) locations
    
    void initialize_connections_for_replaced_neurons(const std::vector<int>& neurons_to_replace); // change connections of neurons with ids in
																								  // neurons_to_replace to new connections based on their new
																								  // coordinates. Previous input and output connections are erased
																								  
    
    
	double p_RA2I(int i, int j); // probability of connection from RA to I neuron
	double p_I2RA(int i, int j); // probability of connection from I to RA neuron

	void sample_weights_based_on_connections(); // populates weights_RA_I and weights_I_RA arrays with
												// weights sampled from distributions read from configuration file.
												// Connections between HVC(RA) and HVC(I) neurons should be read
												// beforehand with function read_network_from_directory

	double sample_Ge2i(); // sample synaptic weight from HVC(RA) to HVC(I) neuron
    double sample_Gi2e(); // sample synaptic weight from HVC(I) to HVC(RA) neuron
    
   
    void set_spatial_parameters(const struct SpatialParameters& spatial_params); // set spatial parameters
    
    // read from files
    void read_coordinates(const char* RA_xy, const char* I_xy); // read coordinates of HVC(RA) and HVC(I) neurons from files
    void read_connections_and_weights(const char* RA_I, const char* I_RA); // read connections and their weights between HVC(RA) 
																		   // and HVC(I) neurons from files
    
    void read_training_neurons(const char* filename); // read training neurons from file
    
    // write to files
    void write_all_coordinates(const char* RA_xy, const char* I_xy); // write coordinates of HVC(RA) and HVC(I) neurons
    void write_coordinates_RA(const char* filename); // write coordinates of HVC(RA) neurons
    void write_coordinates_I(const char* filename); // write coordinates of HVC(I) neurons
    
	void write_invariable_synapses(const char* RA_I, const char* I_RA); // write RA to I and I to RA connections      
	void write_pajek_fixed(const char* filename); // write fixed synapses to a file for pajek
    
    void write_training_neurons(const char* filename); // write ids of training HVC(RA) neuronal
    
    
};

static double distance(double x1, double y1, double x2, double y2)
{
	double d;
	d = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	return d;
}
	
#endif
