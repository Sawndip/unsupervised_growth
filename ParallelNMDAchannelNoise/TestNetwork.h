#ifndef TEST_NETWORK_H
#define TEST_NETWORK_H

#include "poisson_noise.h"
#include "ConfigurationTestNetwork.h"
#include <boost/circular_buffer.hpp>
#include "HHI_buffer.h"
#include "HH2_buffer.h"
#include <vector>

typedef boost::circular_buffer<int> intBuffer;

#define MAX_NUM_THREADS 24

class TestNetwork
{
public:
	TestNetwork();
	
	void create_second_layer_growth_sim(int N_tr, int num_groups, int num_neurons_in_groups, 
					double Gei, double dGie); // prepare network for simulating 2nd layer growth:
											  // training neurons form 1st layer; remaining pool neurons - 2nd.
											  // There is also 1 interneuron in the network
											  // first training neuron connects to interneuron, which in turn 
											  // connects to all HVC-RA neurons in the 2nd layer;
											  // HVC-RA neurons in 2nd layer are split into different groups of
											  // num_neurons_in_groups neurons in each group. Inhibitory strength of  
											  // connection is dGie * group_id

	void simulate(double trial_duration, double timestep, int save_freq, std::string outputDir); // simulate network growth;
													// trials have duration trial_duration; dynamics timestep is timestep;
													// each save_freq trials data is saved to files in directory
													// outputDir
													
													

	// set parameters
	void set_synaptic_parameters(const struct SynapticParameters &syn_params);
	void set_gaba_parameters(const struct GabaParameters &g_params);
	void set_noise_parameters(double mean_s, double sigma_s, double mean_d, double sigma_d);
	
	
private:
	// parameters
	const static double NETWORK_UPDATE_FREQUENCY; // synchronization frequency of neural dynamics in network
	const static double STDP_WINDOW; // burst events outside of this window do not lead to STDP
	const static double WAIT_TIME; // time to wait before current can be injected to training neurons
	const static double TRAINING_KICK; // strength of excitatory conductance pulse delivered to training neurons 

	
	struct SynapticParameters synaptic_params; // synaptic parameters
	struct GabaParameters gaba_params; // maturation parameters

	double white_noise_mean_soma; // mean value of white noise current injected to somatic   compartment
	double white_noise_std_soma; //  std  value of white noise current injected to somatic   compartment
	double white_noise_mean_dend; // mean value of white noise current injected to dendritic compartment
	double white_noise_std_dend; //  std  value of white noise current injected to dendritic compartment

	// data
	int N_RA; // number of HVC-RA neurons
	int N_TR; // number of training HVC-RA neurons
	int N_I; // number of HVC-I neurons
	
	std::vector<HH2_buffer> HVCRA; // array with HVC-RA neurons
	std::vector<HHI_buffer> HVCI;  // array with HVC-I  neurons

	Poisson_noise noise_generators[MAX_NUM_THREADS]; // noise generators
	
	// HVC-RA <-> HVC-I connections
	std::vector<std::vector<int>> i_syn_ID_RA2I; // synaptic ids of HVC-RA -> HVC-I  connections
	std::vector<std::vector<int>> i_syn_ID_I2RA; // synaptic ids of HVC-I  -> HVC-RA connections
	
	std::vector<std::vector<double>> d_weights_RA2I; // weights for HVC-RA -> HVC-I  connections
	std::vector<std::vector<double>> d_weights_I2RA; // weights for HVC-I  -> HVC-RA connections
	
	// HVC-RA <-> HVC-RA connections
	std::vector<std::vector<int>> i_activeSynapses_id; // target ids of active synapses
	std::vector<std::vector<int>> i_superSynapses_id; // target ids of super synapses
	
	std::vector<std::vector<double>> d_weights; // weights between HVC-RA neurons
	
	std::vector<std::vector<bool>> b_active; // indicator that synapse is active
	std::vector<std::vector<bool>> b_super; // indicator that synapse is super
	
	// arrays with spike times in trial
    std::vector<std::vector<double>> d_spikes_in_trial_soma; // somatic spike times of HVC-RA neurons in recent trial
    std::vector<std::vector<double>> d_spikes_in_trial_dend; // dendritic spike times of HVC-RA neurons in recent trial
    std::vector<std::vector<double>> d_spikes_in_trial_interneurons; // spike times of HVC-I neurons in recent trial
    
    // arrays for STDP and maturation processes
    std::vector<double> d_last_dend_spike; // most recent dendritic spike times of HVC-RA neurons
    std::vector<double> d_E_gaba; // current GABA reversal potential of HVC-RA neurons
    std::vector<bool> b_remodeled; // indicator that HVC-RA neuron is in axon-remodeled state
    std::vector<intBuffer> i_num_bursts_in_recent_trials; // vector contains number of dendritic bursts in recent trials within RATE_WINDOW_LONG

    // misc functions
    void resize_arrays(); // resize arrays for HVC-RA and HVC-I neurons
    void reset_after_trial(double trial_duration); // reset dynamics after trial: all stored spikes are deleted,
												   // trial_duration is subtracted from all last dend spikes of HVC-RA neurons
    
    // simulation
    void trial(double trial_duration, double timestep); // run growth trial
    
    // STDP support
    void LTP_burst(double &w, double dt); // LTP rule for bursts
    void LTD_burst(double &w, double dt); // LTD rule for bursts
    void potentiation_decay(); // apply potentiation decay to all synapses of HVC-RA -> HVC-RA connections
    void axon_remodeling(int i); // apply axon remodeling to HVC-RA neuron i
    void update_synapse(int i, int j); // update synapse state of i -> j connection between HVC-RA neurons
    void update_all_synapses(); // update synapse states for all HVC-RA -> HVC-RA connections
    
    // maturation support
    void update_Ei(); // update GABA reversal potential of HVC-RA neurons based on recent firing history
    
    // write data
    void write_time_weights(int trial_number, const char* filename); // write dynamics of synaptic weights 
																	 // for HVC-RA -> HVC-RA connections
																	 
	void write_neuron_states(int trial_number, const char* filename); // write firing rate, GABA potential and remodeled
																	  // state of HVC-RA neurons to a file
    
	
};

#endif
