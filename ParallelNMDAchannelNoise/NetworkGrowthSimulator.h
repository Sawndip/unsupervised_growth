#ifndef NETWORK_GROWTH_SIMULATOR
#define NETWORK_GROWTH_SIMULATOR

#include <vector>
#include <cmath>
#include <mpi.h>
#include <string>

#include "poisson_noise.h"
#include <boost/circular_buffer.hpp>
#include <boost/dynamic_bitset.hpp>
#include "ConfigurationNetworkGrowth.h"
#include "ConfigurationNetworkTopology.h"

#include "HHI_buffer.h"
#include "HH2_test.h"

using std::vector;
using std::pair;

typedef boost::circular_buffer<int> intBuffer;
typedef boost::dynamic_bitset<> bitArray;

class NetworkGrowthSimulator
{
public:
	NetworkGrowthSimulator();
	
	void generate_network_topology(int N_ra, int N_i, int N_tr, 
								const struct TopologyParameters &con_params, std::string networkDir); // generate network topology and write the results to directory networkDir
	
    void read_network_state(std::string dirname, int starting_trial = -1); // read network state from files in the directory dirname

    void initialize_test_connections(int num_RA_targets, int num_RA_target_groups); // initialize test connections: first training neuron is connected to the first
                                                                                    // interneuron; First interneuron in turn is connected to num_RA_targets RA
                                                                                    // neurons. Strength if inhibitory connection increases by Gie_mean in each
                                                                                    // next group of RA neurons
	
	void disable_RA2I_connections(); // disable all connections from HVC(RA) to HVC(I) neurons
	    
	
    void chain_growth(bool training, int save_freq_short, int save_freq_long,
						std::string outputDirectory); // run chain growth  
                                                                       // NOTE: coordinates and connections MUST be initialized before using chain_growth

	
	void chain_growth_with_inhibition_tracking(bool training, int save_freq_long, int num_trials, 
											double time_resolution_conductance, std::string outputDirectory);
											// perform chain growth trials with tracking inhibitory conductance of HVC-RA neurons
	
	void new_chain_growth(const ConfigurationNetworkGrowth &cfg, std::string networkDirectory, std::string fileTraining,  std::string outputDirectory, 
															bool training, int save_freq_short, int save_freq_long);
	
	void new_chain_growth_with_inhibition_tracking(const ConfigurationNetworkGrowth &cfg, std::string networkDirectory, std::string fileTraining,  std::string outputDirectory, 
															bool training, int save_freq_long, int num_trials, double time_resolution_conductance);
															// start new chain growth simulation with inhibitory conductance tracking
															//  at each simulation trial
	
	void continue_chain_growth(std::string dataDir, std::string outDir, int starting_trial, bool training, int save_freq_short, int save_freq_long); // continue chain growth using network state defined by files in 
																						// directory dataDir from trial starting_trial

	void test_inhAndExc_response(const struct ConnectionParameters &con_par,
													const struct NoiseParameters &n_par,
													std::string networkDirectory, std::string fileTraining,
													std::string outputDirectory); // test neuron responses to inhibitory and excitatory inputs

	void test_synthetic_chain(const ConfigurationNetworkGrowth &cfg, int num_layers, int num_group, double probability, double Gee,
										std::string networkDirectory, int num_trials, std::string outputDirectory); // test pruned chain topology with 
																									// "num_layers" layers and "num_group" neurons in each layer
																									// neurons connect to neurons in the next layer with
																									// probability "probability" and stregnth "Gee"

	void test_chain_recovery(std::string dataDir, int starting_trial,
				double fraction, bool training, int save_freq_short, int save_freq_long); // make lesion by killing fraction of chain neurons. Then continue growth to see how network recovers
				
	void test_chain(std::string networkDirectory,  int starting_trial, int num_trials, std::string outputDirectory); 
															// test grown chain
															// data should reside in directory networkDirectory; network state corresponds to trial starting_trial
															// num_trials simulations are made; 
															// somatic, dendritic, interneuron spikes, and jitter are written to directory outputDirectory

	void initialize_network();

    void print_invariable_connections(); // get number of connections for invariable synapses
    void print_received_invariable_connections(); // show received invariable connections

	void print_simulation_parameters(); // print simulation parameters
protected:
		struct ConnectionParameters connection_params; // struct with network parameters: inhibition strength and axonal delays
		struct TopologyParameters topology_params; // structure with network topology parameters;
		struct SynapticParameters synaptic_params; // structure with synaptic STDP parameters;
		struct NoiseParameters noise_params; // structure with noise paramters for HVC-RA neurons
		struct MaturationParameters maturation_params; // structure with maturation parameters for HVC-RA neurons
		
        // number of neurons
		int N_TR; // number of training HVC(RA) neurons
		int N_RA; // number of HVC(RA) neurons
		int N_I; // number of HVC(I) neurons
		int N_RA_local; // number of RA neurons in each process
		int N_I_local; // number of I neurons in each process
		
		vector<HH2_test> HVCRA_local; // array of HVC(RA) neurons
		vector<HHI_buffer> HVCI_local; // array of HVC(I) neurons

		// coordinates info
		vector <double> xx_RA; // x-coordinates of HVC-RA neurons
		vector <double> yy_RA; // y-coordinates of HVC-RA neurons
		vector <double> zz_RA; // z-coordinates of HVC-RA neurons
		
		vector <double> xx_I; // x-coordinates of HVC-I neurons
		vector <double> yy_I; // y-coordinates of HVC-I neurons
		vector <double> zz_I; // z-coordinates of HVC-I neurons
		
		
		double MODEL_INTERNEURON_DISTANCE; // average distance between neighbouring interneurons in the model 
		
		const static double MATURE_SYNAPSE_SCALING; // constant to scale excitatory synapses to mature neurons
		const static double WAITING_TIME; // time in the beginning and end of the trial when no current injection happens
		const static double STDP_WINDOW; // time window in which spikes contribute to STDP change in weight
		const static double TIMESTEP; // time step for dynamics
		const static double TRIAL_DURATION; // duration of trial in milliseconds
		const static double NETWORK_UPDATE_FREQUENCY; // how often network state should be updated in ms
		const static double EXPERIMENTAL_INTERNEURON_DISTANCE; // experimental value for distance between two interneurons in HVC in microns
		const static double EXPERIMENTAL_AXONAL_DELAY_PER_MICRON; // experimental value for axonal delay in ms per micron distance
		
		const static double MIN_INTERSOMATIC_DISTANCE; // minimum allowed distance between somatas of two neurons in the model
		const static double G_TRAINING_KICK; // strength of excitatory conductance kick delivered to training HVC-RA neurons
		const static double MATURATION_SCALE_SPONTANEOUS; // timescale of neuron properties of spontaneous maturation
		const static double MATURATION_SCALE_DRIVEN; // timescale of neuron properties of driven maturation

		double internal_time; // network
        double network_time; // network time
        
		double current_injection_time; // injection time of training current
		int trial_number; // number of simulated trials

		Poisson_noise noise_generator; // poisson noise generator

		vector<vector<double>> Ginh_global; // array of inhibitory conductances of all HVC-RA neurons on master process
		vector<vector<double>> Ginh_local; // array of inhibitory conductances of local HVC-RA neurons

		vector<vector<double>> weights_RA_RA_global; // array of synaptic strength of all connections between RA neurons
        vector<vector<double>> weights_RA_RA_local; // array of synaptic strength of all connections between RA neurons
		
        vector<vector<double>> weights_RA_I_global; // array with synapses from RA to I neurons
		vector<vector<double>> weights_I_RA_global; // array with senapses from I to RA neurons
		vector<vector<double>> weights_RA_I_local; // array with synapses from RA to I neurons
		vector<vector<double>> weights_I_RA_local; // array with senapses from I to RA neurons
        
        vector<vector<int>> syn_ID_RA_I_global; // array with synaptic ID numbers from RA to I neurons
		vector<vector<int>> syn_ID_I_RA_global; // array with synaptic ID numbers from I to RA neurons
		vector<vector<int>> syn_ID_RA_I_local; // array with synaptic ID numbers from RA to I neurons
		vector<vector<int>> syn_ID_I_RA_local; // array with synaptic ID numbers from I to RA neurons

		vector<vector<double>> syn_lengths_RA_I_global; // array with connection lengths for synapses from HVC-RA to HVC-I neurons
		vector<vector<double>> syn_lengths_I_RA_global; // array with connection lengths for synapses from HVC-I to HVC-RA neurons
		vector<vector<double>> syn_lengths_RA_RA_global; // array with connection lengths for synapses from HVC-RA to HVC-RA neurons

		
		vector<vector<double>> axonal_delays_RA_I_global; // array with axonal delays for synapses from HVC-RA to HVC-I neurons
		vector<vector<double>> axonal_delays_RA_I_local; // array with axonal delays for synapses from HVC-RA to HVC-I neurons
		
		vector<vector<double>> axonal_delays_I_RA_global; // array with axonal delays for synapses from HVC-I to HVC-RA neurons
		vector<vector<double>> axonal_delays_I_RA_local; // array with axonal delays for synapses from HVC-I to HVC-RA neurons
		
		vector<vector<double>> axonal_delays_RA_RA_global; // array with axonal delays for synapses from HVC-RA to HVC-RA neurons
		vector<vector<double>> axonal_delays_RA_RA_local; // array with axonal delays for synapses from HVC-RA to HVC-RA neurons


		vector<bitArray> active_indicators_local; // array of HVC(RA) neurons with active synapses
		vector<bitArray> supersynapses_indicators_local; // indicator array for active supersynapses;
		vector<bitArray> rescaled_indicators_local; // indicator array that shows if synapse was rescaled;

		vector<int> training_neurons; // vector with ids of training neurons
		vector<double> training_spread_times; // vector with time spread of training neurons
		
		vector<vector<int>> active_synapses_local; // array of vectors with IDs of active synapses
		vector<vector<int>> supersynapses_local; // array of vectors with IDs of supersynapses
		
        vector<vector<int>> active_synapses_global; // array of vectors with IDs of active synapses
		vector<vector<int>> supersynapses_global; // array of vectors with IDs of supersynapses
		
		// delivery queues
		vector<vector<pair<double,int>>> delivery_queue_RA_RA_soma; // local queue with spike delivery times for HVC-RA -> HVC-RA interactions
		vector<vector<pair<double,int>>> delivery_queue_RA_I; // local queue with spike delivery times for HVC-RA -> HVC-I interactions
		
		vector<vector<pair<double,int>>> delivery_queue_I_RA; // local queue with spike delivery times for HVC-I -> HVC-RA interactions
		
		vector<vector<vector<double>>> delivered_spike_times; // local array with the times when somatic spikes of 
															 // source HVC-RA neuron were delivered to other HVC-RA neurons
			
		//vector<vector<double>> previous_somatic_spike_times_global; // global array with the previous somatic spike times of neurons that are not too old
		vector<vector<double>> previous_dendritic_spike_times_global; // global array with the previous dendritic spike times of neurons that are not too old
		
		vector<vector<double>> previous_somatic_spike_times_local; // local array with the previous somatic spike times of neurons that are not too old
		
		
        vector<vector<double>> spikes_in_trial_soma_global; // array with somatic spike times in the last trial
		vector<vector<double>> spikes_in_trial_dend_global; // array with dendritic spike times in the last trial
		vector<vector<double>> spikes_in_trial_interneuron_global; // array with interneuron spike times in the last trial

        vector<vector<double>> spikes_in_trial_soma_local; // array with somatic spike times in the last trial
		vector<vector<double>> spikes_in_trial_dend_local; // array with dendritic spike times in the last trial
		vector<vector<double>> spikes_in_trial_interneuron_local; // array with interneuron spike times in the last trial
		
		vector<int> remodeled_local; // local indicators if neuron underwent axon remodelling
		vector<int> remodeled_global; // global indicators if neuron underwent axon remodelling
		
		vector<int> mature_local; // local indicators if neuron became mature
		vector<int> mature_global; // global indicators if neuron became mature
		
		vector<int> maturation_scale_local; // local rate of maturation
		vector<int> maturation_scale_global; // global rate of maturation
		
		///////////////////////////////////
		//// Maturation arrays ////////////
		///////////////////////////////////
		
		vector<double> gaba_potential_local; // array with local values of GABA reverse potential
		vector<double> gaba_potential_global; // array with global values of GABA reverse potential
		
		vector<double> rest_potential_local; // array with local values of rest potential
		vector<double> rest_potential_global; // array with global values of rest potential
		
		vector<double> Gk_local; // array with local values of potassium somatic conductance
		vector<double> Gk_global; // array with global values of potassium somatic conductance
		
		vector<double> GNa_local; // array with local values of sodium somatic conductance
		vector<double> GNa_global; // array with global values of sodium somatic conductance
		
		vector<double> GCa_local; // array with local values of calcium dendritic conductance
		vector<double> GCa_global; // global array with calcium dendritic conductances
		
		
		vector<double> Ad_local; // array with local values of dendritic area
		vector<double> Ad_global; // array with global values of dendritic area
		
		vector<double> Rc_local; // array with local values of coupling resistance between soma and dendrite
		vector<double> Rc_global; // array with global values of coupling resistance between soma and dendrite
		
		//////////////////////////////////////////
		//// Keep track of maturation ////////////
		//////////////////////////////////////////
		
		
		vector<double> firing_rate_short_local; // array with local firing rates of HVC(RA) neurons in small window
		vector<double> firing_rate_short_global; // array with global firing rates of HVC(RA) neurons in small window
		vector<double> firing_rate_long_local; // array with local firing rates of HVC(RA) neurons in large window
		vector<double> firing_rate_long_global; // array with global firing rates of HVC(RA) neurons in large window
    
        vector<int> gaba_reached_mature_local; // local indicators for neurons that GABA potential reached mature value

        vector<intBuffer> num_spikes_in_recent_trials_local; // array of circular buffer containing recent neuron rates
		vector<vector<int>> num_spikes_in_recent_trials_global; // global array with the number of spikes occured in previoud RATE_WINDOW_LONG trials

		vector<int> Id_RA_local; // local array with id of RA neurons in each process
		vector<int> Id_RA_global; // global array with ids of RA neurons in all processes
		vector<int> Id_I_local; // Id of I neurons in each process

        // neuron's internal time after previous replacement
        vector<int> num_trials_after_replacement_local;
        vector<int> num_trials_after_replacement_global;


		////////////////
		//// Network testing
		////////////////////
		void generate_synthetic_chain(int num_layers, int num_group, double probability, double Gee); // generate pruned synfire chain topology
		
		void calculate_and_write_inhAndExc(const std::vector<std::vector<double>> &spike_times,
											  const std::vector<std::vector<double>> &relevant_spike_times,
											  const std::vector<std::vector<int>> &num_somatic_spikes_in_trials,
											  const std::vector<std::vector<int>> &num_relevant_spikes_in_trials,
											  const char *filename); // calculate and write results of inhibitory and excitatory inputs test
											  
		void calculate_and_write_jitter(int num_trials, const std::vector<std::vector<double>> &first_soma_spike_times,
											  const std::vector<std::vector<double>> &first_dend_spike_times,
											  const std::vector<std::vector<int>> &num_somatic_spikes_in_trials,
											  const std::vector<std::vector<int>> &num_dendritic_spikes_in_trials,
											  std::string outputDirectory); // calculate jitter and write it to a file in directory
		///////////////////////////////////
		//// Maturation
		///////////////////////////////////
		void check_neuron_activity(std::vector<int>& neurons_to_replace); // check if neurons became driven or silent 
        void update_neuron_properties(); // update all neuron properties based on their age
        void update_neuron_properties_sameDendrite(); // update all neuron properties based on their age for same dendrite size
        void update_neuron_properties_sameDendrite_diffMaturationRate();// update all neuron properties based on their age for same dendrite size;
																		// maturation rate can be different for different neurons
																		
        void set_neuron_properties_sudden_maturation(); // set neuron properties according to mature indicators
        void set_neuron_properties(); // set neuron properties according to values in local arrays
        void set_training_neurons_mature(); // set mature parameters for training  HVC-RA neurons
		void set_all_neurons_immature(); // set immature parameters for all HVC-RA neurons
		void set_neuron_mature(int local_id); // set neuron with local_id mature
		void set_neuron_immature(int local_id); // set neuron witl local_id immature
		
        // current
		void set_training_current(double t); // set current to training neurons. t - current injection time.
        void set_training_current(); // set current to training neurons
		
	    // noise generator	
	    void initialize_generator(); // initialize generator for processes
        
        // trials
        void trial_no_stdp_fixedSpread(std::vector<double>& spread_times); // trial with no stdp and fixed spread of training neurons
        void trial_no_stdp(double training_kick_time); // trial with no stdp, used for testing grown chains 
	    void trial_somatic_stdp_no_delays(bool training); // single trial with STDP rule based on somatic spikes. No axonal delays
	    
	    void trial_event_pre_dend_post_delays_sudden_maturation(bool training); // singe trial with STDP rule; somatic spikes are presynaptic
																				// events; dendritic spikes - postsynaptic. All somatic spikes
																				// occured within certain window are considered as a single event
		
		void trial_burst_pre_dend_post_delays_sudden_maturation_noImmatureOut(bool training); // singe trial with STDP rule; bursts are presynaptic
																				// events; dendritic spikes - postsynaptic. Immature neurons have no output connections
		
		void trial_burst_pre_dend_event_post_delays_sudden_maturation_noImmatureOut(bool training); // singe trial with STDP rule; bursts are presynaptic
																				// events; dendritic events (dendritic spikes within event_window) - postsynaptic. 
																				// Immature neurons have no output connections
		
		
		void trial_burst_pre_dend_event_post_delays_sudden_maturation_noImmatureOut_fixedSpread(bool training, std::vector<double>& spread_times);
														// // singe trial with STDP rule; bursts are presynaptic
																				// events; dendritic spikes - postsynaptic. Immature neurons have no output connections
																				// training neurons have fixed spread
		
		void trial_1stSoma_pre_1stSoma_post_delays_fixedSpread_with_inhibition_tracking(bool training, std::vector<double>& spread_times, 
																					double time_resolution_conductance);
															// single trial with STDP rule; presynaptic event - delivered first somatic spike
															// postsynaptic event - fired first somatic spike
															// each trial inhibitory conductance at time resolution time_resolution_conductance
															// of all HVC-RA neurons is saved to local arrays
															
		
		void trial_1stSoma_pre_1stSoma_post_delays_fixedSpread(bool training, std::vector<double>& spread_times);
													// single trial with STDP rule; presynaptic event - delivered first somatic spike
													// postsynaptic event - fired first somatic spike
		
		void trial_noImmatureOut_fixedSpread(bool training, std::vector<double>& spread_times); // singe trial with STDP rule; bursts are presynaptic
																				// events; dendritic events (dendritic spikes within event_window) - postsynaptic. Immature neurons have no output connections
																				// training neurons have fixed spread
																		
	    void trial_soma_pre_dend_post_stdp_no_delays(bool training); // single trial with STDP rule with somatic spikes as presynaptic
																	 // events and dendritic spikes as postsynaptic. No axonal delays
		
		void trial_soma_pre_dend_post_stdp_delays(bool training); // single trial with STDP rule with somatic spikes as presynaptic
																	 // events and dendritic spikes as postsynaptic. With axonal delays
		
		void trial_soma_pre_dend_post_stdp_delays_sudden_maturation(bool training); // single trial with STDP rule with somatic spikes as presynaptic
																	 // events and dendritic spikes as postsynaptic. With axonal delays
																	 // Neurons instantaneously change properties to mature once they fire robustly
												 
	    void trial_burst_stdp(int training); // make one trial with STDP rules applied to dendritic bursts
		void mature_trial(); // simulation trial without STDP rules
        void test_mature_chain(int num_trials); // test of mature network
	    
        void randomize_after_trial(); // set all neurons to the resting state
        void reset_after_chain_test(); // reset trial after chain test: neuron time is set to zero; neuron parameters to initial values
        void reset_after_trial(); // continue the network activity by reset all neurons
        void setToRest_after_trial_soma_pre_dend_post_delays(); // set neurons to rest and clear all event arrays
        void setToRest_afterEpoch(); // set neurons to rest after epoch
        
        void reset_after_trial_soma_pre_dend_post_no_delays(); // resest network after trial with no axonal delays
        void reset_after_trial_soma_pre_dend_post_delays(); // resest network after trial with axonal delays
		void set_all_mature(); // set all neurons to a mature state

        // supporting STDP rules
        void LTD_burst(double &w, double t); // long-term depression for burst STDP rule
        void LTP_burst(double &w, double t); // long-term potentiation for burst STDP rule
        void LTD(double &w, double t); // long-term depression STDP rule
		void LTP(double &w, double t); // long-term potentiation STDP rule
		void LTP_toRescaled(double &w, double t); // long-term potentiation STDP rule for rescaled
		
		void potentiation_decay(); // apply potentiation decay to all RA-RA synapses
		void potentiation_decay_sudden_maturation(); // apply potentiation decay to all RA-RA synapses with immature targets.
													 // if target neuron is mature, active synapse doesn't decay 
		
		void update_synapse(int i, int j); // update synapse from neuron i to neuron j
		void update_synapse_sudden_maturation(int i, int j); // update synapse from neuron i to neuron j for the case of sudden maturation
		
		void update_all_synapses(); // update synapses between RA neurons potentiation decay
		void update_all_synapses_sudden_maturation(); // update synapses that contact immature targets 
		
		void axon_remodeling(int i); // remove all targets from neuron i except for supersynapses
		
		void rescale_synapses_to_mature(int neuron_id); // rescale all synapses to neuron that got mature
		void rescale_inhibition_to_mature(int neuron_id); // rescale inhibitory synapses to mature neurons
        
        // neurogenesis support
        //void add_new_neurons(int N); // add N immature neurons to network
        
        void replace_neurons(std::vector<int>& neurons_to_replace); // replace all neurons specified by replace arrays
        void remove_neurons_from_network(std::vector<int>& neurons_to_replace); // remove HVC-RA neurons from network
        void resample_neurons(std::vector<int>& neurons_to_replace); // resample coordinates and HVC-RA <-> HVC-I connections for replaced neurons
        void sample_coordinates_for_replaced(std::vector<int>& neurons_to_replace); // sample coordinates for replaced HVC-RA neurons
        void sample_connectionsAndDelays_for_replaced(std::vector<int>& neurons_to_replace); // sample connections and delays for replaced HVC-RA neurons
        void update_replaced_neurons(std::vector<int>& neurons_to_replace); // update properties of replaced neurons
        
        void find_chain_neurons(std::vector<int>& chain_neurons); // find neurons in the chain
        void make_lesion(double fraction); // kill fraction of neurons in the chain and replace them with new neurons
        void kill_neuron(int local_id, int global_id, int process_rank); // erase all outgoing and incoming connections from HVC(RA) 
                                                                         // neurons for replaced neuron. Clean indicator arrays, active and super synapses
	    
	    // create network topology
	    void sample_coordinates(); // sample coordinates of all neurons in the network
	    void sample_connections(); // sample connections between HVC-RA and HVC-I neurons
	    
	    double sample_G(double G_max); // sample strength of synaptic connection
	    void sample_axonal_delays(); // sample axonal delays based on distances between neurons and delay_constant 
									 // parameter in connection_parameters struct
	    
	    void set_delays_RA2RA(double delay); // set all delays between HVC-RA neurons to a fixed value
	    
	    double p_RA2I(double d); // probability of HVC-RA -> HVC-I  connections based on distances between neurons
	    double p_I2RA(double d); // probability of HVC-I  -> HVC-RA connections based on distances between neurons
	    
	    // read network topology and resample weights
	    void resample_weights(); // resample synaptic weights for HVC-RA -> HVC-I and HVC-I -> HVC-RA synapses
								 // based on parameters in connection_params structure
	    
	    
	    void create_second_layer(); // connect training neurons to all pool neurons
		void set_active_synapse_weights(double G); // set all active synapse weights to G
	    // coordinates and connections
	    
	    
	    
	    void initialize_coordinates(); // initialize coordinates of neurons
        
        //void initialize_coordinates_for_clustered_training(); // initialize coordinates so that first N_TR neurons are clustered
		//void initialize_coordinates_for_dispersed_training(); // initialize coordinates so that first N_TR neurons are dispersed

        void initialize_coordinates_for_replaced_neuron(int global_id); // change coordinates of replaced neuron to new coordinates
        void initialize_connections_for_replaced_neuron(int global_id); // erase all previous connections from and onto the replaced neuron. Create new connections.
        
        void initialize_coordinates_for_added_neurons(int n_total_old); // initialize coordinates for added neurons; 
                                                                    // n_total_old - total number of HVC(RA) neurons before new neurons are added;
        
        void initialize_connections_for_added_neurons(int n_total_old); // initialize connections for added neurons; 
                                                                    // n_total_old - total number of HVC(RA) neurons before new neurons are added;
        // spatial initialization
        void initialize_connections(); // initialize connections for neurons based on Gaussian distributions
        
        // testing functions
        void initialize_ideal_chain_connections(int num_layers); // initialize connections like in ideal synfire chain: previous chain layer makes connections on
                                                             // interneurons that in turn connect to the subsequent chain layer

        void initialize_random_chain_connections(int num_layers); // initialize connections like in a real synfire chain, but wire HVC(RA) neurons randomly, ignoring 
                                                                  // any inhibitory structure

		void set_noise(); // set noise for all neurons in the network. White noise is used for HVC(RA) and poisson noise for HVC(I)
		void set_dynamics(); // initialize vector arrays of proper size for dynamics of HVC(RA) and HVC(I) neurons
		
		
		
	
		
		void set_time_for_neurons(double t); // set initial time for all neurons in the network to t. Used when chain growth is restarted

      
        // internal printing functions
        void print_active(); // print active indicators
        void print_super(); // print super indicators

		// reading network topology
		void read_training_neurons(const char *filename); // read training neurons from file
		void read_training_spread(const char *filename); // read spread of training neurons from file
		void read_network_topology(std::string networkDir); // read network topology from files located in directory networkDir
		void read_number_of_neurons(const char *filename); // read number of HVC-RA and HVC-I neurons in the network
		void read_all_coordinates(const char *file_RA_xy, const char *file_I_xy); // read coordinates of HVC-RA and HVC-I neurons
		void read_connections_RAandI(const char* RA_I, const char* I_RA); // read synaptic connections between HVC-RA and HVC-I neurons

        // reading data from files
        void read_remodeled_indicators(const char *filename); // read axon-remodeling indicators for HVC-RA neurons
        void read_mature_indicators(const char *filename); // read maturation indicators for HVC-RA neurons
        void read_maturation_properties(const char *filename); // read maturation properties of HVC-RA neurons
        void read_super_synapses(const char* filename); // read supersynapses from the file
        void read_active_synapses(const char* filename); // read active synapses from the file
        //void read_maturation_info(const char* filename); // read maturation information from the file
        void read_weights(const char* filename); // read all synaptic weights from the file
		//void read_num_bursts_in_recent_trials(const char* filename); // read number of recent bursts produced by HVC(RA) neurons in RATE_WINDOW_LONG previous trials
		void read_replacement_history(const char* filename); // read number of trials since previous neuron replacements
		void read_activity_history(const char* filename); // read number of spikes of HVC-RA neurons produced in previous trials
		void read_axonal_delays(std::vector<std::vector<double>> &axonal_delays, const char* filename); // read axonal time delays between neurons
		
		//void read_last_dendritic_spike_times(const char* filename); // read last dendritic spike information to a file

        // internal write functions

        //void write_weights_time_sequence_from_source_to_target(const std::vector<int>& source, const std::vector<int>&target, const char* filename); // write time
        //void write_maturation_time_sequence(const std::vector<int>& neurons, const char* filename); // write maturation time sequence of neurons to a file
        //void write_replaced_neurons(std::vector<int>& real_id, const char* filename); // write replaced neurons to a file
		
        //void write_chain_test(int num_trials, std::vector<int>& total_num_dend_spikes, std::vector<double>& average_num_dendritic_spikes_in_trials, 
          //                    std::vector<double>& average_num_somatic_spikes_in_trials, std::vector<double>& mean_burst_time, 
            //                  std::vector<double>& std_burst_time, const char* filename); // write results of chain test to file
        
        void write_inhibition_tracking_state(int trial_number, double time_resolution_conductance, std::string outputDirectory); 
																			// write state of the network for simulation
																		   // which tracks inhibitory conductance
																		   // Info written each trial: inhibitory conductances of all neurons
																		   // num active and super synapses; spike times of all neurons
																		   
        void write_inhibitory_conductance(int trial, double time_conductance_resolution, std::string outputDir); // write inhibitory conductance of all HVC-RA neurons during the trial
        
        void write_num_synapses(const char* fileSynapses); // write amount of active synapses and supersynapses
        void write_active_synapses(const char* RA_RA); // write RA to RA active connections
        void write_super_synapses(const char* RA_RA); // write RA to RA supersynapses
        
        //void write_RA(const char * filename, int n); // write dynamics of RA neuron to a file
        //void write_I(const char * filename, int n); // write dynamics of I neuron to a file
        
        void write_activity_history(const char *filename); // write recent activity of HVC-RA neurons
        
        //void write_num_bursts_in_recent_trials(const char* filename); // write number of recent bursts produced by HVC(RA) neurons in RATE_WINDOW_LONG previous trials 
        void write_weights(const char * filename); // write weights of all synapses in network
        void write_weight_statistics(const char * filename); // write mean synaptic weight and synaptic weight standard deviation
        void write_soma_spike_times(const char* filename); // write somatic spike information to a file
        void write_dend_spike_times(const char* filename); // write dendritic spike information to a file
        //void write_last_dend_spike_times(const char* filename); // write last dendritic spike information to a file
        void write_interneuron_spike_times(const char* filename); // write interneuron spike information to a file
        //void write_maturation_info(const char* filename); // write mature neurons
        //void write_time_info(const char* filename); // write simulation time information
        void write_replacement_history(const char* filename); // write when HVC-RA neurons were replaced previous times 
        void write_remodeled_indicators(const char *filename); // write axon-remodeling indicators for HVC-RA neurons
        void write_mature_indicators(const char *filename); // write maturation indicators for HVC-RA neurons
		void write_maturation_properties(const char* filename); // write HVC-RA neuron maturation properties  
		void write_axonal_delays(const std::vector<std::vector<double>> &axonal_delays, const char *filename); // write axonal time delays between neurons
        
        void set_recording(const std::vector<int> &RA_neurons, const std::vector<int> &I_neurons, int trial_number,
																std::string outputDirectory); // record dynamics of HVC-RA and HVC-I neurons to files
																
        void write_number_of_neurons(const char* filename); // write number of HVC-RA and HVC-I neurons in network
        void write_training_neurons(const char *filename); // write training neurons to a file
        void write_training_spread(const char* filename); // write spread of training neurons to a file
        
        void write_coordinates(const vector<double> &x, const vector<double> &y, const vector<double> &z,
										const char *filename); // write coordinates to a file
		
		void write_connections_RAandI(const char *file_RA_I, const char *file_I_RA); // write connections between HVC-RA
																				// and HVC-I neurons to a file
																				
		void write_network_topology(std::string outputDir); // write network topology to a directory
		
		void write_different_training_arrangements(std::string networkDir); // write two different sets of training neurons
																			// 1st is clustered and 2nd is randomly chosen
		
		void write_pajek_topology(const char* filename); // write network topology to a pajek .net file
		
		
		void write_alterable_network(std::string fileEnding, std::string outputDirectory); // write coordinates of HVC-RA neurons and
																		// and connections between HVC-RA and HVC-I neurons to directory outputDirectory
																		// extension will be added to filenames to distinguish them
        
        void write_full_network_state(std::string fileEnding, std::string outputDirectory); // write full state of growth simulation to a directory
        void write_graph_network_state(std::string outputDirectory); // write state of growth simulation needed for gui application to a directory
		void write_chain_test(std::string outputDirectory, int trial_number); // write somatic, dendritic and interneuron spikes to a difrectory
		
		void write_inhAndExc(int num_trials, 
								const std::vector<int> &num_trials_with_relevant_spikes, 
								const std::vector<double> &mean_relevant_spike_time, 
								const std::vector<double> &std_relevant_spike_time, 
								const std::vector<int> &num_trials_with_nonrelevant_spikes, 
								const std::vector<double> &mean_spike_time,
								const std::vector<double> &std_spike_time,
								const char *filename); // write results of inhibitory and excitatory input test
														// to a file
		
		void write_jitter(int num_trials,
							const std::vector<int> &num_trials_with_soma_spikes, 
							const std::vector<double> &average_num_soma_spikes_in_trial,
							const std::vector<double> &mean_first_soma_spike_time,
							const std::vector<double> &std_first_soma_spike_time,
							const std::vector<int> &num_trials_with_dend_spikes, 
							const std::vector<double> &average_num_dend_spikes_in_trial,
							const std::vector<double> &mean_first_dend_spike_time,
							const std::vector<double> &std_first_dend_spike_time,
							const char *filename); // write results of jitter test to a file
		
		void write_replaced_neurons(const std::vector<int>& replaced_neurons, 
												const char* filename); // write neurons that were replaced to a file
		//////////////////////////
		// Check correctness
		//////////////////////////
		int check_bad_values(); // check if neuron membrane potentials are in reasonable range
		
		void disable_RA2I_immature_outputs(); // disable output HVC-RA -> HVC-I connections for immature neurons
        
        ///////////////////////////////
        // Sending and receiving data
        ///////////////////////////////
        void set_synapse_indicators(); // set active and super synapse indicators
        
        void initialize_global_inhibitory_conductances(double time_resolution_conductance); // initialize global array on master process
																							// for inhibitory conductance of HVC-RA neurons
        
        void initialize_local_inhibitory_conductances(double time_resolution_conductance); // initialize local arrays on all processes
																							// for inhibitory conductance of HVC-RA neurons
        
		void resize_arrays_for_master_process(); // resize global data arrays for master process
        void resize_arrays_for_all_processes(); // resize local data arrays for all processes
        
		void send_growth_parameters(); // send parameters of growth simulation to all processes
	    void send_connections_RAandI(); // send connections between HVC-RA and HVC-I neurons to all processes
        void send_connections_RA2RA(); // send connections between HVC-RA neurons to all processes
        void send_axonal_delays_RA2RA(); // send axonal delays between HVC-RA neurons to all processes
        
        void send_replacement_history(); // send replacement history from master to all processes 
        void send_activity_history(); // send activity history from master to all processes 
        void send_active_synapses(); // send active synapses from master to all processes
        void send_super_synapses(); // send super synapses from master to all processes
        void send_remodeled_indicators(); // send remodeled indicators from master to all processes
        void send_mature_indicators(); // send maturation indicators from master to all processes
        void send_maturation_properties(); // send maturation properties from master to all processes
        
        
        void gather_localVectorsToGlobal(const std::vector<int>& vector_local, 
													std::vector<int>& vector_global); // gather elements of local arrays from processes to global array on master process
													
        
        void gather_inhibition(); // gather inhibitory conductance of HVC-RA neurons from local arrays to global on master process
        void gather_mature_data(std::vector<std::vector<double>>& average_dendritic_spike_time); // gather data from all processes in case of mature chain trial
	    void gather_graph_state_data(); // gather data needed for GUI from all processes
	    void gather_full_state_data(); // gather additional data needed for full network state from all processes
        void gather_neurons_2replace(); // gather neurons that are to be replaced from all processes
        void gather_bursts(vector<int>& RA_fired_dend_real_ID, vector<int>& RA_fired_dend_global,
                           vector<double>& spike_times_fired_dend_local); // gather dendritic bursts in the network

		void gather_spiked_or_bursted_neurons(const std::vector<int>& RA_neurons_local, 
									std::vector<int>& RA_neurons_global); // gather id of local HVC-RA neurons to all processes that spiked or bursted

        int MPI_size; // number of processes
        int MPI_rank; // rank of the process

		vector<int> N_RA_sizes; // array with number of RA neurons in each process
		vector<int> N_I_sizes; // array with nubmer of I neurons in each process

	void get_neuronRA_location(int n, int* rank, int* shift); // get location of RA neuron with ID n in terms of process and position in array
	void get_neuronI_location(int n, int* rank, int *shift); // get location of I neuron with ID n in terms of process and position in array



};

#endif
