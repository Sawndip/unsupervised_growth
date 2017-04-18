#ifndef POOLPARALLEL_H_INCLUDED
#define POOLPARALLEL_H_INCLUDED

#include <vector>
#include <cmath>
#include <mpi.h>
#include "poisson_noise.h"
#include <boost/circular_buffer.hpp>
#include <boost/dynamic_bitset.hpp>
#include "Configuration.h"
#include "HHI_final.h"
#include "HH2_final.h"


using std::vector;

typedef boost::circular_buffer<int> intBuffer;
typedef boost::dynamic_bitset<> bitArray;

class PoolParallel
{
public:
	PoolParallel(const Configuration& cfg);
	
    void run_trials_no_save(int num_trials); // run num_trials trials without saving data to files
    void run_trials_with_save(int num_trials); // run num_trials trials with saving sata to files
    
    void read_network_state(std::string dirname); // read network state from files in the directory dirname



    void initialize_test_connections(int num_RA_targets, int num_RA_target_groups); // initialize test connections: first training neuron is connected to the first
                                                                                    // interneuron; First interneuron in turn is connected to num_RA_targets RA
                                                                                    // neurons. Strength if inhibitory connection increases by Gie_mean in each
                                                                                    // next group of RA neurons
	
	void chain_growth_with_no_RA2I_connections(int save_freq_short, int save_freq_long); // run chain growth with disable connections from HVC(RA) to HVC(I)
																						 // neurons. It tests if noisy interneurons can result in neuron recruitment
	
    void chain_growth_default(bool training, int save_freq_short, int save_freq_long); // run chain growth algorithm with default connections and coordinates initialization;
                                                                        // save data for graph update every save_freq_short trials; 
                                                                        // data for analysis every save_freq_long trials

    void chain_growth_with_clustered_training(bool training, int save_freq_short, int save_freq_long); // run chain growth algorithm with default connections and coordinates 
                                                                        // for clustered training neurons;
                                                                        // save data for graph update every save_freq_short trials; 
                                                                        // data for analysis every save_freq_long trials
    
    void chain_growth_manual(bool training, int save_freq_short, int save_freq_long); // run chain growth with manually specified coordinates and connections. 
                                                                       // NOTE: coordinates and connections MUST be initialized before using chain_growth_manual

    void test_grown_chain(int num_trials, std::string dataDir, std::string outputDir); // test grown synfire chain. All data files should 
                                                                            // be located in directory dataDir; output goes to directory outputDir

    void test_random_chain(int num_layers, int num_trials); // test chain with ideal synfire chain connections. Neurons in the chain are wired ignoring
                                                            // any inhibibory structure
    void test_ideal_chain(int num_layers, int num_trials); // run testing trial with ideal chain connections



    void print_invariable_connections(); // get number of connections for invariable synapses
    void print_received_invariable_connections(); // show received invariable connections

	void print_simulation_parameters(); // print simulation parameters
protected:
        std::string outputDirectory; // directory to which write output
		
        // number of neurons
		int N_TR; // number of training HVC(RA) neurons
		int N_RA; // number of HVC(RA) neurons
		int N_I; // number of HVC(I) neurons
		int N_RA_local; // number of RA neurons in each process
		int N_I_local; // number of I neurons in each process
		int Nss; // number of super synapses

		vector<HH2_final> HVCRA_local; // array of HVC(RA) neurons
		vector<HHI_final> HVCI_local; // array of HVC(I) neurons

		// coordinates info
		vector <double> xx_RA; // x-coordinates of RA neurons
		vector <double> yy_RA; // y-coordinates of RA neurons
		vector <double> xx_I; // x-coordinates of I neurons
		vector <double> yy_I; // y-coordinates of I neurons
		vector <double> xx_I_cluster; // x-coordinates of centers of inhibitory clusters
		vector <double> yy_I_cluster; // y-coordinates of centers of inhibitory clusters

        double Gei_mean; // mean synaptic weight from HVC(RA) to HVC(I) neuron
        double Gei_std; // standard deviation of synaptic weight from HVC(RA) to HVC(I) neuron
        double Gie_mean; // mean synaptic weight from HVC(I) to HVC(RA) neuron
        double Gie_std; // standard deviation of synaptic weight from HVC(I) to HVC(RA) neuron

        double MIN_INTERNEURON_DISTANCE; // minimum distance between neurons
		
        double A_RA2I; // constant for nearby connections
		double SIGMA_RA2I; // variance of Gaussian distribution for far away connections

		double SIGMA_I2RA; // spatial scale of probability of connections decay
		double B_I2RA; // constant for I to RA connections
		double SIDE; // length of HVC side
		
		double white_noise_mean_soma; // dc component of white noise to soma
		double white_noise_std_soma; // variance of white noise to soma
		double white_noise_mean_dend; // dc component of white noise to dendrite
		double white_noise_std_dend; // variance of white noise to dendrite

		int size; // time array size for dynamics
		double timeStep; // time step for dynamics
		double trial_duration; // duration of trial in milliseconds
		double internal_time; // internal time of each neuron
        double network_time; // network time
        double network_update_frequency; // how often network state should be updated in ms

		double current_injection_time; // injection time of training current
		int trial_number; // number of simulated trials

		Poisson_noise generator; // poisson noise generator

		vector<vector<double>> weights_global; // array of synaptic strength of all connections between RA neurons
        vector<vector<double>> weights_local; // array of synaptic strength of all connections between RA neurons
		
        vector<vector<double>> weights_RA_I_global; // array with synapses from RA to I neurons
		vector<vector<double>> weights_I_RA_global; // array with senapses from I to RA neurons
		vector<vector<double>> weights_RA_I_local; // array with synapses from RA to I neurons
		vector<vector<double>> weights_I_RA_local; // array with senapses from I to RA neurons
        
        vector<vector<int>> syn_ID_RA_I_global; // array with synaptic ID numbers from RA to I neurons
		vector<vector<int>> syn_ID_I_RA_global; // array with synaptic ID numbers from I to RA neurons
		vector<vector<int>> syn_ID_RA_I_local; // array with synaptic ID numbers from RA to I neurons
		vector<vector<int>> syn_ID_I_RA_local; // array with synaptic ID numbers from I to RA neurons

		vector<bitArray> active_indicators_local; // array of HVC(RA) neurons with active synapses
		vector<bitArray> supersynapses_indicators_local; // indicator array for active supersynapses;
		
		vector<vector<int>> active_synapses_local; // array of vectors with IDs of active synapses
		vector<vector<int>> supersynapses_local; // array of vectors with IDs of supersynapses
		
        vector<vector<int>> active_synapses_global; // array of vectors with IDs of active synapses
		vector<vector<int>> supersynapses_global; // array of vectors with IDs of supersynapses
		
		vector<double> spike_times_dend_global; // array with the most recent dendritic spike times of neurons

        vector<vector<double>> spikes_in_trial_soma_global; // array with somatic spike times in the last trial
		vector<vector<double>> spikes_in_trial_dend_global; // array with dendritic spike times in the last trial
		vector<vector<double>> spikes_in_trial_interneuron_global; // array with interneuron spike times in the last trial

        vector<vector<double>> spikes_in_trial_soma_local; // array with somatic spike times in the last trial
		vector<vector<double>> spikes_in_trial_dend_local; // array with dendritic spike times in the last trial
		vector<vector<double>> spikes_in_trial_interneuron_local; // array with interneuron spike times in the last trial
		
		vector<int> remodeled_local; // local indicators if neuron underwent axon remodelling
		vector<int> remodeled_global; // global indicators if neuron underwent axon remodelling
		
		vector<int> mature_local; // indicator if neuron is mature due to supersynaptic acquisition (neuron matures forever after it acquires all supersynapses)
		vector<int> mature_global; // global array of indicators if neuron is mature due to supersynaptic acquisition 
		
		vector<double> gaba_potential_local; // array with local values of GABA reverse potential
		vector<double> gaba_potential_global; // array with global values of GABA reverse potential

		vector<double> firing_rate_local; // array with local firing rates of HVC(RA) neurons
		vector<double> firing_rate_global; // array with global firing rates of HVC(RA) neurons
    
        vector<int> gaba_reached_mature_local; // local indicators for neurons that GABA potential reached mature value

        vector<intBuffer> num_bursts_in_recent_trials; // array of circular buffer containing recent neuron rates

		vector<int> Id_RA_local; // local array with id of RA neurons in each process
		vector<int> Id_RA_global; // global array with ids of RA neurons in all processes
		vector<int> Id_I_local; // Id of I neurons in each process

		// update conductances and glutamate
		vector<double> update_Ge_RA_local;
		vector<double> update_Gi_RA_local;
		vector<double> update_Ge_I_local;
	
		vector<double> update_Ge_RA_global;
		vector<double> update_Gi_RA_global;
		vector<double> update_Ge_I_global;

        // neurons to be replaced in neurogenesis process
        vector<int> replace_local_id_local; // local id of neurons to be replaced stored locally
        vector<int> replace_real_id_local; // real id of neurons to be replaced stored locally
        vector<int> replace_process_rank_local; // process rank that contains neuron to be replace
        
        vector<int> replace_local_id_global; // all local id of neurons to be replaced
        vector<int> replace_real_id_global; // all real id of neurons to be replaced
        vector<int> replace_process_rank_global; // all process ranks that contain neurons to be replaced

        // neuron's internal time after previous replacement
        vector<int> num_trials_after_replacement_local;

		//const static double ACTIVATION; // activation threshold for synapses
		double p_RA2I(int i, int j); // probability of connection from RA to I neuron
		double p_I2RA(int i, int j); // probability of connection from I to RA neuron

        double sample_Ge2i(); // sample synaptic weight from HVC(RA) to HVC(I) neuron
        double sample_Gi2e(); // sample synaptic weight from HVC(I) to HVC(RA) neuron

		//const static double SUPERSYNAPSE_THRESHOLD; // threshold for supersynaptic connection

		// developmental GABA switch
		void update_Ei(); // update GABA reverse potential after trial
		
        double E_GABA_IMMATURE; // immature reverse GABA potential
		double E_GABA_MATURE; // mature reverse GABA potential
        
        double GABA_RATE_THRESHOLD; // rate threshold for GABA decrease
        double MATURATION_RATE_THRESHOLD; // rate threshold for neuron maturation
        double DEATH_RATE_THRESHOLD; // rate threshold for neuron death

        int RATE_WINDOW_SHORT; // window in which rate for gaba decrease is calculated
        int RATE_WINDOW_LONG; // window in which rate for neuron maturation and death is calculated
		
        double GABA_DOWN; // decrease value for gaba potential

		// constants for STDP-rules
		double STDP_WINDOW; // window for STDP beyond which we ignore all spikes

        double ACTIVATION_THRESHOLD; // threshold for synapse activation
        double SUPERSYNAPSE_THRESHOLD; // threshold for supersynapse activation
        double BETA; // potentiation decay parameter
        double BETA_SUPERSYNAPSE; // potentiation decay parameter for supersynapse
        double WEIGHT_MAX; // maximum synaptic weight
        double A_D; // LTD amplitude 
        double A_P; // LTP amplitude
		double T_P; // time for the most efficient LTP
		double TAU_P; // LTP decay time
		double T_D; // time for the most efficient LTD
		double TAU_D; // LTD decay

		double R; // learning rate
		double F_0; // constant to prevent connections within the same chain group

        double T_P1; // end of silent LTP window
        double T_P2; // time of the best LTP response

        // current
        double WAITING_TIME; // time before training current is injected
		void set_training_current(double t); // set current to training neurons. t - current injection time.
        void set_training_current(); // set current to training neurons
		
	    // noise generator	
	    void initialize_generator(); // initialize generator for processes
        
        // trials
	    void trial(int training); // make one trial
	    void trial_burst_stdp(int training); // make one trial with STDP rules applied to dendritic bursts
		void mature_trial(); // simulation trial without STDP rules
        void test_mature_chain(int num_trials); // test of mature network
	    
        void randomize_after_trial(); // set all neurons to the resting state
		void set_all_mature(); // set all neurons to a mature state

        // supporting STDP rules
        void LTD_burst(double &w, double t); // long-term depression for burst STDP rule
        void LTP_burst(double &w, double t); // long-term potentiation for burst STDP rule
        void LTD(double &w, double t); // long-term depression STDP rule
		void LTP(double &w, double t); // long-term potentiation STDP rule
		void potentiation_decay(); // apply potentiation decay to all RA-RA synapses
		void update_synapse(int i, int j); // update synapse from neuron i to neuron j
		void update_all_synapses(); // update synapses between RA neurons and apply potentiation decay
		void axon_remodeling(int i); // remove all targets from neuron i except for supersynapses

        // neurogenesis support
        void add_new_neurons(int N); // add N immature neurons to network
        void replace_neurons(); // replace all neurons specified by replace arrays
        void kill_neuron(int local_id, int global_id, int process_rank); // erase all outgoing and incoming connections from HVC(RA) 
                                                                         // neurons for replaced neuron. Clean indicator arrays, active and super synapses
	    // coordinates and connections
	    void disable_RA2I_connections(); // disable all connections from HVC(RA) to HVC(I) neurons
	    
	    void initialize_coordinates(); // initialize coordinates of neurons
        
        void initialize_coordinates_for_clustered_training(); // initialize coordinates so that first N_TR neurons are clustered

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

        // setting simulation parameters
        void set_network_parameters(const struct NetworkParameters& network_params); // set network parameters
        void set_noise_parameters(const struct NoiseParameters& noise_params); // set noise parameters
        void set_synaptic_parameters(const struct SynapticParameters& syn_params); // set synaptic parameters
        void set_gaba_parameters(const struct GabaParameters& gaba_params); // set GABA parameters of maturation
        void set_spatial_parameters(const struct SpatialParameters& spatial_params); // set spatial parameters
        void set_time_parameters(const struct TimeParameters& time_params); // set time parameters
    
        // internal printing functions
        void print_active(); // print active indicators
        void print_super(); // print super indicators

        // reading data from files
        void read_global_index_array(const char* filename); // read global array with HVC(RA) neuronal ids
        void read_invariable_synapses(const char* fileRA2I, const char* fileI2RA); // read HVC(RA) -> HVC(I) and HVC(I) -> HVC(RA)
                                                                                                     // connections from the file
        void read_super_synapses(const char* filename); // read supersynapses from the file
        void read_active_synapses(const char* filename); // read active synapses from the file
        void read_maturation_info(const char* filename); // read maturation information from the file
        void read_weights(const char* filename); // read all synaptic weights from the file

        // internal write functions

        void write_global_index_array(const char* filename); // write global array with HVC(RA) neuronal ids to file
        void write_weights_time_sequence_from_source_to_target(const std::vector<int>& source, const std::vector<int>&target, const char* filename); // write time
        void write_maturation_time_sequence(const std::vector<int>& neurons, const char* filename); // write maturation time sequence of neurons to a file
        void write_replaced_neurons(std::vector<int>& real_id, const char* filename); // write replaced neurons to a file
		
        void write_chain_test(int num_trials, std::vector<int>& total_num_dend_spikes, std::vector<double>& average_num_dendritic_spikes_in_trials, 
                              std::vector<double>& average_num_somatic_spikes_in_trials, std::vector<double>& mean_burst_time, 
                              std::vector<double>& std_burst_time, const char* filename); // write results of chain test to file
        
        void write_num_synapses(const char* fileSynapses); // write amount of active synapses and supersynapses
        void write_active_synapses(const char* RA_RA); // write RA to RA active connections
        void write_supersynapses(const char* RA_RA); // write RA to RA supersynapses
        void write_invariable_synapses(const char* RA_I, const char* I_RA); // write RA to I and I to RA connections
        void write_all_coordinates(); // write coordinates of HVC(RA) and HVC(I) neurons
        void write_coordinates_RA(const char* filename); // write coordinates of HVC(RA) neurons
        void write_coordinates_I(const char* filename); // write coordinates of HVC(I) neurons
        void write_RA(const char * filename, int n); // write dynamics of RA neuron to a file
        void write_I(const char * filename, int n); // write dynamics of I neuron to a file
        void write_weights(const char * filename); // write weights of all synapses in network
        void write_weight_statistics(const char * filename); // write mean synaptic weight and synaptic weight standard deviation
        void write_soma_spike_times(const char* filename); // write somatic spike information to a file
        void write_dend_spike_times(const char* filename); // write dendritic spike information to a file
        void write_interneuron_spike_times(const char* filename); // write interneuron spike information to a file
        void write_maturation_info(const char* filename); // write mature neurons
        void write_time_info(const char* filename); // write simulation time information

        void write_pajek_fixed(const char* filename); // write fixed synapses to a file for pajek

        // MPI support
		void resize_arrays_for_I(int n_local, int n_total); // resize data arrays for HVC(I) neurons
        void resize_arrays_for_RA(int n_local, int n_total); // resize data arrays for HVC(RA) neurons
        
	    void send_connections(); // send fixed connections connections to all processes
        
        void gather_mature_data(std::vector<std::vector<double>>& average_dendritic_spike_time); // gather data from all processes in case of mature chain trial
	    void gather_data(); // gather data from all processes
        void gather_neurons_2replace(); // gather neurons that are to be replaced from all processes
        void gather_bursts(vector<int>& RA_fired_dend_real_ID, vector<int>& RA_fired_dend_global,
                           vector<double>& spike_times_fired_dend_local); // gather dendritic bursts in the network

        int MPI_size; // number of processes
        int MPI_rank; // rank of the process

		vector<int> N_RA_sizes; // array with number of RA neurons in each process
		vector<int> N_I_sizes; // array with nubmer of I neurons in each process

	void get_neuronRA_location(int n, int* rank, int* shift); // get location of RA neuron with ID n in terms of process and position in array
	void get_neuronI_location(int n, int* rank, int *shift); // get location of I neuron with ID n in terms of process and position in array



};

static double distance(double x1, double y1, double x2, double y2)
{
	double d;
	d = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	return d;
}



#endif // POOLPARALLEL_H_INCLUDED
