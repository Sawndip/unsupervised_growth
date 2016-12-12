#ifndef POOLPARALLEL_H_INCLUDED
#define POOLPARALLEL_H_INCLUDED

#include <vector>
#include <deque>
#include <cmath>
#include <mpi.h>
#include "poisson_noise.h"
#include <boost/circular_buffer.hpp>

class HH2_final_pool;
class HHI_final_pool;

using std::vector;

typedef boost::circular_buffer<int> intBuffer;

class PoolParallel
{
public:
	PoolParallel(double a, double s_rai, double b, double sigma_ira, double network_update, double Ei,
				 double beta, double beta_s, double Tp, double Td, double tauP, double tauD, double Ap,
				 double Ad, double Ap_super, double Ad_super, double f0, double activation, double super_threshold, 
				 double gaba_down,  double Gmax, int N_ra, int Ni, int N_ss, int N_tr);
	
	~PoolParallel();

	void read_from_file(const char* RA_xy, const char* I_xy, const char* RA_RA_all, const char* RA_RA_active,
						const char* RA_RA_super, const char* RA_I, const char* I_RA, const char* maturation,
						const char* timeInfo); // read network structure from files
	
	void read_connections_from_net(const char *filename, std::vector<int>** target_id, std::vector<double>** target_G); // read connection from .net file
	void read_all_connections_from_net(const char *filename, double*** weights); // read connection from .net file
    void read_all_connections(const char *filename, double** weights); // read all connections between RA neurons from binary file

	void initialize_coordinates(); // initialize coordinates of neurons
	void initialize_generator(); // initialize generator for processes

    void initialize_connections(double Gei_mean, double Gei_var, double Gie_mean, double Gie_var); // initialize connections for neurons

    int get_trial_number(); // get current number of trials performed

	void send_connections(); // send fixed connections connections to all processes
	void send_simulation_parameters(); // send parameters of simulation and RA to RA connections to all processes
	
	void update_synaptic_info(); // update synaptic information after reading network structure from files

	void mature_chain_test(const int num_trials, const char* file_soma_spikes, const char* file_dend_spikes, const char* file_chain_test); // test of mature network
	void trial(int training); // make one trial

	void randomize_after_trial(); // set all neurons to the resting state
	void reset_after_trial(); // reset neurons after trial is over
	void initialize_pool(double Gie_max, double Gei_max); // initialize pools with synapses stored outside neurons
	void visualize_RA_connections(); // visualize RA connections of the network
	void visualize_I_connections(); // visualize I connections of the network
    void print_invariable_connections(); // get number of connections for invariable synapses
	// write to file functions
	void gather_data(); // gather data from all processes

	void write_sim_info(const char* simInfoFile, int synapses_trials_update, int weights_trials_update); // write information about simulation and frequency of updates of synapses and weights
	void write_num_synapses(const char* fileSynapses); // write amount of active synapses and supersynapses
	void write_active_synapses(const char* RA_RA); // write RA to RA active connections
	void write_supersynapses(const char* RA_RA); // write RA to RA supersynapses
	void write_invariable_synapses(const char* RA_I, const char* I_RA); // write RA to I and I to RA connections
	void write_coordinates(const char* xy_RA, const char* xy_I); // write coordinates of neurons
	void write_RA(const char * filename, int n); // write dynamics of RA neuron to a file
	void write_I(const char * filename, int n); // write dynamics of I neuron to a file
	void write_weights(const char * filename); // write weights of all synapses in network
    void write_soma_time_info(const char* filename); // write somatic spike information to a file
    void write_dend_time_info(const char* filename); // write dendritic spike information to a file
    void write_interneuron_time_info(const char* filename); // write interneuron spike information to a file
	void write_maturation_info(const char* filename); // write mature neurons
    void write_time_info(const char* filename); // write simulation time information

    void write_pajek_super(const char* filename); // write supersynapses to a file for pajek
	void write_pajek_active(const char* filename); // write active synapses to a file for pajek
	void write_pajek_all(const char* filename); // write all synapses to a file for pajek
    void write_pajek_fixed(const char* filename); // write fixed synapses to a file for pajek

	// set functions
	void set_simulation_parameters(double interval, double tS, double mu_s, double sigma_s, double mu_d, double sigma_d); // set simulation parameters: trial duration; timestep of dynamics; white noise input to neuron
   
   // current
    void set_training_current(); // set current to training neurons

	void print_simulation_parameters(); // print simulation parameters
protected:
		// number of neurons
		int N_TR; // number of training HVC(RA) neurons
		int N_RA; // number of HVC(RA) neurons
		int N_I; // number of HVC(I) neurons
		int N_RA_local; // number of RA neurons in each process
		int N_I_local; // number of I neurons in each process
		int Nss; // number of super synapses

		HH2_final_pool* HVCRA_local; // array of HVC(RA) neurons
		HHI_final_pool* HVCI_local; // array of HVC(I) neurons

		// coordinates info
		vector <double> xx_RA; // x-coordinates of RA neurons
		vector <double> yy_RA; // y-coordinates of RA neurons
		vector <double> xx_I; // x-coordinates of I neurons
		vector <double> yy_I; // y-coordinates of I neurons
		vector <double> xx_I_cluster; // x-coordinates of centers of inhibitory clusters
		vector <double> yy_I_cluster; // y-coordinates of centers of inhibitory clusters

		const static double CLUSTER_SIZE; // cluster size
        const static double MIN_INTERNEURON_DISTANCE; // minimum distance between neurons
		double A_RA2I; // constant for nearby connections
		double SIGMA_RA2I; // variance of Gaussian distribution for far away connections

		double SIGMA_I2RA; // spatial scale of probability of connections decay
		double B_I2RA; // constant for I to RA connections
		const static double SIDE; // length of HVC side
		
		const static double DEMOTIVATION_WINDOW; // demotivation window in ms
		const static double DELAY_WINDOW; // delay window in which spikes are assumed to be relevant to inhibitory kicks
		double mu_soma; // dc component of white noise to soma
		double sigma_soma; // variance of white noise to soma
		double mu_dend; // dc component of white noise to dendrite
		double sigma_dend; // variance of white noise to dendrite

		int size; // time array size for dynamics
		double timeStep; // time step for dynamics
		double trial_duration; // duration of trial in milliseconds
		double internal_time = 0; // internal time of each neuron
        double network_time = 0; // network time
        double network_update_frequency; // how often network state should be updated in ms

		double current_injection_time; // injection time of training current
		int trial_number = 0; // number of simulated trials

		double delay; // delay between inhibitory kick and dendritic spike
		int num_related_spikes; // number of spikes assumed to be related to the inhibitory kick

		Poisson_noise generator; // poisson noise generator

		bool** active_global; // array of HVC(RA) neurons with active synapses
		bool** supersynapses_global; // indicator array for active supersynapses;
		double** weights_global; // array of synaptic strength of all connections between RA neurons
		vector<double>* weights_RA_I_global; // array with synapses from RA to I neurons
		vector<double>* weights_I_RA_global; // array with senapses from I to RA neurons
		vector<int>* syn_ID_RA_I_global; // array with synaptic ID numbers from RA to I neurons
		vector<int>* syn_ID_I_RA_global; // array with synaptic ID numbers from I to RA neurons

		vector<int>* active_synapses_global; // array of vectors with IDs of active synapses
		vector<int>* active_supersynapses_global; // array of vectors with IDs of supersynapses


        std::deque<double>* last_soma_spikes_local; // last NUM_SOMA_SPIKES spikes in somatic compartment
        std::deque<double>* last_soma_spikes_global; // last NUM_SOMA_SPIKES spikes in somatic compartment

        double* spike_times_soma_local; // local array with neurons somatic spike times
		double* spike_times_dend_local; // local array with neurons dendritic spike times

		double* spike_times_soma_global; // array with the most recent somatic spike times of neurons
		double* spike_times_dend_global; // array with the most recent dendritic spike times of neurons

        vector<double>* spikes_in_trial_soma_global; // array with somatic spike times in the last trial
		vector<double>* spikes_in_trial_dend_global; // array with dendritic spike times in the last trial
		vector<double>* spikes_in_trial_interneuron_global; // array with interneuron spike times in the last trial

        vector<double>* spikes_in_trial_soma_local; // array with somatic spike times in the last trial
		vector<double>* spikes_in_trial_dend_local; // array with dendritic spike times in the last trial
		vector<double>* spikes_in_trial_interneuron_local; // array with interneuron spike times in the last trial

		bool** active_local; // array of HVC(RA) neurons with active synapses
		bool** supersynapses_local; // indicator array for active supersynapses;
		double** weights_local; // array of synaptic strength of all connections between RA neurons
		
		vector<bool> remodeled_local; // indicators if neuron underwent axon remodelling
		
		vector<int> maturation_triggered_local; // local indicators if neuron maturation process is triggered
		vector<int> maturation_triggered_global; // global indicators if neuron maturation process is triggered
		
		vector<int> mature_local; // indicator if neuron is mature due to supersynaptic acquisition (neuron matures forever after it acquires all supersynapses)
		vector<int> mature_global; // global array of indicators if neuron is mature due to supersynaptic acquisition 
		
		vector<double> gaba_potential_local; // array with local values of GABA reverse potential
		vector<double> gaba_potential_global; // array with global values of GABA reverse potential

		vector<intBuffer> num_bursts_in_recent_trials; // array of circular buffer containing recent neuron rates

		vector<double>* spikes_in_trial_local; // local rray with spike times of single trial
		vector<double>* weights_RA_I_local; // array with synapses from RA to I neurons
		vector<double>* weights_I_RA_local; // array with senapses from I to RA neurons
		vector<int>* syn_ID_RA_I_local; // array with synaptic ID numbers from RA to I neurons
		vector<int>* syn_ID_I_RA_local; // array with synaptic ID numbers from I to RA neurons
		
		vector<int>* active_synapses_local; // array of vectors with IDs of active synapses
		vector<int>* active_supersynapses_local; // array of vectors with IDs of supersynapses
		vector<int> Id_RA_local; // Id of RA neurons in each process
		vector<int> Id_I_local; // Id of I neurons in each process

		// update conductances and glutamate
		vector<double> update_Ge_RA_local;
		vector<double> update_Gi_RA_local;
		vector<double> update_Ge_I_local;
	
		vector<double> update_Ge_RA_global;
		vector<double> update_Gi_RA_global;
		vector<double> update_Ge_I_global;

		//const static double ACTIVATION; // activation threshold for synapses
		double p_RA2I(int i, int j); // probability of connection from RA to I neuron
		double p_I2RA(int i, int j); // probability of connection from I to RA neuron

		//const static double SUPERSYNAPSE_THRESHOLD; // threshold for supersynaptic connection

		// developmental GABA switch
		void update_Ei(); // update GABA reverse potential after trial
		double E_GABA_IMMATURE; // immature reverse GABA potential
		const static double E_GABA_MATURE; // mature reverse GABA potential
        
		const static int RATE_WINDOW; // window in which neuron firing rate is estimated
		const static double BURST_RATE_THRESHOLD; // firing rate threshold for D-H shift
		double GABA_DOWN; // decrease value for gaba potential

		// constants for STDP-rules
		const static double STDP_WINDOW; // window for STDP beyond which we ignore all spikes

        double ACTIVATION; // threshold for synapse activation
        double SUPERSYNAPSE_THRESHOLD; // threshold for supersynapse activation
        double BETA; // potentiation decay parameter
        double BETA_SUPERSYNAPSE; // potentiation decay parameter for supersynapse
        double G_MAX; // maximum synaptic weight
        double A_D; // LTD amplitude 
        double A_P; // LTP amplitude
		double A_P_SUPER; // LTP amplitude for supersynapse
		double A_D_SUPER; // LTD amplitude for supersynapse
		double T_P; // time for the most efficient LTP
		double TAU_P; // LTP decay time
		double T_D; // time for the most efficient LTD
		double TAU_D; // LTD decay


		//const static double A_P;
		const static double G_P;

		//const static double A_D;

		const static double R; // learning rate
		double F_0; // constant to prevent connections within the same chain group

		void set_training_current(double t); // set current to training neurons. t - current injection time.
		
		void write_chain_test(int num_trials, std::vector<int>& num_dend_spikes, std::vector<double>& mean_burst_time, std::vector<double>& std_burst_time, const char* filename); // write results of chain test to file
		
		void mature_trial(); // simulation trial without STDP rules
		void LTD(double &w, double t); // long-term depression STDP rule
		void LTP(double &w, double t); // long-term potentiation STDP rule
		void potentiation_decay(); // apply potentiation decay to all RA-RA synapses
		void update_synapse(int i, int j); // update synapse from neuron i to neuron j
		void update_all_synapses(); // update synapses between RA neurons and apply potentiation decay
		void axon_remodeling(int i); // remove all targets from neuron i except for supersynapses

        // MPI support
		void gather_mature_data(std::vector<std::vector<double>>& average_dendritic_spike_time); // gather data from all processes in case of mature chain trial

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
