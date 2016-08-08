#ifndef POOLPARALLEL_H_INCLUDED
#define POOLPARALLEL_H_INCLUDED

#include <vector>
#include <deque>
#include <cmath>
#include <mpi.h>
#include "poisson_noise.h"


class HH2_final_pool;
class HHI_final_pool;

using std::vector;

//typedef std::function<double (double)> DDfunction;

class PoolParallel
{
public:
	PoolParallel(int N_tr, int N_ra, int Nic, int NiInC, int Ns);
	~PoolParallel();

	void read_from_file(const char* RA_xy, const char* I_xy, const char* RA_RA_all, const char* RA_RA_active, const char* RA_RA_super, const char* RA_I, const char* I_RA, const char* mature); // read network structure from files
	void read_connections_from_net(const char *filename, std::vector<unsigned int>** target_id, std::vector<double>** target_G); // read connection from .net file
	void read_all_connections_from_net(const char *filename, double*** weights); // read connection from .net file

	void initialize_inhibitory_clusters(); // initialize coordinates of inhibitory clusters on square lattice
	void initialize_RA_for_inh_clusters(); // initialize coordinates of pool neurons for inhibitory clusters
	void initialize_coordinates(); // initialize coordinates of neurons
	void initialize_clusters(); // initialize clusters of neurons
	void initialize_generator(); // initialize generator for processes
	void initialize_equal_clusters(); // initialize equal clusters of neurons

    void initialize_connections(double Gei_mean, double Gei_var, double Gie_mean, double Gie_var); // initialize connections for neurons
	void initialize_connections_for_clusters(double Gei_mean, double Gei_var, double Gie_mean, double Gie_var); // initialize connections for neurons in equal clusters
	void initialize_connections_for_inhibitory_clusters(double Gei_mean, double Gei_var, double Gie_mean, double Gie_var); // initialize connections for inhibitory clusters
	
	void send_connections(); // send fixed connections connections to all processes
	void send_RARA_connections(); // send RA to RA connections to all processes
	
	void update_synaptic_info(); // update synaptic information after reading network structure from files

	void test_connections(); // function to test how connections are written to file
    void test_initialization(double Gie_max, double Gei_max); // test initialization of the pool

    void initialize_test_allI2RA_connections(double Gie); // initialize all equal I to RA connections except for training neurons
    void initialize_test_allRA2I_connections(double Gei); // initialize all equal RA to I connections except for training neurons

	void trial(bool training); // make one trial
	void ground_state(unsigned N_trials); // make simulation of a ground state;

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
	void write_mature(const char* filename); // write mature neurons

    void write_pajek_super(const char* filename); // write supersynapses to a file for pajek
	void write_pajek_active(const char* filename); // write active synapses to a file for pajek
	void write_pajek_all(const char* filename); // write all synapses to a file for pajek
    void write_pajek_fixed(const char* filename); // write fixed synapses to a file for pajek

	// set functions
	void set_generator(Poisson_noise* g); // set poisson noise generator
	void set_generator4neurons(); // set noise generator for all neurons
	void set_dynamics(double interval, double tS); // set dynamics range for all neurons
	void set_no_noise(); // disable Poisson noise for all neurons
	void set_no_noise_RA(); // disable noise for RA neurons
	void set_no_noise_I(); // disable noise for I neurons
    void set_white_noise_distribution_soma(double mu, double sigma); // set white noise distribution for somatic compartment of RA neurons
    void set_white_noise_distribution_dend(double mu, double sigma); // set white noise distribution for dendritic compartment of RA neurons
    void set_white_noise_soma(); // set white noise for somatic compartment of RA neurons
    void set_white_noise_dend(); // set white noise for dendritic compartment of RA neurons
   
   // motivation
    void enable_motivation_noise(); // switch on dc component of the input noise
    void disable_motivation_noise(); // turn off dc component of the input noise
    
    // statistics
    void statistics();
   // current
    void set_training_current(); // set current to training neurons
    void set_testing_current(); // set testing current for neurons
protected:
		// number of neurons
		int N_TR; // number of training HVC(RA) neurons
		int N_RA; // number of HVC(RA) neurons
		int N_I; // number of HVC(I) neurons
		int N_RA_local; // number of RA neurons in each process
		int N_I_local; // number of I neurons in each process
		int Nss; // number of super synapses
		int num_inh_clusters; // number of inhibitory clusters in a raw or column
		int num_inh_in_cluster; // number of inhibitory neurons in single cluster

		HH2_final_pool* HVCRA_local; // array of HVC(RA) neurons
		HHI_final_pool* HVCI_local; // array of HVC(I) neurons

		// coordinates info
		vector <double> xx_RA; // x-coordinates of RA neurons
		vector <double> yy_RA; // y-coordinates of RA neurons
		vector <double> xx_I; // x-coordinates of I neurons
		vector <double> yy_I; // y-coordinates of I neurons
		vector <double> xx_I_cluster; // x-coordinates of centers of inhibitory clusters
		vector <double> yy_I_cluster; // y-coordinates of centers of inhibitory clusters

		const static double g_KICK; // glutamate kick value
		const static double CLUSTER_SIZE; // cluster size
        const static double MIN_INTERNEURON_DISTANCE; // minimum distance between neurons
		const static double LAMBDA_RA2I; // spatial scale of probability of connections decay (exponential)
		const static double A_RA2I; // constant for nearby connections
		const static double B_RA2I; // constant for far away connections
		const static double MEAN_RA2I; // mean of Gaussian distribution for far away connections
		const static double SIGMA_RA2I; // variance of Gaussian distribution for far away connections

		const static double LAMBDA_I2RA; // spatial scale of probability of connections decay
		const static double CONNECT_CONST_I2RA;
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
		double internal_time = 0; // internal time
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
		vector<unsigned>* syn_ID_RA_I_global; // array with synaptic ID numbers from RA to I neurons
		vector<unsigned>* syn_ID_I_RA_global; // array with synaptic ID numbers from I to RA neurons

		vector<unsigned>* active_synapses_global; // array of vectors with IDs of active synapses
		vector<unsigned>* active_supersynapses_global; // array of vectors with IDs of supersynapses

		// update array for conductances and glutamate

		double* update_Ge_RA_local;
        double* update_Gi_RA_local;
		double* update_g_local;
        double* update_Ge_I_local;

        double* update_Ge_RA_global;
        double* update_Gi_RA_global;
		double* update_g_global;
        double* update_Ge_I_global;

		// spatial info
		int* N_cluster; // number of RA neurons per inhibitory center

        // spike time info
        const static int NUM_SOMA_SPIKES; // number of last somatic spikes to store

        std::deque<double>* last_soma_spikes_local; // last NUM_SOMA_SPIKES spikes in somatic compartment
        std::deque<double>* last_soma_spikes_global; // last NUM_SOMA_SPIKES spikes in somatic compartment

        double* spike_times_soma_local; // local array with neurons somatic spike times
		double* spike_times_dend_local; // local array with neurons sdendritic pike times

		double* spike_times_soma_global; // array with the most recent somatic spike times of neurons
		double* spike_times_dend_global; // array with the most recent dendritic spike times of neurons

        vector<double>* spikes_in_trial_soma_global; // array with somatic spike times in the last trial
		vector<double>* spikes_in_trial_dend_global; // array with dendritic spike times in the last trial

        vector<double>* spikes_in_trial_soma_local; // array with somatic spike times in the last trial
		vector<double>* spikes_in_trial_dend_local; // array with dendritic spike times in the last trial

		bool** active_local; // array of HVC(RA) neurons with active synapses
		bool** supersynapses_local; // indicator array for active supersynapses;
		double** weights_local; // array of synaptic strength of all connections between RA neurons
		bool* remodeled_local; // indicators if neuron underwent axon remodelling
		int* mature_local; // indicator if neuron is mature due to supersynaptic acquisition (neuron matures forever after it acquires all supersynapses)
		int* mature_global; // global array of indicators if neuron is mature due to supersynaptic acquisition 
		
		vector<double>* spikes_in_trial_local; // local rray with spike times of single trial
		vector<double>* weights_RA_I_local; // array with synapses from RA to I neurons
		vector<double>* weights_I_RA_local; // array with senapses from I to RA neurons
		vector<unsigned>* syn_ID_RA_I_local; // array with synaptic ID numbers from RA to I neurons
		vector<unsigned>* syn_ID_I_RA_local; // array with synaptic ID numbers from I to RA neurons

		vector<unsigned>* active_synapses_local; // array of vectors with IDs of active synapses
		vector<unsigned>* active_supersynapses_local; // array of vectors with IDs of supersynapses
		unsigned* Id_RA_local; // Id of RA neurons in each process
		unsigned* Id_I_local; // Id of I neurons in each process

		const static double ACTIVATION; // activation threshold for synapses
		double p_RA2I(int i, int j); // probability of connection from RA to I neuron
		double p_I2RA(int i, int j); // probability of connection from I to RA neuron
		double p_RA2I_distant(int i, int j); // probability of distant connections from RA to I

		const static double SUPERSYNAPSE_THRESHOLD; // threshold for supersynaptic connection

		// developmental GABA switch
		const static double T_GABA; // time scale of maturation
		//double E_GABA(double t); // time-dependent switch
		double E_GABA(int n); // activity-dependent switch
		const static double E_GABA_IMMATURE; // immature reverse GABA potential
        const static double E_GABA_MATURE; // mature reverse GABA potential
        const static int N_MATURATION; // maturation scale for number of active synapses

		// constants for STDP-rules
		const static double LTP_WINDOW; // window for LTP beyond which we ignore all somatic spikes
		const static double LTD_WINDOW; // window for LTD beyond which we ignore all somatic spikes


		const static double G_MAX; // constant for maximum weight value
		const static double BETA; // constant for potentiation decay

		const static double A_P;
		const static double G_P;
		const static double T_P;
		const static double TAU_P;

		const static double A_D;
		const static double T_D;
		const static double TAU_D;

		const static double R;
		const static double F_0;

		void LTD(double &w, double t); // long-term depression STDP rule
		void LTP(double &w, double t); // long-term potentiation STDP rule
		void potentiation_decay(); // apply potentiation decay to all RA-RA synapses
		void STDP(int i, int j, double ti, double tj); // STDP rule for connection from neuron i to neuron j
		void STDP_one2one(unsigned i, unsigned j); // apply STDP rules for connection from neuron i to neuron j
		void STDP_saturated(unsigned n); // apply STDP rules for saturated neuron
		void update_synapse(int i, int j, double w); // update synapse from neuron i to neuron j
		void update_all_synapses(); // update synapses between RA neurons and apply potentiation decay
		void axon_remodeling(int i); // remove all targets from neuron i except for supersynapses

        // MPI support

        int MPI_size; // number of processes
        int MPI_rank; // rank of the process

        int N_RA_process; // number of RA neurons in one process
        int N_I_process; // number of I neurons in one process

	void get_neuronRA_location(unsigned n, int* rank, int* shift); // get location of RA neuron with ID n in terms of process and position in array
	void get_neuronI_location(unsigned n, int* rank, int *shift); // get location of I neuron with ID n in terms of process and position in array



};

static double distance(double x1, double y1, double x2, double y2)
{
	double d;
	d = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	return d;
}



#endif // POOLPARALLEL_H_INCLUDED
