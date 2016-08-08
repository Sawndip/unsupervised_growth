#pragma once
#ifndef POOL
#define POOL

#include <vector>
#include <cmath>

class Poisson_noise;
class HH2_final_pool;
class HHI_final_pool;

using std::vector;

//typedef std::function<double (double)> DDfunction;

class Pool
{
public:
	Pool(int N_tr, int N_ra, int N_i, int Ns);
	~Pool();
    
    void initialize_coordinates(); // initialize coordinates of neurons in pool
	void test_connections(); // function to test how connections are written to file
    void test_initialization(double Gie_max, double Gei_max); // test initialization of the pool
    void test_initialization_no_connections(); // initialize pool with no connections
    void test_initialization_no_connections2I(double Gie_max); // initialize pool with no connections to inhibitory neurons. Connections to excitatory neurons are all the same.
	void trial(); // make one trial
	void ground_state(unsigned N_trials); // make simulation of a ground state;

	void reset_after_trial(); // reset neurons after trial is over
	void initialize_pool(double Gie_max, double Gei_max); // initialize pools with synapses stored outside neurons
	void visualize_RA_connections(); // visualize RA connections of the network
	void visualize_I_connections(); // visualize I connections of the network

	// write to file functions
	void write_variable_synapses(const char* RA_RA); // write RA to RA active connections
	void write_invariable_synapses(const char* RA_I, const char* I_RA); // write RA to I and I to RA connections
	void write_coordinates(const char* xy_RA, const char* xy_I); // write coordinates of neurons
	void write_RA(const char * filename, int n); // write dynamics of RA neuron to a file
	void write_I(const char * filename, int n); // write dynamics of I neuron to a file
	void write_weights(const char * filename); // write weights of all synapses in network
    void write_time_info(const char* filename); // write simulation time to a file
    void write_state(const char* filename); // write state of the network to file

	// set functions
	void set_generator(Poisson_noise* g); // set poisson noise generator
	void set_generator4neurons(); // set noise generator for all neurons
	void set_dynamics(double interval, double tS); // set dynamics range for all neurons
	void set_no_noise(); // disable Poisson noise for all neurons
    void set_no_noise_RA(); // disable noise for RA neurons
    void set_no_noise_I(); // disable noise for I neurons
    void set_training_current(); // set current to training neurons
protected:
		int N_TR; // number of training HVC(RA) neurons
		int N_RA; // number of HVC(RA) neurons
		int N_I; // number of HVC(I) neurons
		int Nss; // number of super synapses

		HH2_final_pool* HVCRA; // array of HVC(RA) neurons
		HHI_final_pool* HVCI; // array of HVC(I) neurons

		vector <double> xx_RA; // x-coordinates of RA neurons
		vector <double> yy_RA; // y-coordinates of RA neurons
		vector <double> xx_I; // x-coordinates of I neurons
		vector <double> yy_I; // y-coordinates of I neurons

		int size; // time array size for dynamics
		double timeStep; // time step for dynamics
		double internal_time = 0; // internal time
		int trial_number = 0; // number of simulated trials

		Poisson_noise* generator; // poisson noise generator

		bool* remodeled; // array indicating axon remodeling state of neurons
		bool** active; // array of HVC(RA) neurons with active synapses
		bool** active_supersynapses; // indicator array for active supersynapses;
		double** weights; // array of synaptic strength of all connections between RA neurons
		double* spike_times; // array with the most recent spike times of neurons
		vector<double>* spikes_in_trial; // array with spike times of single trial
		vector<double>* weights_RA_I; // array with synapses from RA to I neurons
		vector<double>* weights_I_RA; // array with senapses from I to RA neurons
		vector<unsigned>* syn_ID_RA_I; // array with synaptic ID numbers from RA to I neurons
		vector<unsigned>* syn_ID_I_RA; // array with synaptic ID numbers from I to RA neurons

		vector<unsigned>* active_synapses; // array of vectors with IDs of active synapses
		vector<unsigned>* supersynapses; // array of vectors with IDs of supersynapses

        const static double MIN_INTERNEURON_DISTANCE; // minimum distance between neurons
		const static double LAMBDA; // spatial scale of probability of connections decay
		const static double CONNECT_CONST;
		const static double ACTIVATION; // activation threshold for synapses
		const static double SIDE; // length of HVC side
		double p(int i, int j); // probability of connection to another RA neuron
		const static double SUPERSYNAPSE_THRESHOLD; // threshold for supersynaptic connection

		// developmental GABA switch
		const static double T_GABA; // time scale of maturation
		double E_GABA(double t); // time-dependent switch
		double E_GABA(int n); // activity-dependent switch
		const static double E_GABA_IMMATURE; // immature reverse GABA potential
        const static double E_GABA_MATURE; // mature reverse GABA potential
        const static int N_MATURATION; // maturation scale for number of active synapses

		// constants for STDP-rules
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
		void STDP(unsigned n, unsigned index_in_fired_array, std::vector<unsigned>::iterator neurons_fired_iter); // apply spike-dependent plasticity rules to for single neuron fired
		void STDP_one2one(unsigned i, unsigned j); // apply STDP rules for connection from neuron i to neuron j
		void STDP_saturated(unsigned n); // apply STDP rules for saturated neuron
		void update_synapses(); // update synapses between RA neurons and apply potentiation decay
		void axon_remodeling(); // remove all targets except for supersynapses if neuron is saturated
};

static double distance(double x1, double y1, double x2, double y2)
{
	double d;
	d = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	return d;
}

#endif
