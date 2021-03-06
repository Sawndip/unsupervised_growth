#pragma once
#ifndef HH2_BUFFER_H
#define HH2_BUFFER_H
#include <vector>
#include <functional>
#include <cmath>

#define point_distribution_size 30

using std::vector;

typedef std::function<double (double)> DDfunction;

class Poisson_noise;
class HHI_buffer;

class HH2_buffer
{
public:
	HH2_buffer();
	
	int check_bad_values(); // checks if voltage values are in reasonable range
	
	// recording
	void set_recording(std::string filename); // record state of the neuron. Dump recordings to file filename
	void set_recording_full(std::string filename); // record full state of the neuron. Dump recordings to file filename
	
	// print
	void print_param(); // print parameters of the system
	void print_targets(); // print targets of a neuron

	// get internal state
	bool get_fired_soma(); // get neuron's somatic fired state
	bool get_fired_dend(); // get neuron's dendritic fired state
	
	int get_spike_number_soma(); // get number of somatic spikes
	int get_spike_number_dend(); // get number of dendritic spikes
	
	std::vector<double> get_Vs(); // get voltage of somatic compartment
	std::vector<double> get_Vd(); // get voltage of dendritic compartment
    
    // set functions
    void set_original_model_parameters(); // set neuron model parameters as in original model (mature HVC-RA)
    
    void set_c(double c); // set gating variable for Ca dependent K current
    void set_inh_conductance(double G); // set dendritic inhibitory conductance
	void set_Ei(double E){_Egaba = E;}; // set GABA reverse potential
	void set_Erest(double E){EsL = E; EdL = E;}; // set resting potential for somatic and dendritic compartments
	void set_Rc(double R){Rc = R;}; // set coupling resistance between somatic and dendritic compartments
	void set_Ad(double A){Ad = A;}; // set dednritic area
	
	void set_GCa(double G){GdCa = G;}; // set dendritic Ca conductance
	void set_GCaK(double G){GdCaK = G;}; // set dendritic Ca-dependent potassium conductance
	void set_Gk(double G){GsK = G;}; // set somatic potassium conductance
	void set_GNa(double G){GsNa = G;}; // set somatic sodium conductance
	void set_GsL(double G){GsL = G;}; // set somatic leak conductance
	void set_GdL(double G){GdL = G;}; // set dendritic leak conductance
	
	
	
	void reset_time(); // reset time of neuron to zero
	
    void set_soma_current(DDfunction fs); // set current to somatic compartment
	void set_dend_current(DDfunction fd); // set current to dendritic compartment

    // noise
    void set_no_poisson_noise(); // disable Poisson point noisy spikes
	void set_poisson_noise(); // enable Poisson point noisy spikes
	void set_no_white_noise(); // disable white noise
	void set_white_noise(double mu_s, double sigma_s, double mu_d, double sigma_d); // enable white noise

    void set_noise_generator(Poisson_noise* g); // set Poisson noise generator

	// set dynamics
	void set_dynamics(double tS); // set timestep
	void set_to_rest(); // set all variables to the rest state
	void reset(); // reset neuron activity. Last values for all variables are assinged to the first elements 
	
	void Debraband_step(); // one Debraband step; determenistic order 4, weak convergence order 3
	void Debraband_step_no_target_update(); // Debraband step with no update of targets

	void Runge4_step(); // one step of Runge Kutta order 4
	void R4_step_no_target_update(); // Runge Kutta 4 step with no update of targets
	void R4_step_with_target_update(); // Runge Kutta 4 step with update of targets
	
	// change conductance
	void raiseExcWeight(double w) {_Gexc_d += w / Ad;}; // raise excitatory conductance by weight / Ad
	void raiseE(double G); // raise excitatory conductance
	void raiseI(double G); // raise inhibitory conductance

    // targets
	void set_targetI(HHI_buffer* target, int n, double G); // set inhibitory target
	void set_targetRA(HH2_buffer* target, int n, double G); // set excitatory target

protected:
	// buffer
	const static int BUFFER_SIZE; // size of the buffers
	int _buffer_counter; // internal counter used to check buffer overfill
	
	// model parameters

	double cm;	//	membrane capacitance
	double Rc;	// resistance of coupling between soma and dendrite
	double As;	//	soma's membrane area
	double GsL;	//	soma's leak conductance
	double GsNa;	//	soma's sodium channel conductance
	double GsK;	//	soma's potassium channel conductance
	double EsL;	//	soma's leaky reversal potential
	double EsNa;	//	soma's sodium reversal potential
	double EsK;	//	soma's potassium reversal potential
	double Ad;	//	dendrite's membrane area
	double GdL;	//	dendrite's leak conductance
	double GdCa;	//	dendrite's calcium channel conductance
	double GdCaK;	//	dendrite's potassium channel conductance driven by calcium concentration
	double EdL;	//	dendrite's leaky reversal potential
	double EdCa;	//	dendrite's calcium reversal potential
	double EdK;	//	dendrite's potassium reversal potential
	double tExc;	//	time constant for excitatory conductance
	double tInh;	//	time constant for inhibitory conductance

	// dynamic variables
	double _time; // time array
	double _Vs; // soma's membrane potential
	double _Vd; // dendrite's membrane potential
	double _n;	//	gating variable for soma's potassium channel
	double _h;	//	gating variable for soma's sodium channel
	double _r;	//	gating variable for dendrite's calcium channel
	double _c;	//	gating variable for dendrite's potassium channel driven by calcium concentration
	double _Ca;	//	calcium concentration
	double _Is; //	external current to soma
	double _Id;	//	extrenal current to dendrite
	double _Gexc_s; //	excitatory soma conductance
	double _Ginh_s;	//	inhibitory soma conductance
	double _Gexc_d; //	excitatory dendrite conductance
	double _Ginh_d;	//	inhibitory dendrite conductance
	double _Egaba; // reverse GABA potential
	
	std::string _filename; // output filename
	
	// buffers
	std::vector<double> buffer_time;
	std::vector<double> buffer_Vs;
	std::vector<double> buffer_Vd;
	std::vector<double> buffer_Ged;
	std::vector<double> buffer_Gid;
	
	// buffers using for full state recording
	std::vector<double> buffer_n;
	std::vector<double> buffer_h;
	std::vector<double> buffer_r;
	std::vector<double> buffer_c;
	std::vector<double> buffer_Ca;
	
	// constants
	
	//	thresholds
	const static double threshold_spike; // threshold for somatic spike
	const static double threshold_burst; // threshold for dendritic spike
	const static double spike_margin; // New spike is produced only if voltage first goes below threshold - spike_margin
									  // and then crosses the threshold
									  
	const static double burst_margin; // Burst ends only if voltage first goes below threshold_burst - burst_margin
									  // and then crosses the threshold
	
	// functions for external current
	DDfunction _Is_training; // training current to soma
	DDfunction _Id_training; // training current to dendrite
	
	bool _training_dend; // indicator for training current to dendritic compartment
	bool _training_soma; // indicator for training current to somatic compartment

	// internal state
	int _Nspikes;	//	number of spikes occured during dynamics
	bool _fired_soma;

	int _Nbursts; // number of spikes in dendritic compartment
	bool _fired_dend; // state of dendritic compartment

	// functions to check internal state
	void noise_check(double& G, double G_noise, double lambda, double& noise_time); // check if noisy spike came to the neuron's input
	void state_check(); // check if neuron crossed potential threshold
	void state_noise_check(); // check both noise and state

	// dynamics
	double _timestep;	//	time step used to calculate dynamics of neuron
	int _flag_soma; //	flag which stores the state of neuron (fired or not fired yet)
    int _flag_dend; // flag for dendritic compartment
	
	// recording
	bool _recorded; // indicator that neuron is being recorded
	bool _recorded_full; // indicator that all neuron variables are being recorded
	
	// noise
	Poisson_noise *_generator;	//	pointer to Poisson noise generator
	double _noise_exc_soma;	//	time of the next noisy spike in excitatory somatic compartment
	double _noise_inh_soma;	//	time of the next noisy spike in inhibitory somatic compartment
	double _noise_exc_dend;	//	time of the next noisy spike in excitatory dendritic compartment
	double _noise_inh_dend;	//	time of the next noisy spike in inhibitory dendritic compartment

	const static double lambda_inh;	//	parameter of Poisson point process for inhibitory kicks
	const static double lambda_exc; // parameter of Poisson point process for excitatory kicks
	const static double Gs_noise_inh;	//	maximum amplitude of inhibitory noise added to somatic compartment conductances
	const static double Gd_noise_inh;	//	maximum amplitude of inhibitory noise added to dendritic compartment conductances
	const static double Gs_noise_exc;	//	maximum amplitude of excitatory noise added to somatic compartment conductances
	const static double Gd_noise_exc;	//	maximum amplitude of excitatory noise added to dendritic compartment conductances
	bool _poisson_noise;	//	true for noise on, false for noise off

	void initialize_poisson_noise(double& noise_time, double lambda); // initialize noise

	// white noise
	double _mu_soma;
	double _sigma_soma;
	double _mu_dend;
	double _sigma_dend;

	double point_distribution[point_distribution_size]; // array with point distribution for Debraband method

	// targets

	vector<HH2_buffer*> _targets_RA;	//	all neuron's connections to subsequent neurons
	vector<int> _targetsID_RA;	//	vector which contains all ID numbers of target neurons
	vector<double> _targetsG_RA;	//	vector with conductances of all target synapses

	vector<HHI_buffer*> _targets_I;	//	all neuron's connections to inhibitory HVC neurons;
	vector<int> _targetsID_I;	//	vector which contains all ID numbers of target HVCI neurons
	vector<double> _targetsG_I;	//	vector with conductances of all target HVCI synapses

	void postsyn_update();	//	change all postsynaptic conductances


	// write and record
	void write_buffers(); // flush buffer to a file
	void write_buffers_full(); // flush full buffer to a file

	void record_state(); // record neuron state variables
	void record_state_full(); // record neuron full state variables

	//	Additional functions for Runge Kutta's method

	double kVs(double vs, double vd, double n, double h, double t);
	double kVd(double vs, double vd, double r, double c, double ca, double t);
	double kn(double vs,  double n);
	double kh(double vs, double h);
	double kr(double vd, double r);
	double kc(double vd, double c);
	double kCa(double vd, double r, double ca);

	//	Conductance and current functions

	double Ge_s(double t);	//	excitatory soma conductance
	double Ge_d(double t);	//	excitatory dendrite conductance
	double Gi_s(double t);	//	inhibitory soma conductance
	double Gi_d(double t);	//	inhibitory dendrite conductance

	double IsExt(double t); // injected current to somatic compartment
	double IdExt(double t); // injected current to dendritic compartment
	double Is_default(double t){return 0;};	//	default external current for soma
	double Id_default(double t){return 0;};	//	default external current for dendrite

	static double nInf(double v){return 1 / (1 + exp(-(v + 35) / 10));} 
	static double tauN(double v){return 0.1 + 0.5 / (1 + exp((v + 27) / 15));} 
	static double hInf(double v){return 1 / (1 + exp((v + 45) / 7));} 
	static double tauH(double v){return 0.1 + 0.75 / (1 + exp((v + 40.5) / 6));} 
	static double mInf(double v){return 1 / (1 + exp(-(v + 30) / 9.5));} 
	static double rInf(double v){return 1 / (1 + exp(-(v + 5) / 10));} // original +5; spontaneous burst +15
	static double tauR(double v){return 1.0;}  // original 1
	static double cInf(double v){return 1 / (1 + exp(-(v - 10) / 7));} // original -10; spontaneous burst -0
	static double tauC(double v){return 10;} // original 10
};


#endif
