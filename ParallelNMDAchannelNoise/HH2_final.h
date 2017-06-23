#pragma once
#ifndef HH2_FINAL_H
#define HH2_FINAL_H
#include <vector>
#include <functional>
#include <cmath>

#define point_distribution_size 30

using std::vector;

typedef std::function<double (double)> DDfunction;

class Poisson_noise;
class HHI_final;

class HH2_final
{
public:
	HH2_final();
	HH2_final(DDfunction fs, DDfunction fd);

	// print
	void print_param(); // print parameters of the system
	void print_targets(); // print targets of a neuron

	void writeToFile(const char * filename); // write data to a binary file

	// get internal state
	bool get_fired_soma(); // get neuron's somatic fired state
	bool get_fired_dend(); // get neuron's dendritic fired state
	double get_spike_time(); // get time of the last spike
	int get_spike_number_soma(); // get number of somatic spikes
	int get_spike_number_dend(); // get number of dendritic spikes
	std::vector<double> get_Vs(); // get voltage of somatic compartment
	std::vector<double> get_Vd(); // get voltage of dendritic compartment
    
    // set functions
	void set_Ei(double E); // set GABA reverse potential
	void set_time(double t); // set initial time
	
    void set_soma_current(DDfunction fs); // set current to somatic compartment
	void set_dend_current(DDfunction fd); // set current to dendritic compartment

    // noise
    void set_no_poisson_noise(); // disable Poisson point noisy spikes
	void set_poisson_noise(); // enable Poisson point noisy spikes
	void set_no_white_noise(); // disable white noise
	void set_white_noise(double mu_s, double sigma_s, double mu_d, double sigma_d); // enable white noise

    void set_noise_generator(Poisson_noise* g); // set Poisson noise generator

	// set dynamics
	void set_dynamics(double interval, double tS); // set sizes of all arrays
	void set_to_rest(); // set all variables to the rest state
	void renew_neuron(); // renew neuron by setting all variables to steady state values and setting spike numbers to zero
	void reset(); // reset neuron activity. Last values for all variables are assinged to the first elements 
	
	// calculate dynamics
	void Euler_step(); // one Euler step
	void Euler_step_no_target_update(); // Euler step with no update of targets
	
	void EulerMaryama_step(); // one Euler-Maryama step
	void EulerMaryama_step_no_target_update(); // Euler-Maryama step with no update of targets
	
	void Debraband_step(); // one Debraband step; determenistic order 4, weak convergence order 3
	void Debraband_step_no_target_update(); // Debraband step with no update of targets

	void DRI1_step(); // one DRI1 step; determenistic order 3, weak convergence order 3
	void DRI1_step_no_target_update(); // DRI1 step with no update of targets
	
	void Runge4_step(); // one step of Runge Kutta order 4
	void R4_step_no_target_update(); // Runge Kutta 4 step with no update of targets
	void R4_step_with_target_update(); // Runge Kutta 4 step with update of targets

	void Runge6_step(); // one step of Runge Kutta order 6
	void R6_step_no_target_update(); // Runge Kutta 6 step with no update of targets
	
	// change conductance
	void raiseE(double G); // raise excitatory conductance
	void raiseI(double G); // raise inhibitory conductance

    // targets
	void set_targetI(HHI_final* target, int n, double G); // set inhibitory target
	void set_targetRA(HH2_final* target, int n, double G); // set excitatory target

protected:
	// model parameters

	const static double cm;	//	membrane capacitance
	const static double Rc;	// resistance of coupling between soma and dendrite
	const static double As;	//	soma's membrane area
	const static double GsL;	//	soma's leak conductance
	const static double GsNa;	//	soma's sodium channel conductance
	const static double GsK;	//	soma's potassium channel conductance
	const static double EsL;	//	soma's leaky reversal potential
	const static double EsNa;	//	soma's sodium reversal potential
	const static double EsK;	//	soma's potassium reversal potential
	const static double Ad;	//	dendrite's membrane area
	const static double GdL;	//	dendrite's leak conductance
	const static double GdCa;	//	dendrite's calcium channel conductance
	const static double GdCaK;	//	dendrite's potassium channel conductance driven by calcium concentration
	const static double EdL;	//	dendrite's leaky reversal potential
	const static double EdCa;	//	dendrite's calcium reversal potential
	const static double EdK;	//	dendrite's potassium reversal potential
	double Ei;	//	inhibitory reverse potential
	const static double tExc;	//	time constant for excitatory conductance
	const static double tInh;	//	time constant for inhibitory conductance

	// dynamic variables
	vector<double> time; // time array
	vector<double> Vs; // soma's membrane potential
	vector<double> Vd; // dendrite's membrane potential
	vector<double> n;	//	gating variable for soma's potassium channel
	vector<double> h;	//	gating variable for soma's sodium channel
	vector<double> r;	//	gating variable for dendrite's calcium channel
	vector<double> c;	//	gating variable for dendrite's potassium channel driven by calcium concentration
	vector<double> Ca;	//	calcium concentration
	vector<double> Is; //	external current to soma
	vector<double> Id;	//	extrenal current to dendrite
	vector<double> Gexc_s; //	excitatory soma conductance
	vector<double> Ginh_s;	//	inhibitory soma conductance
	vector<double> Gexc_d; //	excitatory dendrite conductance
	vector<double> Ginh_d;	//	inhibitory dendrite conductance
	vector<double> E_gaba; // reverse GABA potential
	// constants
	
	//	thresholds
	const static double threshold; // threshold for somatic spike
	const static double threshold_dend; // threshold for dendritic spike
	const static double spike_margin; // spike margin used to make sure that voltage fluctuations don't lead to
									  // the emergence of non-present dendritic and somatic burst. 
									  // New spike can now be registered only if membrane potential first does below 
									  // spike_threshold - spike_margin and then crosses the spike threshold
	
	
	// functions for external current
	DDfunction Is_training; // training current to soma
	DDfunction Id_training; // training current to dendrite
	
	bool training_dend; // indicator for training current to dendritic compartment
	bool training_soma; // indicator for training current to somatic compartment

	// internal state
	int itime;	//	internal time
	int Nspikes;	//	number of spikes occured during dynamics
	double spike_time; // time of the most recent AP spike
	bool fired_soma;

	int Nspikes_dend; // number of spikes in dendritic compartment
	double spike_time_dend; // time of the last spike in dendritic compartment
    bool fired_dend; // state of dendritic compartment

	// functions to check internal state
	void noise_check(double& G, double G_noise, double lambda, int& noise_time); // check if noisy spike came to the neuron's input
	void state_check(); // check if neuron crossed potential threshold
	void state_noise_check(); // check both noise and state

	// dynamics
	double timeStep;	//	time step used to calculate dynamics of neuron
	int size; // size of dynamic variables' vectors
	vector<int> flag_soma; //	flag which stores the state of neuron (fired or not fired yet)
    vector<int> flag_dend; // flag for dendritic compartment
	// noise

	Poisson_noise * generator;	//	pointer to Poisson noise generator
	int noise_exc_soma;	//	time of the next noisy spike in excitatory somatic compartment
	int noise_inh_soma;	//	time of the next noisy spike in inhibitory somatic compartment
	int noise_exc_dend;	//	time of the next noisy spike in excitatory dendritic compartment
	int noise_inh_dend;	//	time of the next noisy spike in inhibitory dendritic compartment

	double lambda_inh;	//	parameter of Poisson point process for inhibitory kicks
	double lambda_exc; // parameter of Poisson point process for excitatory kicks
	const static double Gs_noise_inh;	//	maximum amplitude of inhibitory noise added to somatic compartment conductances
	const static double Gd_noise_inh;	//	maximum amplitude of inhibitory noise added to dendritic compartment conductances
	const static double Gs_noise_exc;	//	maximum amplitude of excitatory noise added to somatic compartment conductances
	const static double Gd_noise_exc;	//	maximum amplitude of excitatory noise added to dendritic compartment conductances
	bool poisson_noise;	//	true for noise on, false for noise off

	void initialize_poisson_noise(int& noise_time, double lambda); // initialize noise

	// white noise
	double mu_soma;
	double mu_dend;
	double sigma_soma;
	double sigma_dend;

	double point_distribution[point_distribution_size]; // array with point distribution for Debraband method

	// targets

	vector<HH2_final*> targets_RA;	//	all neuron's connections to subsequent neurons
	vector<int> targetsID_RA;	//	vector which contains all ID numbers of target neurons
	vector<double> targetsG_RA;	//	vector with conductances of all target synapses

	vector<HHI_final*> targets_I;	//	all neuron's connections to inhibitory HVC neurons;
	vector<int> targetsID_I;	//	vector which contains all ID numbers of target HVCI neurons
	vector<double> targetsG_I;	//	vector with conductances of all target HVCI synapses

	void postsyn_update();	//	change all postsynaptic conductances

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
	static double rInf(double v){return 1 / (1 + exp(-(v + 15) / 10));}
	static double tauR(double v){return 1;}
	static double cInf(double v){return 1 / (1 + exp(-(v - 0) / 7));}
	static double tauC(double v){return 10;}
};


#endif
