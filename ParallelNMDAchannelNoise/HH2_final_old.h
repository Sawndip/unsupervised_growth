#pragma once
#ifndef HH2_FINAL_H
#define HH2_FINAL_H
#include <vector>
#include <functional>
#include <cmath>

using std::vector;

typedef std::function<double (double)> DDfunction;

class Poisson_noise;
class HHI_final;

class HH2_final
{
public:
	HH2_final();
	HH2_final(bool white_noise_soma, bool white_noise_dendrite);
	HH2_final(DDfunction fs, DDfunction fd);

	// print
	void print_param(); // print parameters of the system
	void print_targets(); // print targets of a neuron

	void writeToFile(const char * filename); // write data to a binary file

	// set external injected current
	void set_soma_current(DDfunction fs); // set current to somatic compartment
	void set_dend_current(DDfunction fd); // set current to dendritic compartment

	// get internal state
	bool get_fired_soma(); // get neuron's fired state of somatic compartment
	double get_spike_time_soma(); // get time of the last spike in somatic compartment
	bool get_fired_dend(); // get neuron's fired state of dendritic compartment
	double get_spike_time_dend(); // get time of the last spike in dendritic compartment

	// noise
	void set_no_noise(); // disable Poisson point noisy spikes
	void set_noise_generator(Poisson_noise* g); // set Poisson noise generator
	void set_white_noise_soma(); // set white noise to somatic compartment
	void set_white_noise_dend(); // set white noise to denndritic compartment
	void set_white_noise_distribution(double mu, double sigma); // set parameters of white noise distribution

	// set dynamics
	void set_dynamics(double interval, double tS); // set sizes of all arrays

	// calculate dynamics
	void Runge4_step(); // one step of Runge Kutta order 4
	void R4_step_no_target_update(); // Runge Kutta step with no update of targets
	void R4_step_with_target_update(); // Runge Kuttta step with update of targets

	// change conductance
	void raiseE(double G); // raise excitatory conductance
	void raiseI(double G); // raise inhibitory conductance

	// targets
	void set_targetI(HHI_final* target, int n, double G); // set inhibitory target
	void set_targetRA(HH2_final* target, int n, double G); // set excitatory target


protected:
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
	double Ei;	//	inhibitory reverse potential
	double tExc;	//	time constant for excitatory conductance
	double tInh;	//	time constant for inhibitory conductance

	double threshold_soma;	//	threshold for action potential in somatic compartment
    double threshold_dendrite; // threshold for action potential in dendritic compartment

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
    vector<double> E_gaba; // reverse gaba potential

	// functions for external current
	DDfunction IsExt; // current to soma
	DDfunction IdExt; // current to dendrite

	// internal state
	int itime;	//	internal time
	int Nspikes;	//	number of spikes occured during dynamics
	double spike_time_soma; // time of the most recent AP spike in somatic compartment
	bool fired_soma;
    double spike_time_dend; // time of the most recent AP spike in dendritic compartment
	bool fired_dendrite;

	// functions to check internal state
	void noise_check(double& G, double G_noise, int& noise_time); // check if noisy spike came to the neuron's input
	void state_check(); // check if neuron crossed potential threshold
	void state_noise_check(); // check both noise and state

	// dynamics
	double timeStep;	//	time step used to calculate dynamics of neuron
	int size; // size of dynamic variables' vectors
	vector<int> flag_soma; //	flag which stores the state of neuron's somatic compartment (fired or not fired yet)
    vector<int> flag_dend; //	flag which stores the state of neuron's dendritic compartment (fired or not fired yet)

	// noise

	Poisson_noise * generator;	//	pointer to Poisson noise generator
	int noise_exc_soma;	//	time of the next noisy spike in excitatory somatic compartment
	int noise_inh_soma;	//	time of the next noisy spike in inhibitory somatic compartment
	int noise_exc_dend;	//	time of the next noisy spike in excitatory dendritic compartment
	int noise_inh_dend;	//	time of the next noisy spike in inhibitory dendritic compartment

	double lambda;	//	parameter of Poisson point process
	double Gs_noise_exc;	//	maximum amplitude of excitatory noise added to somatic compartment conductances
	double Gs_noise_inh;	//	maximum amplitude of inhibitory noise added to somatic compartment conductances
	double Gd_noise_exc;	//	maximum amplitude of excitatory noise added to dendritic compartment conductances
	double Gd_noise_inh;	//	maximum amplitude of inhibitory noise added to dendritic compartment conductances
	bool noise;	//	true for noise on, false for noise off

	void initialize_noise(int& noise_time); // initialize noise

	// white noise
	bool sampled_this_iteration;
	double stored;
	double bin_size;
	int count;

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
	double Gi_s(double t);	//	inhibitory soma conductance
	double Ge_d(double t);	//	excitatory dendrite conductance
	double Gi_d(double t);	//	inhibitory dendrite conductance

	double Is_default(double t){return 0;};	//	default external current for soma
	double Id_default(double t){return 0;};	//	default external current for dendrite
	double I_white_noise(double t); // external white noise current

};

static double nInf(double v){return 1 / (1 + exp(-(v + 35) / 10));}
static double tauN(double v){return 0.1 + 0.5 / (1 + exp((v + 27) / 15));}
static double hInf(double v){return 1 / (1 + exp((v + 45) / 7));}
static double tauH(double v){return 0.1 + 0.75 / (1 + exp((v + 40.5) / 6));}
static double mInf(double v){return 1 / (1 + exp(-(v + 30) / 9.5));}
static double rInf(double v){return 1 / (1 + exp(-(v + 5) / 10));}
static double tauR(double v){return 1;}
static double cInf(double v){return 1 / (1 + exp(-(v - 10) / 7));}
static double tauC(double v){return 10;}

#endif
