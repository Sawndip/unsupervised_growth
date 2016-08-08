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
	bool get_fired_soma(); // get neuron's somatic fired state
	bool get_fired_dend(); // get neuron's dendritic fired state
	double get_spike_time(); // get time of the last spike
	int get_spike_number_soma(); // get number of somatic spikes
	int get_spike_number_dend(); // get number of dendritic spikes

	// noise
	void set_noise(double Gs_exc, double Gs_inh, double Gd_exc, double Gd_inh); // set noise for the neuron
	void set_Ei(double E); // set GABA reverse potential
	void set_glutamate(double g); // set dendritic glutamate
	void set_no_noise(); // disable Poisson point noisy spikes
	void set_noise_generator(Poisson_noise* g); // set Poisson noise generator
	void set_white_noise_soma(); // set white noise to somatic compartment
	void set_white_noise_dend(); // set white noise to denndritic compartment
	void set_white_noise_distribution_soma(double mu, double sigma); // set parameters of white noise distribution for somatic compartment
    void set_white_noise_distribution_dend(double mu, double sigma); // set parameters of white noise distribution for dendritic compartment

    void set_kick(double G); // set the amplitude of the inhibitory kick conductance

	// set dynamics
	void set_dynamics(double interval, double tS); // set sizes of all arrays
	void set_to_rest(); // set all variables to the rest state
	
	// calculate dynamics
	void Runge4_step(); // one step of Runge Kutta order 4
	void R4_step_no_target_update(); // Runge Kutta step with no update of targets
	void R4_step_with_target_update(); // Runge Kutta step with update of targets
    void R4_dend_kick(); // Runge Kutta step with kicks in inhibitory dendritic conductance every kick_time
    void R4_soma_kick(); // Runge Kutta step with kicks in inhibitory somatic conductance every kick_time

	// change conductance
	void raiseE(double G); // raise excitatory conductance
	void raiseI(double G); // raise inhibitory conductance
	void raiseGlutamate(double g); // raise dendritic glutamate
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
	const static double tNMDA; // time constant for NMDA excitatory conductance

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
	vector<double> G_NMDA; // NMDA conductance of dendritic compartment
	vector<double> G_AMPA; // AMPA conductance of dendritic compartment
	vector<double> E_gaba; // reverse GABA potential
	vector<double> s_soma; // fluctuations in fraction of open somatic NMDA channels
	vector<double> s_dend; // fluctuations in fraction of open dendritic NMDA channels
	// constants
	
	const static double alpha; // transition rate to open state 
	const static double beta; // transition rate to closed state
	const static double T; // glutamate extracellular concentration
	const static double G_channel; // conductance of open NMDA channel
	//const static double p0; // probability for a channel to be open
	const static double s0; // fraction of open NMDA channels
	const static double C; // magnesium extracellular concentration
	const static double bin_size; // bin size for current white noise stimulus
	const static double bin_size_nmda; // bin size for nmda fraction of open channels simulation
	const static int nmda_soma; // number of NMDA extrasynaptic receptors in somatic compartment
	const static int nmda_dend; // number of NMDA extrasynaptic receptors in dendritic compartment

	//	thresholds
	const static double threshold; // threshold for somatic spike
	const static double threshold_dend; // threshold for dendritic spike
	
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
	void glutamate_noise_check(double& g, double mean, double sigma, double lambda, int& noise_time); // check if glutamate noisy spikes arrived
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
	bool noise;	//	true for noise on, false for noise off
    double G_kick; // amplitude of inhibitory kick conductance

	double G_SPIKE; // kick to NMDA receptors in dendritic compartment after spike of presynaptic neuron


	void initialize_noise(int& noise_time, double lambda); // initialize noise

	// white noise
	double stored_dend;
	double stored_soma;
	double stored_nmda_soma;
	double stored_nmda_dend;
	double mu_soma;
	double mu_dend;
	double sigma_soma;
	double sigma_dend;
	int count;
	bool white_noise_soma; // indicator for white noise to somatic compartment
	bool white_noise_dend; // indicator for white noise to dendritic compartment

	// NMDA channel noise
	double mu;
	double A_soma;
	double A_dend;

	// targets

	vector<HH2_final*> targets_RA;	//	all neuron's connections to subsequent neurons
	vector<int> targetsID_RA;	//	vector which contains all ID numbers of target neurons
	vector<double> targetsG_RA;	//	vector with conductances of all target synapses

	vector<HHI_final*> targets_I;	//	all neuron's connections to inhibitory HVC neurons;
	vector<int> targetsID_I;	//	vector which contains all ID numbers of target HVCI neurons
	vector<double> targetsG_I;	//	vector with conductances of all target HVCI synapses

	void postsyn_update();	//	change all postsynaptic conductances

    // kick
    void kick_dend_Ginh();  // add G_kick to inhibitory dendritic conductance
    void kick_soma_Ginh();  // add G_kick to inhibitory somatic conductance

    int kick_time; // time to kick
	//	Additional functions for Runge Kutta's method

	double kVs(double vs, double vd, double n, double h, double t);
	double kVd(double vs, double vd, double r, double c, double ca, double t);
	double kn(double vs,  double n);
	double kh(double vs, double h);
	double kr(double vd, double r);
	double kc(double vd, double c);
	double kCa(double vd, double r, double ca);

	//	Conductance and current functions

	double Gs_NMDA(double vs, int timeInd);
	double Gd_NMDA(double vd, int timeInd);

	double Gi_s(double t);	//	inhibitory soma conductance
	double Gi_d(double t);	//	inhibitory dendrite conductance
	double G_ampa(double t); // excitatory AMPA conductance of dendritic compartment

	double IsExt(double t); // injected current to somatic compartment
	double IdExt(double t); // injected current to dendritic compartment
	double Is_default(double t){return 0;};	//	default external current for soma
	double Id_default(double t){return 0;};	//	default external current for dendrite
	double I_white_noise_soma(double t); // external white noise current to somatic compartment
    double I_white_noise_dend(double t); // external white noise current to dendritic compartment

	static double B(double v){return 1.0 / (1 + exp(-0.062 * v) * C / 3.57);}
	static double nInf(double v){return 1 / (1 + exp(-(v + 35) / 10));}
	static double tauN(double v){return 0.1 + 0.5 / (1 + exp((v + 27) / 15));}
	static double hInf(double v){return 1 / (1 + exp((v + 45) / 7));}
	static double tauH(double v){return 0.1 + 0.75 / (1 + exp((v + 40.5) / 6));}
	static double mInf(double v){return 1 / (1 + exp(-(v + 30) / 9.5));}
	static double rInf(double v){return 1 / (1 + exp(-(v + 5) / 10));}
	static double tauR(double v){return 1;}
	static double cInf(double v){return 1 / (1 + exp(-(v - 10) / 7));}
	static double tauC(double v){return 10;}
};


#endif
