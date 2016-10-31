#include "poisson_noise.h"
#include "HH2_final.h"
#include "HHI_final.h"
#include <iostream>
#include <fstream>

using namespace std::placeholders;

// NMDA extrasynaptic channel kinetics
const double HH2_final::beta = 6.6; // s^-1

// extracellular magnesium concentration
const double HH2_final::C = 1; // mM

// time bins for white noise
const double HH2_final::bin_size = 1.0; // bin size for current white noise stimulus
const double HH2_final::bin_size_nmda = 0.01; // bin size for simulating fraction of NMDA open channels

// noise conductances
const double HH2_final::Gs_noise_inh = 0.0;
const double HH2_final::Gd_noise_inh = 0.0;

const double HH2_final::cm = 1;
const double HH2_final::Rc = 0.055;
const double HH2_final::As = 50;
const double HH2_final::GsL = 0.1;
const double HH2_final::GsNa = 60;
const double HH2_final::GsK = 8;
const double HH2_final::EsL = -80;
const double HH2_final::EsNa = 55;
const double HH2_final::EsK = -90;
const double HH2_final::Ad = 100;
const double HH2_final::GdL = 0.1;
const double HH2_final::GdCa = 55;
const double HH2_final::GdCaK = 150;
const double HH2_final::EdL = -80;
const double HH2_final::EdCa = 120;
const double HH2_final::EdK = -90;
const double HH2_final::tExc = 5;
const double HH2_final::tInh = 7;

const double HH2_final::threshold = 0;
const double HH2_final::threshold_dend = 0;

HH2_final::HH2_final()
{
	// parameters of the model
	Ei = -80;

	//internal state
	itime = 0;
	Nspikes = 0;
	fired_soma = false;

	Nspikes_dend = 0;
	fired_dend = false;

	G_kick = 0.3;
	kick_time = 45000;

	// white_noise

	stored_soma = 0; // stored somatic current value
	stored_dend = 0; // stored dendritic current value

	stored_nmda_soma = 0; // stored somatic nmda open fraction value
	stored_nmda_dend = 0; // stored dendritic nmda open fraction value

	// noise
	noise = true;
	lambda_exc = 75;
	lambda_inh = 100;

	training_dend = false;
	training_soma = false;

	white_noise_soma = false;
	white_noise_dend = false;

	// external current
	//IsExt = std::bind(&HH2_final::Is_default, this, _1);
	//IdExt = std::bind(&HH2_final::Id_default, this, _1);
}


void HH2_final::set_soma_current(DDfunction fs)
{
	training_soma = true;
	Is_training = fs;
}

void HH2_final::set_dend_current(DDfunction fd)
{
	training_dend = true;
	Id_training = fd;
}

void HH2_final::set_no_noise()
{
	noise = false;
}

void HH2_final::initialize_noise(int& noise_time, double lambda)
{
	noise_time = (int) round(1000 * generator->get_spike_time(lambda) / timeStep);
	while (noise_time == 0)
		noise_time = (int) round(1000 * generator->get_spike_time(lambda) / timeStep);
}

void HH2_final::set_kick(double G)
{
    G_kick = G;
}


void HH2_final::set_Ei(double E)
{
    Ei = E;
}

void HH2_final::set_dynamics(double interval, double tS)
{
	timeStep = tS;
	size = (int) round(interval / timeStep) + 1;

	//	Resize all arrays

	time.resize(size);
	Vs.resize(size);
	Is.resize(size);
	n.resize(size);
	h.resize(size);
	Vd.resize(size);
	Id.resize(size);
	Ca.resize(size);
	r.resize(size);
	c.resize(size);
	Gexc_s.resize(size);
	Ginh_s.resize(size);
	G_AMPA.resize(size);
	G_NMDA.resize(size);
	Gexc_d.resize(size);
	Ginh_d.resize(size);
	flag_soma.resize(size);
    flag_dend.resize(size);
    E_gaba.resize(size);

	// set initial values
	time[0] = 0;
	Vs[0] = -79.97619025;
	Vd[0] = -79.97268759;
	n[0] = 0.01101284;
	h[0] = 0.9932845;
	r[0] = 0.00055429;
	c[0] = 0.00000261762353;
	Ca[0] = 0.01689572;
	Ginh_s[0] = 0.0;
	Gexc_s[0] = 0.0;
	Ginh_d[0] = 0.0;
	Gexc_d[0] = 0.0;
	G_AMPA[0] = 0.0;
	G_NMDA[0] = 0.0;
	E_gaba[0] = Ei;

	flag_soma[0] = 0;
 	flag_dend[0] = 0;

	Is[0] = IsExt(time[0]);
	Id[0] = IdExt(time[0]);
	
	// set constants for NMDA noise
	/*
	double tau = 1.0 / (alpha*T + beta); // time constant
	double sInf = alpha*T*tau; // fraction of open channels at time = infinity
	
	mu = exp(-timeStep / tau);

	if (nmda_soma > 0)
		A_soma = sqrt(sInf*(1-sInf)/nmda_soma) * sqrt(1 - mu*mu);
	else
		A_soma = 0.0;
	
	if (nmda_dend > 0)
		A_dend = sqrt(sInf*(1-sInf)/nmda_dend) * sqrt(1 - mu*mu);
	else
		A_dend = 0.0;
	*/

	//	initialize up noise

	this->initialize_noise(noise_exc_soma, lambda_exc);
	this->initialize_noise(noise_inh_soma, lambda_inh);
	this->initialize_noise(noise_exc_dend, lambda_exc);
	this->initialize_noise(noise_inh_dend, lambda_inh);

}

void HH2_final::set_to_rest()
{
	itime = 0;

	time[0] = 0;
	Vs[0] = -79.97619025;
	Vd[0] = -79.97268759;
	n[0] = 0.01101284;
	h[0] = 0.9932845;
	r[0] = 0.00055429;
	c[0] = 0.00000261762353;
	Ca[0] = 0.01689572;
	Ginh_s[0] = 0.0;
	Gexc_s[0] = 0.0;
	Ginh_d[0] = 0.0;
	Gexc_d[0] = 0.0;
	G_AMPA[0] = 0.0;
	G_NMDA[0] = 0.0;
	E_gaba[0] = Ei;

	flag_soma[0] = 0;
    flag_dend[0] = 0;

	Is[0] = IsExt(time[0]);
	Id[0] = IdExt(time[0]);

	//	initialize up noise

	this->initialize_noise(noise_exc_soma, lambda_exc);
	this->initialize_noise(noise_inh_soma, lambda_inh);
	this->initialize_noise(noise_exc_dend, lambda_exc);
	this->initialize_noise(noise_inh_dend, lambda_inh);

}

void HH2_final::set_noise_generator(Poisson_noise* g)
{
	generator = g;
}

void HH2_final::set_white_noise_distribution_soma(double mu, double sigma)
{
	mu_soma = mu;
	sigma_soma = sigma;
}

void HH2_final::set_white_noise_distribution_dend(double mu, double sigma)
{
	mu_dend = mu;
	sigma_dend = sigma;
}

void HH2_final::set_white_noise_soma()
{
	white_noise_soma = true;
}

void HH2_final::set_white_noise_dend()
{
	white_noise_dend = true;
}

double HH2_final::IdExt(double t)
{
	if ((white_noise_dend)&&(training_dend))
	{
		return I_white_noise_dend(t) + Id_training(t);
	}
	else
	{
		if (white_noise_dend)
			return I_white_noise_dend(t);
		else
			return Id_default(t);
	}
}


double HH2_final::IsExt(double t)
{
	if ((white_noise_soma)&&(training_soma))
	{
		return I_white_noise_soma(t) + Is_training(t);
	}
	else
	{
		if (white_noise_soma)
			return I_white_noise_soma(t);
		else
			return Is_default(t);
	}
}

void HH2_final::kick_dend_Ginh()
{
    Ginh_d[itime] += G_kick;
}

void HH2_final::kick_soma_Ginh()
{
    Ginh_s[itime] += G_kick;
}

void HH2_final::set_targetRA(HH2_final *target, int n, double G)
{
	targets_RA.push_back(target);
	targetsID_RA.push_back(n);
	targetsG_RA.push_back(G);
}

void HH2_final::set_silent_targetRA(HH2_final *target, int n, double G)
{
	silent_targets_RA.push_back(target);
	silent_targetsID_RA.push_back(n);
	silent_targetsG_RA.push_back(G);
}

void HH2_final::set_targetI(HHI_final *target, int n, double G)
{
	targets_I.push_back(target);
	targetsID_I.push_back(n);
	targetsG_I.push_back(G);
}

std::vector<double> HH2_final::get_Vs(){return Vs;}
std::vector<double> HH2_final::get_Vd(){return Vd;}

double HH2_final::get_spike_time()
{
	return spike_time;
}

bool HH2_final::get_fired_soma()
{
	return fired_soma;
}

bool HH2_final::get_fired_dend()
{
	return fired_dend;
}

int HH2_final::get_spike_number_soma()
{
    return Nspikes;
}

int HH2_final::get_spike_number_dend()
{
    return Nspikes_dend;

}

void HH2_final::print_targets()
{
	std::cout << "Neuron has "<< targetsID_I.size()
			  << "connections to inhibitory neurons and " << targetsID_RA.size()
			  << "connections to excitatory neurons" << std::endl;
	std::cout << std::endl << "Inhibitory targets: " << std::endl;

	for (unsigned i = 0; i < targetsID_I.size(); i++)
		std::cout << "Inh. target " << targetsID_I[i] << "with synaptic conductance "
				   << targetsG_I[i] << std::endl;

	std::cout << std::endl << "Excitatory targets: " << std::endl;

	for (unsigned i = 0; i < targetsID_RA.size(); i++)
		std::cout << "Exc. target " << targetsID_RA[i] << "with synaptic conductance "
				   << targetsG_RA[i] << std::endl;
}


void HH2_final::print_param()
{
	std::cout << "cm = " << cm << std::endl;
	std::cout << "Rc = " << Rc << std::endl;
	std::cout << "As = " << As << std::endl;
	std::cout << "GsL = " << GsL << std::endl;
	std::cout << "GsNa = " << GsNa << std::endl;
	std::cout << "GsK = " << GsK << std::endl;
	std::cout << "EsL = " << EsL << std::endl;
	std::cout << "EsNa = " << EsNa << std::endl;
	std::cout << "EsK = " << EsK << std::endl;
	std::cout << "Ad = " << Ad << std::endl;
	std::cout << "GdL = " << GdL << std::endl;
	std::cout << "GdCa = " << GdCa << std::endl;
	std::cout << "GdCaK = " << GdCaK << std::endl;
	std::cout << "EdL = " << EdL << std::endl;
	std::cout << "EdCa = " << EdCa << std::endl;
	std::cout << "EdK = " << EdK << std::endl;
}

void HH2_final::writeToFile(const char * filename)
{
	std::ofstream output;

	output.open(filename, std::ios::out | std::ios::binary); //	open file to write binary data

	//	write all parameters of the neuron

	output.write(reinterpret_cast<char*>(&size), sizeof(size));
	output.write(reinterpret_cast<char*>(&Nspikes), sizeof(Nspikes));

	output.write(reinterpret_cast<char*>(&Nspikes_dend), sizeof(Nspikes_dend));
	output.write(reinterpret_cast<char*>(&timeStep), sizeof(timeStep));

	//	Now let's form data array

	for (unsigned i = 0; i < time.size(); ++i)
	{
		output.write(reinterpret_cast<char*>(&time[i]), sizeof(time[i]));
		output.write(reinterpret_cast<char*>(&Vs[i]), sizeof(Vs[i]));
		output.write(reinterpret_cast<char*>(&Is[i]), sizeof(Is[i]));
		output.write(reinterpret_cast<char*>(&n[i]), sizeof(n[i]));
		output.write(reinterpret_cast<char*>(&h[i]), sizeof(h[i]));
		output.write(reinterpret_cast<char*>(&Vd[i]), sizeof(Vd[i]));
		output.write(reinterpret_cast<char*>(&Id[i]), sizeof(Id[i]));
		output.write(reinterpret_cast<char*>(&r[i]), sizeof(r[i]));
		output.write(reinterpret_cast<char*>(&c[i]), sizeof(c[i]));
		output.write(reinterpret_cast<char*>(&Ca[i]), sizeof(Ca[i]));
		output.write(reinterpret_cast<char*>(&Gexc_d[i]), sizeof(Gexc_d[i]));
		output.write(reinterpret_cast<char*>(&Ginh_d[i]), sizeof(Ginh_d[i]));
		output.write(reinterpret_cast<char*>(&Gexc_s[i]), sizeof(Gexc_s[i]));
		output.write(reinterpret_cast<char*>(&Ginh_s[i]), sizeof(Ginh_s[i]));
		output.write(reinterpret_cast<char*>(&G_NMDA[i]), sizeof(G_NMDA[i]));
		output.write(reinterpret_cast<char*>(&E_gaba[i]), sizeof(E_gaba[i]));
		output.write(reinterpret_cast<char*>(&flag_soma[i]), sizeof(flag_soma[i]));
	}

	output.close(); //	close file

}

void HH2_final::raise_AMPA(double G)
{
	G_AMPA[itime] = G_AMPA[itime] + G;
}

void HH2_final::raise_NMDA(double G)
{
	G_NMDA[itime] = G_NMDA[itime] + G * HH2_final::B(Vd[itime]);
}

void HH2_final::raiseI(double G)
{
	Ginh_d[itime] = Ginh_d[itime] + G;
}

void HH2_final::postsyn_update()
{
	for (int i = 0; i < targets_RA.size(); i++)
	{
		targets_RA[i]->raise_AMPA(targetsG_RA[i]);
	}

	for (int i = 0; i < silent_targets_RA.size(); i++)
	{
		silent_targets_RA[i]->raise_NMDA(silent_targetsG_RA[i]);
	}


	for (int i = 0; i < targets_I.size(); i++)
	{
		targets_I[i]->raiseE(targetsG_I[i]);
	}
}


void HH2_final::R4_dend_kick()
{
	if (itime % kick_time == 0)
        this->kick_dend_Ginh();

	this->state_noise_check();

	this->Runge4_step();
}

void HH2_final::R4_soma_kick()
{
	if (itime % kick_time == 0)
        this->kick_soma_Ginh();

	this->state_noise_check();

	this->Runge4_step();
}

double HH2_final::I_white_noise_soma(double t)
{
	if (itime % (int) round(bin_size/timeStep) == 0)
	{
		//printf("(itime+1)*timeStep = %f\tbin size * count = %f\n", (itime + 1)*timeStep, bin_size*count);
		//printf("bin_size = %f\tcount = %d\n", bin_size, count);

		stored_soma = mu_soma + sigma_soma * generator->normal_distribution();
		//printf("stored = %f\n", stored);
		//printf("Sample new\titime = %d\n", itime);
		return stored_soma;
	}
	else
	{
		//printf("stored = %f\n", stored);
		//printf("Sample old\titime = %d\n", itime);
		return stored_soma;

	}
}

double HH2_final::I_white_noise_dend(double t)
{
    //printf("(int) round(bin_size/timeStep) = %d\n",(int) round(bin_size/timeStep));
	if (itime % (int) round(bin_size/timeStep) == 0)
	{
		//printf("(itime+1)*timeStep = %f\tbin size * count = %f\n", (itime + 1)*timeStep, bin_size*count);
		//printf("bin_size = %f\tcount = %d\n", bin_size, count);

		stored_dend = mu_dend + sigma_dend * generator->normal_distribution();
		//printf("stored = %f\n", stored);
		//printf("Sample new\titime = %d\n", itime);
		return stored_dend;
	}
	else
	{
		//printf("stored = %f\n", stored);
		//printf("Sample old\titime = %d\n", itime);
		return stored_dend;

	}
}

void HH2_final::R4_step_with_target_update()
{
	this->state_noise_check();

	if (this->get_fired_soma())
		this->postsyn_update();

	this->Runge4_step();
}

void HH2_final::R4_step_no_target_update()
{
	this->state_noise_check();
	this->Runge4_step();
}

void HH2_final::state_noise_check()
{
	this->state_check();

	if (noise)
	{
		//this->noise_check(gs[itime], MEAN_gs, SIGMA_gs, lambda_exc, noise_exc_soma);
		this->noise_check(Ginh_s[itime], Gs_noise_inh, lambda_inh, noise_inh_soma);
		//this->glutamate_noise_check(gd[itime], MEAN_gd, SIGMA_gd, lambda_exc, noise_exc_dend);
		this->noise_check(Ginh_d[itime], Gd_noise_inh, lambda_inh, noise_inh_dend);
	}
}

void HH2_final::state_check()
{
	if ((flag_soma[itime] == 1) && (Vs[itime] < threshold))
	{
		spike_time = time[itime];
		flag_soma[itime + 1] = 0;
		fired_soma = true;
		Nspikes = Nspikes + 1;
	}
	else
	{
		//	check if we should change the state of neuron (voltage crossed the threshold)
		if ((flag_soma[itime] == 0) && (Vs[itime] > threshold))
			flag_soma[itime + 1] = 1;
		else
		{
			fired_soma = false;
			flag_soma[itime + 1] = flag_soma[itime];	//	otherwise state hasn't changed
		}
	}

	if ((flag_dend[itime] == 1) && (Vd[itime] < threshold_dend))
	{
		spike_time_dend = time[itime];
		flag_dend[itime + 1] = 0;
		fired_dend = true;
		Nspikes_dend = Nspikes_dend + 1;
	}
	else
	{
		//	check if we should change the state of neuron (voltage crossed the threshold)
		if ((flag_dend[itime] == 0) && (Vd[itime] > threshold_dend))
			flag_dend[itime + 1] = 1;
		else
		{
			fired_dend = false;
			flag_dend[itime + 1] = flag_dend[itime];	//	otherwise state hasn't changed
		}
	}
}

void HH2_final::noise_check(double& G, double G_noise, double lambda, int& noise_time)
{
	if (itime == noise_time)
		{
			G += generator->random(G_noise);

			int random = round(1000 * generator->get_spike_time(lambda) / timeStep);
			while (random == 0)
			{
					random = round(1000 * generator->get_spike_time(lambda) / timeStep);
					G += generator->random(G_noise);
			}
			noise_time = noise_time + random;
		}

}


void HH2_final::glutamate_noise_check(double& g, double mean, double sigma, double lambda, int& noise_time)
{

	if (itime == noise_time)
		{
			double random_g = mean + sigma * generator->normal_distribution();
			int random_time = round(1000 * generator->get_spike_time(lambda) / timeStep);
			

			while (random_time == 0)
			{
					random_time = round(1000 * generator->get_spike_time(lambda) / timeStep);
					random_g = mean + sigma * generator->normal_distribution();
			}

			noise_time = noise_time + random_time;

			if (random_g > g)
				g = random_g;
		}


}

void HH2_final::Runge4_step()
{
	double n1, h1, c1, r1;	//	temporary values of gating variables
	double vts, vtd, Cat;	//	temporary values of variables
	double k1Vs, k2Vs, k3Vs, k4Vs;
	double k1Vd, k2Vd, k3Vd, k4Vd;
	double k1Ca, k2Ca, k3Ca, k4Ca;
	double k1n, k2n, k3n, k4n;
	double k1h, k2h, k3h, k4h;
	double k1c, k2c, k3c, k4c;
	double k1r, k2r, k3r, k4r;
	double t;

	vts = Vs[itime];
	n1 = n[itime];
	h1 = h[itime];
	vtd = Vd[itime];
	r1 = r[itime];
	c1 = c[itime];
	Cat = Ca[itime];
	t = time[itime];

	k1Vs = kVs(vts, vtd, n1, h1, t);
	k1n = kn(vts, n1);
	k1h = kh(vts, h1);
	k1Vd = kVd(vts, vtd, r1, c1, Cat, t);
	k1r = kr(vtd, r1);
	k1c = kc(vtd, c1);
	k1Ca = kCa(vtd, r1, Cat);

	vts = Vs[itime] + timeStep * k1Vs / 3;
	n1 = n[itime] + timeStep * k1n / 3;
	h1 = h[itime] + timeStep * k1h / 3;
	vtd = Vd[itime] + timeStep * k1Vd / 3;
	r1 = r[itime] + timeStep * k1r / 3;
	c1 = c[itime] + timeStep * k1c / 3;
	Cat = Ca[itime] + timeStep * k1Ca / 3;

	t = time[itime] + timeStep / 3;

	k2Vs = kVs(vts, vtd, n1, h1, t);
	k2n = kn(vts, n1);
	k2h = kh(vts, h1);
	k2Vd = kVd(vts, vtd, r1, c1, Cat, t);
	k2r = kr(vtd, r1);
	k2c = kc(vtd, c1);
	k2Ca = kCa(vtd, r1, Cat);

	vts = Vs[itime] + timeStep * (-k1Vs / 3 + k2Vs);
	n1 = n[itime] + timeStep * (-k1n / 3 + k2n);
	h1 = h[itime] + timeStep * (-k1h / 3 + k2h);
	vtd = Vd[itime] + timeStep * (-k1Vd / 3 + k2Vd);
	r1 = r[itime] + timeStep * (-k1r / 3 + k2r);
	c1 = c[itime] + timeStep * (-k1c / 3 + k2c);
	Cat = Ca[itime] + timeStep * (-k1Ca / 3 + k2Ca);

	t = time[itime] + 2 * timeStep / 3;

	k3Vs = kVs(vts, vtd, n1, h1, t);
	k3n = kn(vts, n1);
	k3h = kh(vts, h1);
	k3Vd = kVd(vts, vtd, r1, c1, Cat, t);
	k3r = kr(vtd, r1);
	k3c = kc(vtd, c1);
	k3Ca = kCa(vtd, r1, Cat);

	vts = Vs[itime] + timeStep * (k1Vs - k2Vs + k3Vs);
	n1 = n[itime] + timeStep * (k1n - k2n + k3n);
	h1 = h[itime] + timeStep * (k1h - k2h + k3h);
	vtd = Vd[itime] + timeStep * (k1Vd - k2Vd + k3Vd);
	r1 = r[itime] + timeStep * (k1r - k2r + k3r);
	c1 = c[itime] + timeStep * (k1c - k2c + k3c);
	Cat = Ca[itime] + timeStep * (k1Ca - k2Ca + k3Ca);

	t = time[itime] + timeStep;

	k4Vs = kVs(vts, vtd, n1, h1, t);
	k4n = kn(vts, n1);
	k4h = kh(vts, h1);
	k4Vd = kVd(vts, vtd, r1, c1, Cat, t);
	k4r = kr(vtd, r1);
	k4c = kc(vtd, c1);
	k4Ca = kCa(vtd, r1, Cat);

	time[itime + 1] = time[itime] + timeStep;
	
	Vs[itime+1] = Vs[itime] + timeStep * (k1Vs + 3 * k2Vs + 3 * k3Vs + k4Vs) / 8;
	n[itime + 1] = n[itime] + timeStep * (k1n + 3 * k2n + 3 * k3n + k4n) / 8;
	h[itime + 1] = h[itime] + timeStep * (k1h + 3 * k2h + 3 * k3h + k4h) / 8;
	Ginh_s[itime + 1] = Gi_s(time[itime + 1]);
	Vd[itime + 1] = Vd[itime] + timeStep * (k1Vd + 3 * k2Vd + 3 * k3Vd + k4Vd) / 8;
	r[itime + 1] = r[itime] + timeStep * (k1r + 3 * k2r + 3 * k3r + k4r) / 8;
	c[itime + 1] = c[itime] + timeStep * (k1c + 3 * k2c + 3 * k3c + k4c) / 8;
	Ca[itime + 1] = Ca[itime] + timeStep * (k1Ca + 3 * k2Ca + 3 * k3Ca + k4Ca) / 8;
	
	G_NMDA[itime + 1] = G_NMDA_future(Vd[itime + 1], time[itime + 1]);
	G_AMPA[itime + 1] = G_ampa(time[itime + 1]);
	
	Gexc_d[itime + 1] = G_AMPA[itime + 1] + G_NMDA[itime + 1];
	Ginh_d[itime + 1] = Gi_d(time[itime + 1]);
	
	Id[itime + 1] = IdExt(time[itime + 1]);
	Is[itime + 1] = IsExt(time[itime + 1]);
	E_gaba[itime + 1] = Ei;

	itime = itime + 1;
}


double HH2_final::G_ampa(double t){return G_AMPA[itime] * exp(-(t - time[itime]) / tExc);}
double HH2_final::Gi_s(double t){return Ginh_s[itime] * exp(-(t - time[itime]) / tInh);}
double HH2_final::Gi_d(double t){return Ginh_d[itime] * exp(-(t - time[itime]) / tInh);}


//double HH2_final::B(double v){return 1.0 / (1 + exp(-0.062 * v) * C / 3.57);}

double HH2_final::G_NMDA_future(double vd, double t)
{
	return G_NMDA[itime] * HH2_final::B(vd) * exp(- beta * (t - time[itime]) / 1000.0) / HH2_final::B(Vd[itime]);
}

/*double HH2_final::Gs_NMDA(double vs)
{
	
	if (itime % (int) round(bin_size_nmda/timeStep) == 0)
	{

		stored_nmda_soma = G_channel * (nmda_soma*p0 + generator->normal_distribution()*sqrt(nmda_soma*p0*(1-p0))) / ((1 + exp(-0.062 * vs) * C / 3.57) * (1000 * As));
		return stored_nmda_soma;
	}
	else
	{
		return stored_nmda_soma;
	}
}

double HH2_final::Gd_NMDA(double vd)
{
	
	if (itime % (int) round(bin_size_nmda/timeStep) == 0)
	{

		stored_nmda_dend = G_channel * (nmda_dend*p0 + generator->normal_distribution()*sqrt(nmda_dend*p0*(1-p0))) / ((1 + exp(-0.062 * vd) * C / 3.57) * (1000*Ad));
		return stored_nmda_dend;
	}
	else
	{
		return stored_nmda_dend;
	}
}*/
double HH2_final::kVs(double vs, double vd, double n, double h, double t)
{
	double m3, n4;

	m3 = mInf(vs) * mInf(vs) * mInf(vs);
	n4 = n * n * n * n;
	return (-GsL * (vs - EsL) - GsNa * m3 * h * (vs - EsNa) - GsK * n4 * (vs - EsK)
		 - Gi_s(t) * (vs - Ei) + IsExt(t) / As + (vd - vs) / (Rc * As)) / cm;
}
double HH2_final::kVd(double vs, double vd, double r, double c, double ca, double t)
{
	double r2;

	r2 = r * r;
	return (-GdL * (vd - EdL) - GdCa * r2 * (vd - EdCa) - GdCaK * (vd - EdK) * c / (1 + 6 / ca)
		- (G_NMDA_future(vd, t) + G_ampa(t)) * vd - Gi_d(t) * (vd - Ei) + IdExt(t) / Ad + (vs - vd) / (Rc * Ad)) / cm;
}

double HH2_final::kn(double vs, double n){return (HH2_final::nInf(vs) - n) / HH2_final::tauN(vs);}
double HH2_final::kh(double vs, double h){return (HH2_final::hInf(vs) - h) / HH2_final::tauH(vs);}
double HH2_final::kr(double vd, double r){return (HH2_final::rInf(vd) - r) / HH2_final::tauR(vd);}
double HH2_final::kc(double vd, double c){return (HH2_final::cInf(vd) - c) / HH2_final::tauC(vd);}
double HH2_final::kCa(double vd, double r, double ca)
	{return (-0.1 * GdCa * r * r * (vd - EdCa) - 0.02 * ca);}
