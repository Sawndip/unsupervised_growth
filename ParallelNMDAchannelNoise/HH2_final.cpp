#include "poisson_noise.h"
#include "HH2_final.h"
#include "HHI_final.h"
#include "exception.h"
#include <iostream>
#include <algorithm>
#include <fstream>

using namespace std::placeholders;

// noise conductances
const double HH2_final::Gs_noise_inh = 0.010;
const double HH2_final::Gd_noise_inh = 0.010;
const double HH2_final::Gs_noise_exc = 0.010;
const double HH2_final::Gd_noise_exc = 0.010;

const double HH2_final::cm = 1;
const double HH2_final::Rc = 0.055; // original = 0.055
const double HH2_final::As = 50;
const double HH2_final::GsL = 0.1;
const double HH2_final::GsNa = 60;
const double HH2_final::GsK = 8;
const double HH2_final::EsL = -80;
const double HH2_final::EsNa = 55;
const double HH2_final::EsK = -90;
const double HH2_final::Ad = 250;
const double HH2_final::GdL = 0.1;
const double HH2_final::GdCa = 55;
const double HH2_final::GdCaK = 150;
const double HH2_final::EdL = -80;
const double HH2_final::EdCa = 120;
const double HH2_final::EdK = -90;
const double HH2_final::tExc = 5;
const double HH2_final::tInh = 5;

const double HH2_final::threshold = 0;
const double HH2_final::threshold_dend = 0;

HH2_final::HH2_final() : mu_soma(0.0), sigma_soma(0.0), mu_dend(0.0), sigma_dend(0.0), point_distribution{sqrt(6), -sqrt(6), 1, 1, 1, 1, 1, 1, 1, 1, 1, 
																										  -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0,
																										  0, 0, 0, 0, 0, 0, 0}
{
	// parameters of the model
	Ei = -80;

	//internal state
	itime = 0;
	Nspikes = 0;
	fired_soma = false;

	Nspikes_dend = 0;
	fired_dend = false;

	// noise
	poisson_noise = false;

	lambda_exc = 75;
	lambda_inh = 100;

	training_dend = false;
	training_soma = false;

    generator = nullptr;
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

void HH2_final::set_white_noise(double mu_s, double sigma_s, double mu_d, double sigma_d)
{
    try
    {
        if (generator == nullptr)
            throw NoGenerator("Noise generator is not set up!\n");
        else
        {
            mu_soma = mu_s;
            sigma_soma = sigma_s;
            mu_dend = mu_d;
            sigma_dend = sigma_d;
        }
    }

    catch (NoGenerator const& e)
    {
        std::cerr << "NoGenerator Exception: " << e.what() << std::endl;

    }

}

void HH2_final::set_no_white_noise()
{
    mu_soma = 0.0;
	sigma_soma = 0.0;
    mu_dend = 0.0;
	sigma_dend = 0.0;
}

void HH2_final::set_poisson_noise()
{
    poisson_noise = true;
    
    this->initialize_poisson_noise(noise_exc_soma, lambda_exc);
	this->initialize_poisson_noise(noise_inh_soma, lambda_inh);
	this->initialize_poisson_noise(noise_exc_dend, lambda_exc);
	this->initialize_poisson_noise(noise_inh_dend, lambda_inh);
 
}

void HH2_final::set_no_poisson_noise()
{
	poisson_noise = false;
}

void HH2_final::initialize_poisson_noise(int& noise_time, double lambda)
{
    try
    {
        if (generator == nullptr)
            throw NoGenerator("Noise generator is not set up!\n");
        else
        {
	        noise_time = (int) round(1000 * generator->get_spike_time(lambda) / timeStep);
	        while (noise_time == 0)
		        noise_time = (int) round(1000 * generator->get_spike_time(lambda) / timeStep);
        }
    }

    catch (NoGenerator const& e)
    {
        std::cerr << "NoGenerator Exception: " << e.what() << std::endl;

    }

}

void HH2_final::set_time(double t)
{
	time[0] = t;
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
	E_gaba[0] = Ei;

	flag_soma[0] = 0;
 	flag_dend[0] = 0;

	Is[0] = IsExt(time[0]);
	Id[0] = IdExt(time[0]);
	
	//	initialize up noise
    if (poisson_noise)
    {
	    this->initialize_poisson_noise(noise_exc_soma, lambda_exc);
	    this->initialize_poisson_noise(noise_inh_soma, lambda_inh);
	    this->initialize_poisson_noise(noise_exc_dend, lambda_exc);
	    this->initialize_poisson_noise(noise_inh_dend, lambda_inh);
    }

}

void HH2_final::reset()
{
	itime = 0;

	time[0] = time.back();
	Vs[0] = Vs.back();
	Vd[0] = Vd.back();
	n[0] = n.back();
	h[0] = h.back();
	r[0] = r.back();
	c[0] = c.back();
	Ca[0] = Ca.back();
	Ginh_s[0] = Ginh_s.back();
	Gexc_s[0] = Gexc_s.back();
	Ginh_d[0] = Ginh_d.back();
	Gexc_d[0] = Gexc_d.back();
	E_gaba[0] = E_gaba.back();

	flag_soma[0] = 0;
    flag_dend[0] = 0;

	if (poisson_noise)
	{
		noise_exc_soma -= size - 1;
		noise_inh_soma -= size - 1;
		noise_exc_dend -= size - 1;
		noise_inh_dend -= size - 1;
			
		if (noise_exc_soma < 0)
		{
			std::cerr << "Noise time is negative! noise_exc_soma = " << noise_exc_soma << std::endl;
			this->initialize_poisson_noise(noise_exc_soma, lambda_exc);
		}
			
		if (noise_inh_soma < 0)
		{
			std::cerr << "Noise time is negative! noise_inh_soma = " << noise_inh_soma << std::endl;
			this->initialize_poisson_noise(noise_inh_soma, lambda_inh);
		}
		
		if (noise_exc_dend < 0)
		{
			std::cerr << "Noise time is negative! noise_exc_dend = " << noise_exc_dend << std::endl;
			this->initialize_poisson_noise(noise_exc_dend, lambda_exc);
		}
		
		if (noise_inh_dend < 0)
		{
			std::cerr << "Noise time is negative! noise_inh_dend = " << noise_inh_dend << std::endl;
			this->initialize_poisson_noise(noise_inh_dend, lambda_inh);
		}
		
	}
}

void HH2_final::renew_neuron()
{
	itime = 0;

	time[0] = time.back();
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
	E_gaba[0] = Ei;

	flag_soma[0] = 0;
    flag_dend[0] = 0;
    
    Nspikes = 0;
    Nspikes_dend = 0;

	//	initialize up noise
    if (poisson_noise)
    {
	    this->initialize_poisson_noise(noise_exc_soma, lambda_exc);
	    this->initialize_poisson_noise(noise_inh_soma, lambda_inh);
	    this->initialize_poisson_noise(noise_exc_dend, lambda_exc);
	    this->initialize_poisson_noise(noise_inh_dend, lambda_inh);
    }

	
}

void HH2_final::set_to_rest()
{
	itime = 0;

	time[0] = time.back();
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
	E_gaba[0] = Ei;

	flag_soma[0] = 0;
    flag_dend[0] = 0;

	//	initialize up noise
    if (poisson_noise)
    {
	    this->initialize_poisson_noise(noise_exc_soma, lambda_exc);
	    this->initialize_poisson_noise(noise_inh_soma, lambda_inh);
	    this->initialize_poisson_noise(noise_exc_dend, lambda_exc);
	    this->initialize_poisson_noise(noise_inh_dend, lambda_inh);
    }

}

void HH2_final::set_noise_generator(Poisson_noise* g)
{
	generator = g;
}

double HH2_final::IdExt(double t)
{
	if (training_dend)
		return Id_training(t);
	else
		return Id_default(t);
}

double HH2_final::IsExt(double t)
{
	if (training_soma)
		return Is_training(t);
	else
		return Is_default(t);
}

void HH2_final::set_targetRA(HH2_final *target, int n, double G)
{
	targets_RA.push_back(target);
	targetsID_RA.push_back(n);
	targetsG_RA.push_back(G);
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
		output.write(reinterpret_cast<char*>(&E_gaba[i]), sizeof(E_gaba[i]));
		output.write(reinterpret_cast<char*>(&flag_soma[i]), sizeof(flag_soma[i]));
	}

	output.close(); //	close file

}

void HH2_final::raiseE(double G)
{
    Gexc_d[itime] += G;
}

void HH2_final::raiseI(double G)
{
	Ginh_d[itime] = Ginh_d[itime] + G;
}

void HH2_final::postsyn_update()
{
	for (int i = 0; i < targets_RA.size(); i++)
	{
		targets_RA[i]->raiseE(targetsG_RA[i]);
	}

	for (int i = 0; i < targets_I.size(); i++)
	{
		targets_I[i]->raiseE(targetsG_I[i]);
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

void HH2_final::Debraband_step_no_target_update()
{
	this->state_noise_check();
	this->Debraband_step();
}

void HH2_final::DRI1_step_no_target_update()
{
	this->state_noise_check();
	this->DRI1_step();
}

void HH2_final::R6_step_no_target_update()
{
	this->state_noise_check();
	this->Runge6_step();
}

void HH2_final::Euler_step_no_target_update()
{
	this->state_noise_check();
	this->Euler_step();
}

void HH2_final::EulerMaryama_step_no_target_update()
{
	this->state_noise_check();
	this->EulerMaryama_step();
}

void HH2_final::state_noise_check()
{
	this->state_check();

	if (poisson_noise)
	{
		this->noise_check(Ginh_s[itime], Gs_noise_inh, lambda_inh, noise_inh_soma);
		this->noise_check(Ginh_d[itime], Gd_noise_inh, lambda_inh, noise_inh_dend);
		this->noise_check(Gexc_s[itime], Gs_noise_exc, lambda_exc, noise_exc_soma);
		this->noise_check(Gexc_d[itime], Gd_noise_exc, lambda_exc, noise_exc_dend);
	}
}

void HH2_final::state_check()
{
    // somatic spike is defined as membrane potential second crossing of th threshold (when potential gows down)
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

    // dendritic spike is defined as membrane potential first crossing of the threshold (when potential gows up)
	if ((flag_dend[itime] == 0) && (Vd[itime] >= threshold_dend))
	{
		spike_time_dend = time[itime];
		flag_dend[itime + 1] = 1;
		fired_dend = true;
		Nspikes_dend = Nspikes_dend + 1;
	}
	else
	{
		//	check if we should change the state of neuron (voltage crossed the threshold)
		if ((flag_dend[itime] == 1) && (Vd[itime] < threshold_dend))
			flag_dend[itime + 1] = 0;
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

void HH2_final::Euler_step()
{
	double m3;
	double n4;

	m3 = mInf(Vs[itime]) * mInf(Vs[itime]) * mInf(Vs[itime]);
	n4 = n[itime] * n[itime] * n[itime] * n[itime];

	Vs[itime + 1] = Vs[itime] + timeStep * ( (-GsL * (Vs[itime] - EsL) - GsNa * m3 * h[itime] * (Vs[itime] - EsNa) - GsK * n4 * (Vs[itime] - EsK)
		 - Ginh_s[itime] * (Vs[itime] - Ei) - Gexc_s[itime] * Vs[itime] + IsExt(time[itime]) / As + (Vd[itime] - Vs[itime]) / (Rc * As)) / cm );

	Vd[itime + 1] = Vd[itime] + timeStep * ( (-GdL * (Vd[itime] - EdL) - GdCa * r[itime] * r[itime] * (Vd[itime] - EdCa) - 
					GdCaK * (Vd[itime] - EdK) * c[itime] / (1 + 6 / Ca[itime]) - Gexc_d[itime] * Vd[itime] - 
					Ginh_d[itime] * (Vd[itime] - Ei) + IdExt(time[itime]) / Ad + (Vs[itime] - Vd[itime]) / (Rc * Ad)) / cm );

	n[itime + 1] = nInf(Vs[itime]);
	h[itime + 1] = hInf(Vs[itime]);
	r[itime + 1] = rInf(Vd[itime]);
	c[itime + 1] = cInf(Vd[itime]);
	
	Ca[itime + 1] = Ca[itime] + timeStep * (-0.1 * GdCa * r[itime] * r[itime] * (Vd[itime] - EdCa) - 0.02 * Ca[itime]);

	time[itime + 1] = time[itime] + timeStep;

	Gexc_d[itime + 1] = Ge_d(time[itime + 1]);
    Gexc_s[itime + 1] = Ge_s(time[itime + 1]);

	Ginh_d[itime + 1] = Gi_d(time[itime + 1]);
	Ginh_s[itime + 1] = Gi_s(time[itime + 1]);
	
	Id[itime + 1] = IdExt(time[itime + 1]);
	Is[itime + 1] = IsExt(time[itime + 1]);
	E_gaba[itime + 1] = Ei;

	itime = itime + 1;
}

void HH2_final::EulerMaryama_step()
{
	double m3;
	double n4;

	m3 = mInf(Vs[itime]) * mInf(Vs[itime]) * mInf(Vs[itime]);
	n4 = n[itime] * n[itime] * n[itime] * n[itime];

	double Is_random = generator->normal_distribution();
	double Id_random = generator->normal_distribution();

	Id[itime] += sigma_dend * Id_random + mu_dend;
	Is[itime] += sigma_soma * Is_random + mu_soma;
	
	Vs[itime + 1] = Vs[itime] + timeStep * ( (-GsL * (Vs[itime] - EsL) - GsNa * m3 * h[itime] * (Vs[itime] - EsNa) - GsK * n4 * (Vs[itime] - EsK)
		 		  - Ginh_s[itime] * (Vs[itime] - Ei) - Gexc_s[itime] * Vs[itime] + (IsExt(time[itime]) + mu_soma) / As + (Vd[itime] - Vs[itime]) / (Rc * As)) / cm )
		 		  + sqrt(timeStep) * sigma_soma * Is_random / (As * cm);

	Vd[itime + 1] = Vd[itime] + timeStep * ( (-GdL * (Vd[itime] - EdL) - GdCa * r[itime] * r[itime] * (Vd[itime] - EdCa) - 
					GdCaK * (Vd[itime] - EdK) * c[itime] / (1 + 6 / Ca[itime]) - Gexc_d[itime] * Vd[itime] - 
					Ginh_d[itime] * (Vd[itime] - Ei) + (IdExt(time[itime]) + mu_dend) / Ad + (Vs[itime] - Vd[itime]) / (Rc * Ad)) / cm )
					+ sqrt(timeStep) * sigma_dend * Id_random / (Ad * cm);

	n[itime + 1] = nInf(Vs[itime]);
	h[itime + 1] = hInf(Vs[itime]);
	r[itime + 1] = rInf(Vd[itime]);
	c[itime + 1] = cInf(Vd[itime]);
	
	Ca[itime + 1] = Ca[itime] + timeStep * (-0.1 * GdCa * r[itime] * r[itime] * (Vd[itime] - EdCa) - 0.02 * Ca[itime]);

	time[itime + 1] = time[itime] + timeStep;

	Gexc_d[itime + 1] = Ge_d(time[itime + 1]);
    Gexc_s[itime + 1] = Ge_s(time[itime + 1]);

	Ginh_d[itime + 1] = Gi_d(time[itime + 1]);
	Ginh_s[itime + 1] = Gi_s(time[itime + 1]);
	
	E_gaba[itime + 1] = Ei;

	Id[itime + 1] = IdExt(time[itime + 1]);
	Is[itime + 1] = IsExt(time[itime + 1]);
	
	itime = itime + 1;

}

void HH2_final::DRI1_step()
{

	double H1Vs, H2Vs, H3Vs;	
	double H1Vd, H2Vd, H3Vd;	
	double H1Ca, H2Ca, H3Ca;	
	double H1n, H2n, H3n;	
	double H1h, H2h, H3h;	
	double H1c, H2c, H3c;	
	double H1r, H2r, H3r;	

	double a1Vs, a2Vs, a3Vs;	
	double a1Vd, a2Vd, a3Vd;	
	double a1Ca, a2Ca, a3Ca;	
	double a1n, a2n, a3n;	
	double a1h, a2h, a3h;	
	double a1c, a2c, a3c;	
	double a1r, a2r, a3r;	
	
	double Js = generator->normal_distribution();
	double Jd = generator->normal_distribution();

	Id[itime] += mu_dend + Jd * sigma_dend;
	Is[itime] += mu_soma + Js * sigma_soma;
	
	H1Vs = Vs[itime];
	H1Vd = Vd[itime];
	H1Ca = Ca[itime];
	H1n = n[itime];
	H1h = h[itime];
	H1c = c[itime];
	H1r = r[itime];

	a1Vs = kVs(H1Vs, H1Vd, H1n, H1h, time[itime]);	
	a1Vd = kVd(H1Vs, H1Vd, H1r, H1c, H1Ca, time[itime]);
	a1n = kn(H1Vs, H1n);
	a1h = kh(H1Vs, H1h);
	a1r = kr(H1Vd, H1r);
	a1c = kc(H1Vd, H1c);
	a1Ca = kCa(H1Vd, H1r, H1Ca);

	H2Vs = Vs[itime] + timeStep * 0.5 * a1Vs + sqrt(timeStep) * sigma_soma * (6 - sqrt(6) ) * Js / (10 * As * cm);
	H2Vd = Vd[itime] + timeStep * 0.5 * a1Vd + sqrt(timeStep) * sigma_dend * (6 - sqrt(6) ) * Jd / (10 * Ad * cm);
	H2Ca = Ca[itime] + timeStep * 0.5 * a1Ca;
	H2n = n[itime] + timeStep * 0.5 * a1n;
	H2h = h[itime] + timeStep * 0.5 * a1h;
	H2c = c[itime] + timeStep * 0.5 * a1c;
	H2r = r[itime] + timeStep * 0.5 * a1r;

	a2Vs = kVs(H2Vs, H2Vd, H2n, H2h, time[itime] + 0.5 * timeStep);	
	a2Vd = kVd(H2Vs, H2Vd, H2r, H2c, H2Ca, time[itime] + 0.5 * timeStep);
	a2n = kn(H2Vs, H2n);
	a2h = kh(H2Vs, H2h);
	a2r = kr(H2Vd, H2r);
	a2c = kc(H2Vd, H2c);
	a2Ca = kCa(H2Vd, H2r, H2Ca);

	H3Vs = Vs[itime] + timeStep * (- a1Vs + 2 * a2Vs) + sqrt(timeStep) * sigma_soma * (3 + 2 * sqrt(6)) * Js / (5 * As * cm);
	H3Vd = Vd[itime] + timeStep * (- a1Vd + 2 * a2Vd) + sqrt(timeStep) * sigma_dend * (3 + 2 * sqrt(6)) * Jd / (5 * Ad * cm);
	H3Ca = Ca[itime] + timeStep * (- a1Ca + 2 * a2Ca);
	H3n = n[itime] + timeStep * (- a1n + 2 * a2n);
	H3h = h[itime] + timeStep * (- a1h + 2 * a2h);
	H3c = c[itime] + timeStep * (- a1c + 2 * a2c);
	H3r = r[itime] + timeStep * (- a1r + 2 * a2r);

	a3Vs = kVs(H3Vs, H3Vd, H3n, H3h, time[itime] + timeStep);	
	a3Vd = kVd(H3Vs, H3Vd, H3r, H3c, H3Ca, time[itime] + timeStep);
	a3n = kn(H3Vs, H3n);
	a3h = kh(H3Vs, H3h);
	a3r = kr(H3Vd, H3r);
	a3c = kc(H3Vd, H3c);
	a3Ca = kCa(H3Vd, H3r, H3Ca);

	time[itime + 1] = time[itime] + timeStep;
	
	Vs[itime + 1] = Vs[itime] + timeStep * (a1Vs + 4 * a2Vs + a3Vs) / 6
				  + sqrt(timeStep) * sigma_soma * Js / (As * cm);
	
	Vd[itime + 1] = Vd[itime] + timeStep * (a1Vd + 4 * a2Vd + a3Vd) / 6
				  + sqrt(timeStep) * sigma_dend * Jd / (Ad * cm);
	
	Ca[itime + 1] = Ca[itime] + timeStep * (a1Ca + 4 * a2Ca + a3Ca) / 6;
	n[itime + 1] = n[itime] + timeStep * (a1n + 4 * a2n + a3n) / 6;
	h[itime + 1] = h[itime] + timeStep * (a1h + 4 * a2h + a3h) / 6;
	r[itime + 1] = r[itime] + timeStep * (a1r + 4 * a2r + a3r) / 6;
	c[itime + 1] = c[itime] + timeStep * (a1c + 4 * a2c + a3c) / 6;

	Gexc_d[itime + 1] = Ge_d(time[itime + 1]);
    Gexc_s[itime + 1] = Ge_s(time[itime + 1]);

	Ginh_d[itime + 1] = Gi_d(time[itime + 1]);
	Ginh_s[itime + 1] = Gi_s(time[itime + 1]);
	
	Id[itime + 1] = IdExt(time[itime + 1]);
	Is[itime + 1] = IsExt(time[itime + 1]);
	E_gaba[itime + 1] = Ei;

	itime = itime + 1;


}

void HH2_final::Debraband_step()
{
	double H1Vs, H2Vs, H3Vs, H4Vs;	
	double H1Vd, H2Vd, H3Vd, H4Vd;	
	double H1Ca, H2Ca, H3Ca, H4Ca;	
	double H1n, H2n, H3n, H4n;	
	double H1h, H2h, H3h, H4h;	
	double H1c, H2c, H3c, H4c;	
	double H1r, H2r, H3r, H4r;	

	double a1Vs, a2Vs, a3Vs, a4Vs;	
	double a1Vd, a2Vd, a3Vd, a4Vd;	
	double a1Ca, a2Ca, a3Ca, a4Ca;	
	double a1n, a2n, a3n, a4n;	
	double a1h, a2h, a3h, a4h;	
	double a1c, a2c, a3c, a4c;	
	double a1r, a2r, a3r, a4r;	
	
	//double J1s = generator->normal_distribution();
	//double J2s = generator->normal_distribution();
	//double J1d = generator->normal_distribution();
	//double J2d = generator->normal_distribution();

	int ind1s, ind2s, ind1d, ind2d;

	ind1s = generator->sample_index_for_point_distribution();
	ind2s = generator->sample_index_for_point_distribution();
	ind1d = generator->sample_index_for_point_distribution();
	ind2d = generator->sample_index_for_point_distribution();

	double J1s = point_distribution[ind1s];
	double J2s = point_distribution[ind2s];
	double J1d = point_distribution[ind1d];
	double J2d = point_distribution[ind2d];
	
	Id[itime] += mu_dend + J1d * sigma_dend;
	Is[itime] += mu_soma + J1s * sigma_soma;
	
	H1Vs = Vs[itime] + sqrt(timeStep) * sigma_soma * (-0.01844540496323970 * J1s - 0.1866426386543421 * J2s) / (As * cm);
	H1Vd = Vd[itime] + sqrt(timeStep) * sigma_dend * (-0.01844540496323970 * J1d - 0.1866426386543421 * J2d) / (Ad * cm);
	H1Ca = Ca[itime];
	H1n = n[itime];
	H1h = h[itime];
	H1c = c[itime];
	H1r = r[itime];

	a1Vs = kVs(H1Vs, H1Vd, H1n, H1h, time[itime]);	
	a1Vd = kVd(H1Vs, H1Vd, H1r, H1c, H1Ca, time[itime]);
	a1n = kn(H1Vs, H1n);
	a1h = kh(H1Vs, H1h);
	a1r = kr(H1Vd, H1r);
	a1c = kc(H1Vd, H1c);
	a1Ca = kCa(H1Vd, H1r, H1Ca);

	H2Vs = Vs[itime] + timeStep * a1Vs + sqrt(timeStep) * sigma_soma * (0.8017012756521233 * J1s - 0.8575745885712401 * J2s) / (As * cm);
	H2Vd = Vd[itime] + timeStep * a1Vd + sqrt(timeStep) * sigma_dend * (0.8017012756521233 * J1d - 0.8575745885712401 * J2d) / (Ad * cm);
	H2Ca = Ca[itime] + timeStep * a1Ca;
	H2n = n[itime] + timeStep * a1n;
	H2h = h[itime] + timeStep * a1h;
	H2c = c[itime] + timeStep * a1c;
	H2r = r[itime] + timeStep * a1r;

	a2Vs = kVs(H2Vs, H2Vd, H2n, H2h, time[itime] + timeStep);	
	a2Vd = kVd(H2Vs, H2Vd, H2r, H2c, H2Ca, time[itime] + timeStep);
	a2n = kn(H2Vs, H2n);
	a2h = kh(H2Vs, H2h);
	a2r = kr(H2Vd, H2r);
	a2c = kc(H2Vd, H2c);
	a2Ca = kCa(H2Vd, H2r, H2Ca);

	H3Vs = Vs[itime] + timeStep * (3 * a1Vs + a2Vs) / 8 + sqrt(timeStep) * sigma_soma * (0.5092227024816198 * J1s - 0.4723392695015512 * J2s) / (As * cm);
	H3Vd = Vd[itime] + timeStep * (3 * a1Vd + a2Vd) / 8 + sqrt(timeStep) * sigma_dend * (0.5092227024816198 * J1d - 0.4723392695015512 * J2d) / (Ad * cm);
	H3Ca = Ca[itime] + timeStep * (3 * a1Ca + a2Ca) / 8;
	H3n = n[itime] + timeStep * (3 * a1n + a2n) / 8;
	H3h = h[itime] + timeStep * (3 * a1h + a2h) / 8;
	H3c = c[itime] + timeStep * (3 * a1c + a2c) / 8;
	H3r = r[itime] + timeStep * (3 * a1r + a2r) / 8;

	a3Vs = kVs(H3Vs, H3Vd, H3n, H3h, time[itime] + 0.5 * timeStep);	
	a3Vd = kVd(H3Vs, H3Vd, H3r, H3c, H3Ca, time[itime] + 0.5 * timeStep);
	a3n = kn(H3Vs, H3n);
	a3h = kh(H3Vs, H3h);
	a3r = kr(H3Vd, H3r);
	a3c = kc(H3Vd, H3c);
	a3Ca = kCa(H3Vd, H3r, H3Ca);

	H4Vs = Vs[itime] + timeStep * (-0.4526683126055039 * a1Vs - 0.4842227708685013 * a2Vs + 1.9368910834740051 * a3Vs)
		 + sqrt(timeStep) * sigma_soma * (0.9758794209767762 * J1s + 0.3060354860326548 * J2s) / (As * cm);
	
	H4Vd = Vd[itime] + timeStep * (-0.4526683126055039 * a1Vd - 0.4842227708685013 * a2Vd + 1.9368910834740051 * a3Vd)
		 + sqrt(timeStep) * sigma_dend * (0.9758794209767762 * J1d + 0.3060354860326548 * J2d) / (Ad * cm);

	H4Ca = Ca[itime] + timeStep * (-0.4526683126055039 * a1Ca - 0.4842227708685013 * a2Ca + 1.9368910834740051 * a3Ca);
	H4n = n[itime] + timeStep * (-0.4526683126055039 * a1n - 0.4842227708685013 * a2n + 1.9368910834740051 * a3n);
	H4h = h[itime] + timeStep * (-0.4526683126055039 * a1h - 0.4842227708685013 * a2h + 1.9368910834740051 * a3h);
	H4c = c[itime] + timeStep * (-0.4526683126055039 * a1c - 0.4842227708685013 * a2c + 1.9368910834740051 * a3c);
	H4r = r[itime] + timeStep * (-0.4526683126055039 * a1r - 0.4842227708685013 * a2r + 1.9368910834740051 * a3r);

	a4Vs = kVs(H4Vs, H4Vd, H4n, H4h, time[itime] + timeStep);	
	a4Vd = kVd(H4Vs, H4Vd, H4r, H4c, H4Ca, time[itime] + timeStep);
	a4n = kn(H4Vs, H4n);
	a4h = kh(H4Vs, H4h);
	a4r = kr(H4Vd, H4r);
	a4c = kc(H4Vd, H4c);
	a4Ca = kCa(H4Vd, H4r, H4Ca);

	time[itime + 1] = time[itime] + timeStep;
	
	Vs[itime + 1] = Vs[itime] + timeStep * (a1Vs / 6 - 0.005430430675258792 * a2Vs + 2 * a3Vs / 3 + 0.1720970973419255 * a4Vs)
				  + sqrt(timeStep) * sigma_soma * J1s / (As * cm);
	
	Vd[itime + 1] = Vd[itime] + timeStep * (a1Vd / 6 - 0.005430430675258792 * a2Vd + 2 * a3Vd / 3 + 0.1720970973419255 * a4Vd)
				  + sqrt(timeStep) * sigma_dend * J1d / (Ad * cm);
	
	Ca[itime + 1] = Ca[itime] + timeStep * (a1Ca / 6 - 0.005430430675258792 * a2Ca + 2 * a3Ca / 3 + 0.1720970973419255 * a4Ca);
	n[itime + 1] = n[itime] + timeStep * (a1n / 6 - 0.005430430675258792 * a2n + 2 * a3n / 3 + 0.1720970973419255 * a4n);
	h[itime + 1] = h[itime] + timeStep * (a1h / 6 - 0.005430430675258792 * a2h + 2 * a3h / 3 + 0.1720970973419255 * a4h);
	r[itime + 1] = r[itime] + timeStep * (a1r / 6 - 0.005430430675258792 * a2r + 2 * a3r / 3 + 0.1720970973419255 * a4r);
	c[itime + 1] = c[itime] + timeStep * (a1c / 6 - 0.005430430675258792 * a2c + 2 * a3c / 3 + 0.1720970973419255 * a4c);

	Gexc_d[itime + 1] = Ge_d(time[itime + 1]);
    Gexc_s[itime + 1] = Ge_s(time[itime + 1]);

	Ginh_d[itime + 1] = Gi_d(time[itime + 1]);
	Ginh_s[itime + 1] = Gi_s(time[itime + 1]);
	
	Id[itime + 1] = IdExt(time[itime + 1]);
	Is[itime + 1] = IsExt(time[itime + 1]);
	E_gaba[itime + 1] = Ei;

	itime = itime + 1;
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
	Vd[itime + 1] = Vd[itime] + timeStep * (k1Vd + 3 * k2Vd + 3 * k3Vd + k4Vd) / 8;
	r[itime + 1] = r[itime] + timeStep * (k1r + 3 * k2r + 3 * k3r + k4r) / 8;
	c[itime + 1] = c[itime] + timeStep * (k1c + 3 * k2c + 3 * k3c + k4c) / 8;
	Ca[itime + 1] = Ca[itime] + timeStep * (k1Ca + 3 * k2Ca + 3 * k3Ca + k4Ca) / 8;
	
	Gexc_d[itime + 1] = Ge_d(time[itime + 1]);
    Gexc_s[itime + 1] = Ge_s(time[itime + 1]);

	Ginh_d[itime + 1] = Gi_d(time[itime + 1]);
	Ginh_s[itime + 1] = Gi_s(time[itime + 1]);
	
	Id[itime + 1] = IdExt(time[itime + 1]);
	Is[itime + 1] = IsExt(time[itime + 1]);
	E_gaba[itime + 1] = Ei;

	itime = itime + 1;
}

void HH2_final::Runge6_step()
{
	double n1, h1, c1, r1;	//	temporary values of gating variables
	double vts, vtd, Cat;	//	temporary values of variables
	double k1Vs, k2Vs, k3Vs, k4Vs, k5Vs, k6Vs, k7Vs;
	double k1Vd, k2Vd, k3Vd, k4Vd, k5Vd, k6Vd, k7Vd;
	double k1Ca, k2Ca, k3Ca, k4Ca, k5Ca, k6Ca, k7Ca;
	double k1n, k2n, k3n, k4n, k5n, k6n, k7n;
	double k1h, k2h, k3h, k4h, k5h, k6h, k7h;
	double k1c, k2c, k3c, k4c, k5c, k6c, k7c;
	double k1r, k2r, k3r, k4r, k5r, k6r, k7r;
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
	
	vts = Vs[itime] + 2 * timeStep * k2Vs / 2;
	n1 = n[itime] + 2 * timeStep * k2n / 3;
	h1 = h[itime] + 2 * timeStep * k2h / 3;
	vtd = Vd[itime] + 2 * timeStep * k2Vd / 2;
	r1 = r[itime] + 2 * timeStep * k2r / 3;
	c1 = c[itime] + 2 * timeStep * k2c / 3;
	Cat = Ca[itime] + 2 * timeStep * k2Ca / 3;
	t = time[itime] + 2 * timeStep / 3;

	k3Vs = kVs(vts, vtd, n1, h1, t);
	k3n = kn(vts, n1);
	k3h = kh(vts, h1);
	k3Vd = kVd(vts, vtd, r1, c1, Cat, t);
	k3r = kr(vtd, r1);
	k3c = kc(vtd, c1);
	k3Ca = kCa(vtd, r1, Cat);

	vts = Vs[itime] + timeStep * (k1Vs / 12 + k2Vs / 3 - k3Vs / 12);
	n1 = n[itime] + timeStep * (k1n / 12 + k2n / 3 - k3n / 12);
	h1 = h[itime] + timeStep * (k1h / 12 + k2h / 3 - k3h / 12);
	vtd = Vd[itime] + timeStep * (k1Vd / 12 + k2Vd / 3 - k3Vd / 12);
	r1 = r[itime] + timeStep * (k1r / 12 + k2r / 3 - k3r / 12);
	c1 = c[itime] + timeStep * (k1c / 12 + k2c / 3 - k3c / 12);
	Cat = Ca[itime] + timeStep * (k1Ca / 12 + k2Ca / 3 - k3Ca / 12);
	t = time[itime] + timeStep / 3;

	k4Vs = kVs(vts, vtd, n1, h1, t);
	k4n = kn(vts, n1);
	k4h = kh(vts, h1);
	k4Vd = kVd(vts, vtd, r1, c1, Cat, t);
	k4r = kr(vtd, r1);
	k4c = kc(vtd, c1);
	k4Ca = kCa(vtd, r1, Cat);
	
	vts = Vs[itime] + timeStep * (25 * k1Vs / 48 - 55 * k2Vs / 24 + 35 * k3Vs / 48 + 15 * k4Vs / 8);
	n1 = n[itime] + timeStep * (25 * k1n / 48 - 55 * k2n / 24 + 35 * k3n / 48 + 15 * k4n / 8);
	h1 = h[itime] + timeStep * (25 * k1h / 48 - 55 * k2h / 24 + 35 * k3h / 48 + 15 * k4h / 8);
	vtd = Vd[itime] + timeStep * (25 * k1Vd / 48 - 55 * k2Vd / 24 + 35 * k3Vd / 48 + 15 * k4Vd / 8);
	r1 = r[itime] + timeStep * (25 * k1r / 48 - 55 * k2r / 24 + 35 * k3r / 48 + 15 * k4r / 8);
	c1 = c[itime] + timeStep * (25 * k1c / 48 - 55 * k2c / 24 + 35 * k3c / 48 + 15 * k4c / 8);
	Cat = Ca[itime] + timeStep * (25 * k1Ca / 48 - 55 * k2Ca / 24 + 35 * k3Ca / 48 + 15 * k4Ca / 8);
	t = time[itime] + 5 * timeStep / 6;

	k5Vs = kVs(vts, vtd, n1, h1, t);
	k5n = kn(vts, n1);
	k5h = kh(vts, h1);
	k5Vd = kVd(vts, vtd, r1, c1, Cat, t);
	k5r = kr(vtd, r1);
	k5c = kc(vtd, c1);
	k5Ca = kCa(vtd, r1, Cat);
	
	vts = Vs[itime] + timeStep * (3 * k1Vs / 20 - 11 * k2Vs / 24 - k3Vs / 8 + k4Vs / 2 + k5Vs / 10);
	n1 = n[itime] + timeStep * (3 * k1n / 20 - 11 * k2n / 24 - k3n / 8 + k4n / 2 + k5n / 10);
	h1 = h[itime] + timeStep * (3 * k1h / 20 - 11 * k2h / 24 - k3h / 8 + k4h / 2 + k5h / 10);
	vtd = Vd[itime] + timeStep * (3 * k1Vd / 20 - 11 * k2Vd / 24 - k3Vd / 8 + k4Vd / 2 + k5Vd / 10);
	r1 = r[itime] + timeStep * (3 * k1r / 20 - 11 * k2r / 24 - k3r / 8 + k4r / 2 + k5r / 10);
	c1 = c[itime] + timeStep * (3 * k1c / 20 - 11 * k2c / 24 - k3c / 8 + k4c / 2 + k5c / 10);
	Cat = Ca[itime] + timeStep * (3 * k1Ca / 20 - 11 * k2Ca / 24 - k3Ca / 8 + k4Ca / 2 + k5Ca / 10);
	t = time[itime] + timeStep / 6;

	k6Vs = kVs(vts, vtd, n1, h1, t);
	k6n = kn(vts, n1);
	k6h = kh(vts, h1);
	k6Vd = kVd(vts, vtd, r1, c1, Cat, t);
	k6r = kr(vtd, r1);
	k6c = kc(vtd, c1);
	k6Ca = kCa(vtd, r1, Cat);
	
	vts = Vs[itime] + timeStep * (- 261 * k1Vs / 260 + 33 * k2Vs / 13 + 43 * k3Vs / 156 - 118 * k4Vs / 39 + 32 * k5Vs / 195 + 80 * k6Vs / 39);
	n1 = n[itime] + timeStep * (- 261 * k1n / 260 + 33 * k2n / 13 + 43 * k3n / 156 - 118 * k4n / 39 + 32 * k5n / 195 + 80 * k6n / 39);
	h1 = h[itime] + timeStep * (- 261 * k1h / 260 + 33 * k2h / 13 + 43 * k3h / 156 - 118 * k4h / 39 + 32 * k5h / 195 + 80 * k6h / 39);
	vtd = Vd[itime] + timeStep * (- 261 * k1Vd / 260 + 33 * k2Vd / 13 + 43 * k3Vd / 156 - 118 * k4Vd / 39 + 32 * k5Vd / 195 + 80 * k6Vd / 39);
	r1 = r[itime] + timeStep * (- 261 * k1r / 260 + 33 * k2r / 13 + 43 * k3r / 156 - 118 * k4r / 39 + 32 * k5r / 195 + 80 * k6r / 39);
	c1 = c[itime] + timeStep * (- 261 * k1c / 260 + 33 * k2c / 13 + 43 * k3c / 156 - 118 * k4c / 39 + 32 * k5c / 195 + 80 * k6c / 39);
	Cat = Ca[itime] + timeStep * (- 261 * k1Ca / 260 + 33 * k2Ca / 13 + 43 * k3Ca / 156 - 118 * k4Ca / 39 + 32 * k5Ca / 195 + 80 * k6Ca / 39);
	t = time[itime] + timeStep;

	k7Vs = kVs(vts, vtd, n1, h1, t);
	k7n = kn(vts, n1);
	k7h = kh(vts, h1);
	k7Vd = kVd(vts, vtd, r1, c1, Cat, t);
	k7r = kr(vtd, r1);
	k7c = kc(vtd, c1);
	k7Ca = kCa(vtd, r1, Cat);

	//	update all values for next time point

	time[itime + 1] = time[itime] + timeStep;
	
	Vs[itime + 1] = Vs[itime] + timeStep * (13 * (k1Vs + k7Vs) / 200 + 11 * (k3Vs + k4Vs) / 40 + 4 * (k5Vs + k6Vs) / 25);
	n[itime + 1] = n[itime] + timeStep * (13 * (k1n + k7n) / 200 + 11 * (k3n + k4n) / 40 + 4 * (k5n + k6n) / 25);
	h[itime + 1] = h[itime] + timeStep * (13 * (k1h + k7h) / 200 + 11 * (k3h + k4h) / 40 + 4 * (k5h + k6h) / 25);
	
	Vd[itime + 1] = Vd[itime] + timeStep * (13 * (k1Vd + k7Vd) / 200 + 11 * (k3Vd + k4Vd) / 40 + 4 * (k5Vd + k6Vd) / 25);
	r[itime + 1] = r[itime] + timeStep * (13 * (k1r + k7r) / 200 + 11 * (k3r + k4r) / 40 + 4 * (k5r + k6r) / 25);
	c[itime + 1] = c[itime] + timeStep * (13 * (k1c + k7c) / 200 + 11 * (k3c + k4c) / 40 + 4 * (k5c + k6c) / 25);
	Ca[itime + 1] = Ca[itime] + timeStep * (13 * (k1Ca + k7Ca) / 200 + 11 * (k3Ca + k4Ca) / 40 + 4 * (k5Ca + k6Ca) / 25);
	
	Gexc_d[itime + 1] = Ge_d(time[itime + 1]);
    Gexc_s[itime + 1] = Ge_s(time[itime + 1]);

	Ginh_d[itime + 1] = Gi_d(time[itime + 1]);
	Ginh_s[itime + 1] = Gi_s(time[itime + 1]);
	
	Id[itime + 1] = IdExt(time[itime + 1]);
	Is[itime + 1] = IsExt(time[itime + 1]);
	E_gaba[itime + 1] = Ei;
	
	//	update internal time

	itime = itime + 1;
}


double HH2_final::Gi_s(double t){return Ginh_s[itime] * exp(-(t - time[itime]) / tInh);}
double HH2_final::Gi_d(double t){return Ginh_d[itime] * exp(-(t - time[itime]) / tInh);}
double HH2_final::Ge_d(double t){return Gexc_d[itime] * exp(-(t - time[itime]) / tExc);}
double HH2_final::Ge_s(double t){return Gexc_s[itime] * exp(-(t - time[itime]) / tExc);}
//double HH2_final::Gi_d(double t){return Ginh_d[itime];}

double HH2_final::kVs(double vs, double vd, double n, double h, double t)
{
	double m3, n4;

	m3 = mInf(vs) * mInf(vs) * mInf(vs);
	n4 = n * n * n * n;
	return (-GsL * (vs - EsL) - GsNa * m3 * h * (vs - EsNa) - GsK * n4 * (vs - EsK)
		 - Gi_s(t) * (vs - Ei) - Ge_s(t) * vs + (IsExt(t) + mu_soma) / As + (vd - vs) / (Rc * As)) / cm;
}
double HH2_final::kVd(double vs, double vd, double r, double c, double ca, double t)
{
	double r2;

	r2 = r * r;
	return (-GdL * (vd - EdL) - GdCa * r2 * (vd - EdCa) - GdCaK * (vd - EdK) * c / (1 + 6 / ca)
		- Ge_d(t) * vd - Gi_d(t) * (vd - Ei) + (IdExt(t) + mu_dend) / Ad + (vs - vd) / (Rc * Ad)) / cm;
}

double HH2_final::kn(double vs, double n){return (HH2_final::nInf(vs) - n) / HH2_final::tauN(vs);}
double HH2_final::kh(double vs, double h){return (HH2_final::hInf(vs) - h) / HH2_final::tauH(vs);}
double HH2_final::kr(double vd, double r){return (HH2_final::rInf(vd) - r) / HH2_final::tauR(vd);}
double HH2_final::kc(double vd, double c){return (HH2_final::cInf(vd) - c) / HH2_final::tauC(vd);}
double HH2_final::kCa(double vd, double r, double ca)
	{return (-0.1 * GdCa * r * r * (vd - EdCa) - 0.02 * ca);}
