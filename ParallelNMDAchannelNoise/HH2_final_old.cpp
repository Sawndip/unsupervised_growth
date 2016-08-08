#include "poisson_noise.h"
#include "HH2_final.h"
#include "HHI_final.h"
#include <iostream>
#include <fstream>

using namespace std::placeholders;

HH2_final::HH2_final()
{
	// parameters of the model
	cm = 1;
	Rc = 0.055;
	As = 50;
	GsL = 0.1;
	GsNa = 60;
	GsK = 8;
	EsL = -80;
	EsNa = 55;
	EsK = -90;
	Ad = 100;
	GdL = 0.1;
	GdCa = 55;
	GdCaK = 150;
	EdL = -80;
	EdCa = 120;
	EdK = -90;
	Ei = -80;
	tExc = 5;
	tInh = 5;

	threshold_soma = 0;
	threshold_dendrite = -20;

	//internal state
	itime = 0;
	Nspikes = 0;
	fired_soma = false;
    fired_dendrite = false;

	// noise
	Gs_noise_exc = 0.025;
	Gs_noise_inh = 0.025;
	Gd_noise_exc = 0.065;
	Gd_noise_inh = 0.065;

	noise = true;
	lambda = 100;

	// white noise
	bin_size = 0.25;
	count = 1;

	// external current
	IsExt = std::bind(&HH2_final::Is_default, this, _1);
	IdExt = std::bind(&HH2_final::Id_default, this, _1);
}

HH2_final::HH2_final(bool white_noise_soma, bool white_noise_dendrite) : HH2_final()
{
	if (white_noise_soma)
		IsExt = std::bind(&HH2_final::I_white_noise, this, _1);

	if (white_noise_dendrite)
		IdExt = std::bind(&HH2_final::I_white_noise, this, _1);
}

HH2_final::HH2_final(DDfunction fs, DDfunction fd) : HH2_final()
{
	IsExt = fs;
	IdExt = fd;
}

void HH2_final::set_soma_current(DDfunction fs){IsExt = fs;}
void HH2_final::set_dend_current(DDfunction fd){IdExt = fd;}

void HH2_final::set_no_noise()
{
	noise = false;
}

void HH2_final::initialize_noise(int& noise_time)
{
	noise_time = (int) round(1000 * generator->get_spike_time(lambda) / timeStep);
	while (noise_time == 0)
		noise_time = (int) round(1000 * generator->get_spike_time(lambda) / timeStep);
}

void HH2_final::set_dynamics(double interval, double tS)
{
	timeStep = tS;
	size = (int) round(interval / timeStep);

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

    E_gaba[0] = Ei;
    time[0] = 0;
	Vs[0] = -80;
	Vd[0] = -80;
	n[0] = 0.0110128650344;
	h[0] = 0.993284481138;
	r[0] = 0.000554291044048;
	c[0] = 0.0000025762773059;
	Ca[0] = 0.0168;
	Ginh_s[0] = 0;
	Gexc_s[0] = 0;
	Ginh_d[0] = 0;
	Gexc_d[0] = 0;

	flag_soma[0] = 0;
    flag_dend[0] = 0;

	Is[0] = IsExt(time[0]);
	Id[0] = IdExt(time[0]);

	//	initialize up noise

	this->initialize_noise(noise_exc_soma);
	this->initialize_noise(noise_inh_soma);
	this->initialize_noise(noise_exc_dend);
	this->initialize_noise(noise_inh_dend);

}

void HH2_final::set_noise_generator(Poisson_noise* g)
{
	generator = g;
}

void HH2_final::set_white_noise_distribution(double mu, double sigma)
{
	generator->set_normal_distribution(mu, sigma);
}

void HH2_final::set_white_noise_soma()
{
	IsExt = std::bind(&HH2_final::I_white_noise, this, _1);
}

void HH2_final::set_white_noise_dend()
{
	IdExt = std::bind(&HH2_final::I_white_noise, this, _1);
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

double HH2_final::get_spike_time_soma()
{
	return spike_time_soma;
}

bool HH2_final::get_fired_soma()
{
	return fired_soma;
}

bool HH2_final::get_fired_dend()
{
    return fired_dendrite;
}

double HH2_final::get_spike_time_dend()
{
    return spike_time_dend;
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
	Gexc_d[itime] = Gexc_d[itime] + G;
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

	sampled_this_iteration = false;
	this->Runge4_step();
}

void HH2_final::R4_step_no_target_update()
{
	this->state_noise_check();
	sampled_this_iteration = false;
	this->Runge4_step();
}

void HH2_final::state_noise_check()
{
	this->state_check();

	if (noise)
	{
		this->noise_check(Gexc_s[itime], Gs_noise_exc, noise_exc_soma);
		this->noise_check(Ginh_s[itime], Gs_noise_inh, noise_inh_soma);
		this->noise_check(Gexc_d[itime], Gd_noise_exc, noise_exc_dend);
		this->noise_check(Ginh_d[itime], Gd_noise_inh, noise_inh_dend);
	}
}

void HH2_final::state_check()
{
	if ((flag_soma[itime] == 1) && (Vs[itime] < threshold_soma))
	{
		spike_time_soma = time[itime];
		flag_soma[itime + 1] = 0;
		fired_soma = true;
		Nspikes = Nspikes + 1;
	}
	else
	{
		//	check if we should change the state of neuron (voltage crossed the threshold)
		if ((flag_soma[itime] == 0) && (Vs[itime] > threshold_soma))
			flag_soma[itime + 1] = 1;
		else
		{
			fired_soma = false;
			flag_soma[itime + 1] = flag_soma[itime];	//	otherwise state hasn't changed
		}
	}

	if ((flag_dend[itime] == 1) && (Vd[itime] < threshold_dendrite))
	{
		spike_time_dend = time[itime];
		flag_dend[itime + 1] = 0;
		fired_dendrite = true;
		//Nspikes = Nspikes + 1;
	}
	else
	{
		//	check if we should change the state of neuron (voltage crossed the threshold)
		if ((flag_dend[itime] == 0) && (Vd[itime] > threshold_dendrite))
			flag_dend[itime + 1] = 1;
		else
		{
			fired_dendrite = false;
			flag_dend[itime + 1] = flag_dend[itime];	//	otherwise state hasn't changed
		}
	}
}

void HH2_final::noise_check(double& G, double G_noise, int& noise_time)
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

double HH2_final::I_white_noise(double t)
{
	if (((itime + 1) * timeStep > bin_size * count)&&(!sampled_this_iteration))
	{
		count++;
		stored = generator->normal_distribution();
		sampled_this_iteration = true;
		return stored;
	}
	else
		return stored;
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
	Gexc_s[itime + 1] = Ge_s(time[itime + 1]);
	Ginh_s[itime + 1] = Gi_s(time[itime + 1]);
	Vd[itime + 1] = Vd[itime] + timeStep * (k1Vd + 3 * k2Vd + 3 * k3Vd + k4Vd) / 8;
	r[itime + 1] = r[itime] + timeStep * (k1r + 3 * k2r + 3 * k3r + k4r) / 8;
	c[itime + 1] = c[itime] + timeStep * (k1c + 3 * k2c + 3 * k3c + k4c) / 8;
	Ca[itime + 1] = Ca[itime] + timeStep * (k1Ca + 3 * k2Ca + 3 * k3Ca + k4Ca) / 8;
	Gexc_d[itime + 1] = Ge_d(time[itime + 1]);
	Ginh_d[itime + 1] = Gi_d(time[itime + 1]);

	Id[itime + 1] = IdExt(time[itime + 1]);
	Is[itime + 1] = IsExt(time[itime + 1]);
    E_gaba[itime + 1] = Ei;
	itime = itime + 1;
}

double HH2_final::Ge_s(double t){return Gexc_s[itime] * exp(-(t - time[itime]) / tExc);}
double HH2_final::Gi_s(double t){return Ginh_s[itime] * exp(-(t - time[itime]) / tInh);}
double HH2_final::Ge_d(double t){return Gexc_d[itime] * exp(-(t - time[itime]) / tExc);}
double HH2_final::Gi_d(double t){return Ginh_d[itime] * exp(-(t - time[itime]) / tInh);}

double HH2_final::kVs(double vs, double vd,
	double n, double h, double t)
{
	double m3, n4;

	m3 = mInf(vs) * mInf(vs) * mInf(vs);
	n4 = n * n * n * n;
	return (-GsL * (vs - EsL) - GsNa * m3 * h * (vs - EsNa) - GsK * n4 * (vs - EsK)
		- Ge_s(t) * vs - Gi_s(t) * (vs - Ei) + IsExt(t) / As + (vd - vs) / (Rc * As)) / cm;
}
double HH2_final::kVd(double vs, double vd,
	double r, double c, double ca, double t)
{
	double r2;

	r2 = r * r;
	return (-GdL * (vd - EdL) - GdCa * r2 * (vd - EdCa) - GdCaK * (vd - EdK) * c / (1 + 6 / ca)
		- Ge_d(t) * vd - Gi_d(t) * (vd - Ei) + IdExt(t) / Ad + (vs - vd) / (Rc * Ad)) / cm;
}

double HH2_final::kn(double vs, double n){return (nInf(vs) - n) / tauN(vs);}
double HH2_final::kh(double vs, double h){return (hInf(vs) - h) / tauH(vs);}
double HH2_final::kr(double vd, double r){return (rInf(vd) - r) / tauR(vd);}
double HH2_final::kc(double vd, double c){return (cInf(vd) - c) / tauC(vd);}

double HH2_final::kCa(double vd, double r, double ca)
	{return (-0.1 * GdCa * r * r * (vd - EdCa) - 0.02 * ca);}
