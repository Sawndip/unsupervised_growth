#include "HHI_final.h"
#include "HH2_final.h"
#include "poisson_noise.h"
#include <fstream>
#include <iostream>
#include "exception.h"

using namespace std::placeholders;

HHI_final::HHI_final()
{
	// parameters of the model

	cm = 1;	// micro F / cm2
	A = 0.006;	// mm2

	Ena = 55;	// mV
	Ek = -80;
	El = -65;
	Ei = -75;

	gNa = 100;	// mS/cm2
	gKdr = 20;
	gKHT = 500;
	gL = 0.1;

	tExc = 2;	//	ms
	tInh = 5;

	threshold = -20;	//	mV

	// internal state

	itime = 0;
	Nspikes = 0;

	// noise
	generator = nullptr;

    // Poisson noise
	G_noise = 0.2;	//	maximum noise conductance
	poisson_noise = false;	//	turn on the noise

	lambda = 250; // intensity parameter for Poisson noise

	// white noise
	bin_size = 1;
	stored = 0;
    white_noise = false;

	// external current
	Iext = std::bind(&HHI_final::I_default, this, _1); // set external current to default function
}

HHI_final::HHI_final(DDfunction f) : HHI_final()
{
	Iext = f;
}

void HHI_final::print_param()
{
	std::cout << "Cm = " << cm << " microF / cm^2" << std::endl;

	std::cout << "EL = " << El << " mV" << std::endl;
	std::cout << "EK = " << Ek << " mV" << std::endl;
	std::cout << "ENa = " << Ena << " mV" << std::endl;
	std::cout << "Ei = " << Ei << " mV" << std::endl;

	std::cout << "gL = " << gL << " mS / cm^2" << std::endl;
	std::cout << "gKdr = " << gKdr << " mS / cm^2" << std::endl;
	std::cout << "gKHT = " << gKHT << " mS / cm^2" << std::endl;
	std::cout << "gNa = " << gNa << " mS / cm^2" << std::endl;

	std::cout << "tExc = " << tExc << " ms" << std::endl;
	std::cout << "tInh = " << tInh << " ms" << std::endl;
}

bool HHI_final::get_fired()
{
	return fired;
}

void HHI_final::print_targets()
{
	printf("Connections to HVC(RA) neurons:\n");
	for (int i = 0; i < targets.size(); i++)
	{
		printf("Target neuron %d \t Synaptic conductance %f",
			targetsID[i], targetsG[i]);
	}

}

vector<int> HHI_final::get_targetsID()
{
	return targetsID;
}

vector<double> HHI_final::get_targetsG()
{
	return targetsG;
}

void HHI_final::writeToFile(const char * filename)
{
	std::ofstream outputFile;

	outputFile.open(filename, std::ios::out | std::ios::binary); //	open file to write binary data

	// write parameters of the model and number of data points

	outputFile.write(reinterpret_cast<char *>(&size), sizeof(size)); // total number of data points
	outputFile.write(reinterpret_cast<char *>(&Nspikes), sizeof(Nspikes));
	outputFile.write(reinterpret_cast<char *>(&timeStep), sizeof(timeStep));

	//	Now let's write time data

	for (vector<double>::size_type i = 0; i < voltage.size(); ++i)
	{
		outputFile.write(reinterpret_cast<char *>(&time[i]), sizeof(time[i]));
		outputFile.write(reinterpret_cast<char *>(&voltage[i]), sizeof(voltage[i]));
		outputFile.write(reinterpret_cast<char *>(&I[i]), sizeof(I[i]));
		outputFile.write(reinterpret_cast<char *>(&n[i]), sizeof(n[i]));
		outputFile.write(reinterpret_cast<char *>(&m[i]), sizeof(m[i]));
		outputFile.write(reinterpret_cast<char *>(&h[i]), sizeof(h[i]));
		outputFile.write(reinterpret_cast<char *>(&w[i]), sizeof(w[i]));
		outputFile.write(reinterpret_cast<char *>(&Gexc[i]), sizeof(Gexc[i]));
		outputFile.write(reinterpret_cast<char *>(&Ginh[i]), sizeof(Ginh[i]));
		outputFile.write(reinterpret_cast<char *>(&flag[i]), sizeof(flag[i]));
	}

	outputFile.close(); //	close file

}

int HHI_final::get_spike_number()
{
	return Nspikes;
}

vector<int> HHI_final::get_activity()
{
	return flag;
}

void HHI_final::set_injected_current(DDfunction f)
{
	Iext = f;
}

void HHI_final::set_noise_generator(Poisson_noise* g)
{
	generator = g;
}

void HHI_final::set_white_noise(double m, double s)
{
    try
    {
        if (generator == nullptr)
        {
            throw NoGenerator("Noise generator is not set for HHI neuron!\n");
        }
        else
        {
            white_noise = true;
            
            mu = m;
            sigma = s;
            
            Iext = std::bind(&HHI_final::I_white_noise, this, _1);
        }
    }

    catch (NoGenerator const& e)
    {
        std::cerr << "NoGenerator Exception: " << e.what() << std::endl;
    }
}

void HHI_final::set_poisson_noise()
{
    poisson_noise = true;
            
    this->initialize_noise(noise_exc);
	this->initialize_noise(noise_inh);
}

void HHI_final::set_no_white_noise()
{
    white_noise = false;
}

void HHI_final::set_no_poisson_noise()
{
	poisson_noise = false;
}

void HHI_final::set_dynamics(double interval, double tS)
{

	timeStep = tS;
	size = (int) round(interval / timeStep) + 1;

	//	Resize all arrays

	time.resize(size);
	voltage.resize(size);
	I.resize(size);
	n.resize(size);
	m.resize(size);
	h.resize(size);
	w.resize(size);
	Gexc.resize(size);
	Ginh.resize(size);
	flag.resize(size);

	time[0] = 0;
	voltage[0] = -66;
	n[0] = 0.125;
	m[0] = 0;
	h[0] = 0.99;
	w[0] = 0;
	Gexc[0] = 0;
	Ginh[0] = 0;
	flag[0] = 0;
	I[0] = Iext(time[0]);

	//	set up noise
    if (poisson_noise)
    {
	    this->initialize_noise(noise_exc);
	    this->initialize_noise(noise_inh);
    }
}

void HHI_final::set_to_rest()
{

	itime = 0;

	time[0] = 0;
	voltage[0] = -66;
	n[0] = 0.125;
	m[0] = 0;
	h[0] = 0.99;
	w[0] = 0;
	Gexc[0] = 0;
	Ginh[0] = 0;
	flag[0] = 0;
	I[0] = Iext(time[0]);

	//	set up noise
    if (poisson_noise)
    {
	    this->initialize_noise(noise_exc);
	    this->initialize_noise(noise_inh);
    }
}

void HHI_final::initialize_noise(int& noise_time)
{
    try
    {
        if (generator == nullptr)
        {
            throw NoGenerator("Noise generator is not set for HHI neuron!\n");
        }
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

void HHI_final::set_target(HH2_final* target, int ID, double G)
{
	targets.push_back(target);
	targetsID.push_back(ID);
	targetsG.push_back(G);
}

void HHI_final::state_noise_check()
{
	this->state_check();

	if (poisson_noise)
	{
		this->noise_check(Gexc[itime], noise_exc);
		this->noise_check(Ginh[itime], noise_inh);
	}

}

void HHI_final::state_check()
{
	if ((flag[itime] == 1) && (voltage[itime] < threshold))
	{
		flag[itime + 1] = 0;
		Nspikes = Nspikes + 1;
		fired = true;
	}
	else
	{
		fired = false;
		//	check if we should change the state of neuron (voltage crossed the threshold)
		if ((flag[itime] == 0) && (voltage[itime] > threshold))
			flag[itime + 1] = 1;
		else
			flag[itime + 1] = flag[itime];	//	otherwise state hasn't changed
	}

}

void HHI_final::noise_check(double& G, int& noise_time)
{
	if (itime == noise_time)
	{
		G += generator->random(G_noise);

		int random = round(1000 * generator->get_spike_time(lambda) / timeStep);
		while (random == 0)
		{
				random = (int) round(1000 * generator->get_spike_time(lambda) / timeStep);
				G += generator->random(G_noise);
		}
		noise_time = noise_time + random;
	}

}

void HHI_final::R4_step_with_target_update()
{
	this->state_noise_check();

	if (this->get_fired())
		this->postsyn_update();
	
    this->Runge4_step();

}

void HHI_final::R4_step_no_target_update()
{
	this->state_noise_check();
	//printf("I[%d] = %f\n", itime, I[itime]);
	this->Runge4_step();
}


void HHI_final::Runge4_step()
{
	double m1, n1, h1, w1, m3, n4;
	double v, t;
	double k1V, k2V, k3V, k4V;
	double k1n, k2n, k3n, k4n;
	double k1m, k2m, k3m, k4m;
	double k1h, k2h, k3h, k4h;
	double k1w, k2w, k3w, k4w;

	m1 = m[itime];
	n1 = n[itime];
	h1 = h[itime];
	w1 = w[itime];
	v = voltage[itime];
	t = time[itime];
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k1V = kV(v, t, h1, w1, m3, n4);
	k1n = kn(v, n1);
	k1m = km(v, m1);
	k1h = kh(v, h1);
	k1w = kw(v, w1);

	m1 = m[itime] + timeStep * k1m / 3;
	n1 = n[itime] + timeStep * k1n / 3;
	h1 = h[itime] + timeStep * k1h / 3;
	w1 = w[itime] + timeStep * k1w / 3;
	v = voltage[itime] + timeStep * k1V / 3;
	t = time[itime] + timeStep / 3;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k2V = kV(v, t, h1, w1, m3, n4);
	k2n = kn(v, n1);
	k2m = km(v, m1);
	k2h = kh(v, h1);
	k2w = kw(v, w1);

	m1 = m[itime] + timeStep * (-k1m / 3 + k2m);
	n1 = n[itime] + timeStep * (-k1n / 3 + k2n);
	h1 = h[itime] + timeStep * (-k1h / 3 + k2h);
	w1 = w[itime] + timeStep * (-k1w / 3 + k2w);
	v = voltage[itime] + timeStep * (-k1V / 3 + k2V);
	t = time[itime] + 2 * timeStep / 3;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k3V = kV(v, t, h1, w1, m3, n4);
	k3n = kn(v, n1);
	k3m = km(v, m1);
	k3h = kh(v, h1);
	k3w = kw(v, w1);

	m1 = m[itime] + timeStep * (k1m - k2m + k3m);
	n1 = n[itime] + timeStep * (k1n - k2n + k3n);
	h1 = h[itime] + timeStep * (k1h - k2h + k3h);
	w1 = w[itime] + timeStep * (k1w - k2w + k3w);
	v = voltage[itime] + timeStep * (k1V - k2V + k3V);
	t = time[itime] + timeStep;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k4V = kV(v, t, h1, w1, m3, n4);
	k4n = kn(v, n1);
	k4m = km(v, m1);
	k4h = kh(v, h1);
	k4w = kw(v, w1);

	//	update all values for next time point

	voltage[itime + 1] = voltage[itime] + timeStep * (k1V + 3 * k2V + 3 * k3V + k4V) / 8;
	n[itime + 1] = n[itime] + timeStep * (k1n + 3 * k2n + 3 * k3n + k4n) / 8;
	m[itime + 1] = m[itime] + timeStep * (k1m + 3 * k2m + 3 * k3m + k4m) / 8;
	h[itime + 1] = h[itime] + timeStep * (k1h + 3 * k2h + 3 * k3h + k4h) / 8;
	w[itime + 1] = w[itime] + timeStep * (k1w + 3 * k2w + 3 * k3w + k4w) / 8;
	time[itime + 1] = time[itime] + timeStep;
	Gexc[itime + 1] = Ge(time[itime + 1]);
	Ginh[itime + 1] = Gi(time[itime + 1]);
	I[itime + 1] = Iext(time[itime + 1]);
	//	update internal time

	itime = itime + 1;
}


void HHI_final::raiseE(double G)
{
	Gexc[itime] = Gexc[itime] + G;
}

void HHI_final::raiseI(double G)
{
	Ginh[itime] = Ginh[itime] + G;
}


void HHI_final::postsyn_update()
{
	for (unsigned i = 0; i < targets.size(); i++)
	{
		//targets[i]->raiseI(targetsG[i]);
	}
}

double HHI_final::I_white_noise(double t)
{
	if (itime % (int) round((bin_size / timeStep)) == 0)
	{
		//printf("(itime+1)*timeStep = %f\tbin size * count = %f\n", (itime + 1)*timeStep, bin_size*count);
		//printf("bin_size = %f\tcount = %d\n", bin_size, count);

		stored = mu + sigma * generator->normal_distribution();
		//printf("stored = %f\n", stored);
		//printf("Sample new\titime = %d\n", itime);
		return stored;
	}
	else
	{
		//printf("stored = %f\n", stored);
		//printf("Sample old\titime = %d\n", itime);
		return stored;

	}
}


double HHI_final::Ge(double t){return Gexc[itime] * exp(-(t - time[itime]) / tExc);}
double HHI_final::Gi(double t){return Ginh[itime] * exp(-(t - time[itime]) / tInh);}
double HHI_final::kV(double v, double t, double h, double w, double m3, double n4){
	return (-gL * (v - El) - gNa * h * m3 * (v - Ena) - gKdr * n4 * (v - Ek)
		- gKHT * w * (v - Ek) - Ge(t) * v - Gi(t) * (v - Ei) + Iext(t)) / cm;}
double HHI_final::kn(double v, double n){return HHI_final::an(v)*(1 - n) - HHI_final::bn(v)*n;}
double HHI_final::km(double v, double m){return HHI_final::am(v)*(1 - m) - HHI_final::bm(v)*m;}
double HHI_final::kh(double v, double h){return HHI_final::ah(v)*(1 - h) - HHI_final::bh(v)*h;}
double HHI_final::kw(double v, double w){return (HHI_final::wInf(v) - w) / HHI_final::tauW(v);}
double HHI_final::kexc(double gexc){return -gexc / tExc;}
double HHI_final::kinh(double ginh){return -ginh / tInh;}

