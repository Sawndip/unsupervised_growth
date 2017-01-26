#include "HHI_final.h"
#include "HH2_final.h"
#include "poisson_noise.h"
#include <fstream>
#include <iostream>
#include "exception.h"
#include <algorithm>
#include <cmath>

using namespace std::placeholders;

HHI_final::HHI_final() : mu(0.0), sigma(0.0)
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
	injected_current = false;

	// noise
	generator = nullptr;

    // Poisson noise
	G_noise = 0.2;	//	maximum noise conductance
	poisson_noise = false;	//	turn on the noise

	lambda = 250; // intensity parameter for Poisson noise

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
	I_injected = f;
	injected_current = true;
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
            mu = m;
            sigma = s;
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
    mu = 0.0;
    sigma = 0.0;
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
	
	fm.resize(size);
	fn.resize(size);
	fh.resize(size);
	fw.resize(size);
	fv.resize(size);

	I[0] = I_total(time[0]);

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
	
	I[0] = I_total(time[0]);

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
        	noise_time = static_cast<int>(1000 * generator->get_spike_time(lambda) / timeStep);
	        
            while (noise_time == 0)
		        noise_time = static_cast<int>(1000 * generator->get_spike_time(lambda) / timeStep);
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

		int random = static_cast<int>(1000 * generator->get_spike_time(lambda) / timeStep);
		while (random == 0)
		{
				random = static_cast<int>(1000 * generator->get_spike_time(lambda) / timeStep);
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

void HHI_final::Euler_step_no_target_update()
{
	this->state_noise_check();
	this->Euler_step();
}

void HHI_final::iEuler_step_no_target_update()
{
	this->state_noise_check();
	this->iEuler_step();
}

void HHI_final::iTrapezoid_step_no_target_update()
{
	this->state_noise_check();
	this->iTrapezoid_step();
}

void HHI_final::R4_step_no_target_update()
{
	this->state_noise_check();
	this->Runge4_step();
}

void HHI_final::DP5_step_no_target_update()
{
	this->state_noise_check();
	this->DP5_step();
}

void HHI_final::PC4_step_no_target_update()
{
	this->state_noise_check();
	this->AB4_step();
	itime -= 1;
	this->AM4_step();
}

void HHI_final::AB4_step_no_target_update()
{
	this->state_noise_check();
	this->AB4_step();
}

void HHI_final::R6_step_no_target_update()
{
	this->state_noise_check();
	this->Runge6_step();
}

void HHI_final::BDF4_step_no_target_update()
{
	this->state_noise_check();
	this->BDF4_step();
}

void HHI_final::BDF6_step_no_target_update()
{
	this->state_noise_check();
	this->BDF6_step();
}

void HHI_final::DP8_step_no_target_update()
{
	this->state_noise_check();
	this->DP8_step();
}

void HHI_final::BDF6_step()
{
	
	// if there are enough initial steps
	if (itime >= 5)
	{
		time[itime + 1] = time[itime] + timeStep;
		Gexc[itime + 1] = Ge(time[itime + 1]);
		Ginh[itime + 1] = Gi(time[itime + 1]);
		I[itime + 1] = I_total(time[itime + 1]);
		
		double eps = 1e-4;
		double vmax = 100.0;
		int Nmax = 100;

		DDfunction f = std::bind(&HHI_final::f_bdf6, this, _1);
		voltage[itime + 1] = bisection(f, voltage[itime], vmax, eps, Nmax);

		m[itime + 1] = m_bdf6(voltage[itime + 1]);
		n[itime + 1] = n_bdf6(voltage[itime + 1]);
		h[itime + 1] = h_bdf6(voltage[itime + 1]);
		w[itime + 1] = w_bdf6(voltage[itime + 1]);

		itime += 1;
	}
	else // if not make Runge Kutta order 4 step
		//Runge4_step();
		//iEuler_step();
		DP8_step();
}


void HHI_final::AB4_step()
{
	// if there are enough initial steps
	if (itime >= 3)
	{
		double n4 = n[itime] * n[itime] * n[itime] * n[itime];
		double m3 = m[itime] * m[itime] * m[itime];
		
		fv[itime] = kV(voltage[itime], time[itime], h[itime], w[itime], m3, n4); // value of right handside of equation for voltage 
		fn[itime] = kn(voltage[itime], n[itime]); // value of right handside of equation for n
		fm[itime] = km(voltage[itime], m[itime]); // value of right handside of equation for m
		fh[itime] = kh(voltage[itime], h[itime]); // value of right handside of equation for h
		fw[itime] = kw(voltage[itime], w[itime]); // value of right handside of equation for w
		
		time[itime + 1] = time[itime] + timeStep;
		Gexc[itime + 1] = Ge(time[itime + 1]);
		Ginh[itime + 1] = Gi(time[itime + 1]);
		I[itime + 1] = I_total(time[itime + 1]);
		
		// predicted values
		voltage[itime + 1] = voltage[itime] + timeStep * (55 * fv[itime] - 59 * fv[itime - 1] + 37 * fv[itime - 2] - 9 * fv[itime - 3]) / 24;
		n[itime + 1] = n[itime] + timeStep * (55 * fn[itime] - 59 * fn[itime - 1] + 37 * fn[itime - 2] - 9 * fn[itime - 3]) / 24;
		m[itime + 1] = m[itime] + timeStep * (55 * fm[itime] - 59 * fm[itime - 1] + 37 * fm[itime - 2] - 9 * fm[itime - 3]) / 24;
		h[itime + 1] = h[itime] + timeStep * (55 * fh[itime] - 59 * fh[itime - 1] + 37 * fh[itime - 2] - 9 * fh[itime - 3]) / 24;
		w[itime + 1] = w[itime] + timeStep * (55 * fw[itime] - 59 * fw[itime - 1] + 37 * fw[itime - 2] - 9 * fw[itime - 3]) / 24;
	
		itime += 1;
	}
	else // if not make i Euler step
	{	
		double n4 = n[itime] * n[itime] * n[itime] * n[itime];
		double m3 = m[itime] * m[itime] * m[itime];

		fv[itime] = kV(voltage[itime], time[itime], h[itime], w[itime], m3, n4); // value of right handside of equation for voltage 
		fn[itime] = kn(voltage[itime], n[itime]); // value of right handside of equation for n
		fm[itime] = km(voltage[itime], m[itime]); // value of right handside of equation for m
		fh[itime] = kh(voltage[itime], h[itime]); // value of right handside of equation for h
		fw[itime] = kw(voltage[itime], w[itime]); // value of right handside of equation for w
		
		//iEuler_step();
		Runge4_step();
	}
}

void HHI_final::AM4_step()
{
	// if there are enough initial steps
	if (itime >= 3)
	{
		// predicted values of right handside	
		double n4 = n[itime + 1] * n[itime + 1] * n[itime + 1] * n[itime + 1];
		double m3 = m[itime + 1] * m[itime + 1] * m[itime + 1];
		
		fv[itime + 1] = kV(voltage[itime + 1], time[itime + 1], h[itime + 1], w[itime + 1], m3, n4); // value of right handside of equation for voltage 
		fn[itime + 1] = kn(voltage[itime + 1], n[itime + 1]); // value of right handside of equation for n
		fm[itime + 1] = km(voltage[itime + 1], m[itime + 1]); // value of right handside of equation for m
		fh[itime + 1] = kh(voltage[itime + 1], h[itime + 1]); // value of right handside of equation for h
		fw[itime + 1] = kw(voltage[itime + 1], w[itime + 1]); // value of right handside of equation for w

		// corrected values
		voltage[itime + 1] = voltage[itime] + timeStep * (9 * fv[itime + 1] + 19 * fv[itime] - 5 * fv[itime - 1] + fv[itime - 2]) / 24;
		n[itime + 1] = n[itime] + timeStep * (9 * fn[itime + 1] + 19 * fn[itime] - 5 * fn[itime - 1] + fn[itime - 2]) / 24;
		m[itime + 1] = m[itime] + timeStep * (9 * fm[itime + 1] + 19 * fm[itime] - 5 * fm[itime - 1] + fm[itime - 2]) / 24;
		h[itime + 1] = h[itime] + timeStep * (9 * fh[itime + 1] + 19 * fh[itime] - 5 * fh[itime - 1] + fh[itime - 2]) / 24;
		w[itime + 1] = w[itime] + timeStep * (9 * fw[itime + 1] + 19 * fw[itime] - 5 * fw[itime - 1] + fw[itime - 2]) / 24;
		
	}
	itime += 1;
}

void HHI_final::BDF4_step()
{
	
	// if there are enough initial steps
	if (itime >= 3)
	{
		time[itime + 1] = time[itime] + timeStep;
		Gexc[itime + 1] = Ge(time[itime + 1]);
		Ginh[itime + 1] = Gi(time[itime + 1]);
		I[itime + 1] = I_total(time[itime + 1]);
		
		double eps = 1e-4;
		double vmax = 100.0;
		int Nmax = 100;

		DDfunction f = std::bind(&HHI_final::f_bdf4, this, _1);
		voltage[itime + 1] = bisection(f, voltage[itime], vmax, eps, Nmax);

		m[itime + 1] = m_bdf4(voltage[itime + 1]);
		n[itime + 1] = n_bdf4(voltage[itime + 1]);
		h[itime + 1] = h_bdf4(voltage[itime + 1]);
		w[itime + 1] = w_bdf4(voltage[itime + 1]);

		itime += 1;
	}
	else // if not make Runge Kutta order 4 step
		//Runge4_step();
		//iEuler_step();
		DP8_step();
}

void HHI_final::Euler_step()
{
	voltage[itime + 1] = voltage[itime] + timeStep * ( (-gL * (voltage[itime] - El) - gNa * h[itime] * m[itime] * m[itime] * m[itime] * (voltage[itime] - Ena) - 
						gKdr * n[itime] * n[itime] * n[itime] * n[itime] * (voltage[itime] - Ek) - gKHT * w[itime] * (voltage[itime] - Ek) - 
						Gexc[itime] * voltage[itime] - Ginh[itime] * (voltage[itime] - Ei) + I[itime]) / cm);

	m[itime + 1] = m[itime] + timeStep * (am(voltage[itime]) * (1 - m[itime]) - bm(voltage[itime]) * m[itime]);
	n[itime + 1] = n[itime] + timeStep * (an(voltage[itime]) * (1 - n[itime]) - bn(voltage[itime]) * n[itime]);
	h[itime + 1] = h[itime] + timeStep * (ah(voltage[itime]) * (1 - h[itime]) - bh(voltage[itime]) * h[itime]);
	w[itime + 1] = w[itime] + timeStep * (wInf(voltage[itime]) - w[itime]) / tauW(voltage[itime]);


	time[itime + 1] = time[itime] + timeStep;
	Gexc[itime + 1] = Ge(time[itime + 1]);
	Ginh[itime + 1] = Gi(time[itime + 1]);
	I[itime + 1] = I_total(time[itime + 1]);
	
	//	update internal time

	itime = itime + 1;

}

void HHI_final::iEuler_step()
{

	time[itime + 1] = time[itime] + timeStep;
	Gexc[itime + 1] = Ge(time[itime + 1]);
	Ginh[itime + 1] = Gi(time[itime + 1]);
	I[itime + 1] = I_total(time[itime + 1]);
		
	double eps = 1e-4;
	double vmax = 100.0;
	int Nmax = 100;

	DDfunction f = std::bind(&HHI_final::f_iEuler, this, _1);
	voltage[itime + 1] = bisection(f, voltage[itime], vmax, eps, Nmax);

	m[itime + 1] = m_iEuler(voltage[itime + 1]);
	n[itime + 1] = n_iEuler(voltage[itime + 1]);
	h[itime + 1] = h_iEuler(voltage[itime + 1]);
	w[itime + 1] = w_iEuler(voltage[itime + 1]);

	itime += 1;

}

void HHI_final::iTrapezoid_step()
{
	time[itime + 1] = time[itime] + timeStep;
	Gexc[itime + 1] = Ge(time[itime + 1]);
	Ginh[itime + 1] = Gi(time[itime + 1]);
	I[itime + 1] = I_total(time[itime + 1]);
		
	double eps = 1e-4;
	double vmax = 100.0;
	int Nmax = 100;


	double n4 = n[itime] * n[itime] * n[itime] * n[itime];
	double m3 = m[itime] * m[itime] * m[itime];

	double fv = kV(voltage[itime], time[itime], h[itime], w[itime], m3, n4); // value of right handside of equation for voltage 
	double fn = kn(voltage[itime], n[itime]); // value of right handside of equation for n
	double fm = km(voltage[itime], m[itime]); // value of right handside of equation for m
	double fh = kh(voltage[itime], h[itime]); // value of right handside of equation for h
	double fw = kw(voltage[itime], w[itime]); // value of right handside of equation for w

	DDfunction f = std::bind(&HHI_final::f_iTrapezoid, this, fv, fn, fm, fh, fw, _1);
	
	voltage[itime + 1] = bisection(f, voltage[itime], vmax, eps, Nmax);

	m[itime + 1] = m_iTrapezoid(fm, voltage[itime + 1]);
	n[itime + 1] = n_iTrapezoid(fn, voltage[itime + 1]);
	h[itime + 1] = h_iTrapezoid(fh, voltage[itime + 1]);
	w[itime + 1] = w_iTrapezoid(fw, voltage[itime + 1]);

	itime += 1;

}

void HHI_final::Runge6_step()
{
	double m1, n1, h1, w1, m3, n4;
	double v, t;
	double k1V, k2V, k3V, k4V, k5V, k6V, k7V;
	double k1n, k2n, k3n, k4n, k5n, k6n, k7n;
	double k1m, k2m, k3m, k4m, k5m, k6m, k7m;
	double k1h, k2h, k3h, k4h, k5h, k6h, k7h;
	double k1w, k2w, k3w, k4w, k5w, k6w, k7w;

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

	m1 = m[itime] + 2 * timeStep * k2m / 3;
	n1 = n[itime] + 2 * timeStep * k2n / 3;
	h1 = h[itime] + 2 * timeStep * k2h / 3;
	w1 = w[itime] + 2 * timeStep * k2w / 3;
	v = voltage[itime] + 2 * timeStep * k2V / 2;
	t = time[itime] + 2 * timeStep / 3;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k3V = kV(v, t, h1, w1, m3, n4);
	k3n = kn(v, n1);
	k3m = km(v, m1);
	k3h = kh(v, h1);
	k3w = kw(v, w1);

	m1 = m[itime] + timeStep * (k1m / 12 + k2m / 3 - k3m / 12);
	n1 = n[itime] + timeStep * (k1n / 12 + k2n / 3 - k3n / 12);
	h1 = h[itime] + timeStep * (k1h / 12 + k2h / 3 - k3h / 12);
	w1 = w[itime] + timeStep * (k1w / 12 + k2w / 3 - k3w / 12);
	v = voltage[itime] + timeStep * (k1V / 12 + k2V / 3 - k3V / 12);
	t = time[itime] + timeStep / 3;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k4V = kV(v, t, h1, w1, m3, n4);
	k4n = kn(v, n1);
	k4m = km(v, m1);
	k4h = kh(v, h1);
	k4w = kw(v, w1);

	m1 = m[itime] + timeStep * (25 * k1m / 48 - 55 * k2m / 24 + 35 * k3m / 48 + 15 * k4m / 8);
	n1 = n[itime] + timeStep * (25 * k1n / 48 - 55 * k2n / 24 + 35 * k3n / 48 + 15 * k4n / 8);
	h1 = h[itime] + timeStep * (25 * k1h / 48 - 55 * k2h / 24 + 35 * k3h / 48 + 15 * k4h / 8);
	w1 = w[itime] + timeStep * (25 * k1w / 48 - 55 * k2w / 24 + 35 * k3w / 48 + 15 * k4w / 8);
	v = voltage[itime] + timeStep * (25 * k1V / 48 - 55 * k2V / 24 + 35 * k3V / 48 + 15 * k4V / 8);
	t = time[itime] + 5 * timeStep / 6;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k5V = kV(v, t, h1, w1, m3, n4);
	k5n = kn(v, n1);
	k5m = km(v, m1);
	k5h = kh(v, h1);
	k5w = kw(v, w1);
	
	m1 = m[itime] + timeStep * (3 * k1m / 20 - 11 * k2m / 24 - k3m / 8 + k4m / 2 + k5m / 10);
	n1 = n[itime] + timeStep * (3 * k1n / 20 - 11 * k2n / 24 - k3n / 8 + k4n / 2 + k5n / 10);
	h1 = h[itime] + timeStep * (3 * k1h / 20 - 11 * k2h / 24 - k3h / 8 + k4h / 2 + k5h / 10);
	w1 = w[itime] + timeStep * (3 * k1w / 20 - 11 * k2w / 24 - k3w / 8 + k4w / 2 + k5w / 10);
	v = voltage[itime] + timeStep * (3 * k1V / 20 - 11 * k2V / 24 - k3V / 8 + k4V / 2 + k5V / 10);
	t = time[itime] + timeStep / 6;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k6V = kV(v, t, h1, w1, m3, n4);
	k6n = kn(v, n1);
	k6m = km(v, m1);
	k6h = kh(v, h1);
	k6w = kw(v, w1);
	
	m1 = m[itime] + timeStep * (- 261 * k1m / 260 + 33 * k2m / 13 + 43 * k3m / 156 - 118 * k4m / 39 + 32 * k5m / 195 + 80 * k6m / 39);
	n1 = n[itime] + timeStep * (- 261 * k1n / 260 + 33 * k2n / 13 + 43 * k3n / 156 - 118 * k4n / 39 + 32 * k5n / 195 + 80 * k6n / 39);
	h1 = h[itime] + timeStep * (- 261 * k1h / 260 + 33 * k2h / 13 + 43 * k3h / 156 - 118 * k4h / 39 + 32 * k5h / 195 + 80 * k6h / 39);
	w1 = w[itime] + timeStep * (- 261 * k1w / 260 + 33 * k2w / 13 + 43 * k3w / 156 - 118 * k4w / 39 + 32 * k5w / 195 + 80 * k6w / 39);
	v = voltage[itime] + timeStep * (- 261 * k1V / 260 + 33 * k2V / 13 + 43 * k3V / 156 - 118 * k4V / 39 + 32 * k5V / 195 + 80 * k6V / 39);
	t = time[itime] + timeStep;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k7V = kV(v, t, h1, w1, m3, n4);
	k7n = kn(v, n1);
	k7m = km(v, m1);
	k7h = kh(v, h1);
	k7w = kw(v, w1);

	//	update all values for next time point

	voltage[itime + 1] = voltage[itime] + timeStep * (13 * (k1V + k7V) / 200 + 11 * (k3V + k4V) / 40 + 4 * (k5V + k6V) / 25);
	n[itime + 1] = n[itime] + timeStep * (13 * (k1n + k7n) / 200 + 11 * (k3n + k4n) / 40 + 4 * (k5n + k6n) / 25);
	m[itime + 1] = m[itime] + timeStep * (13 * (k1m + k7m) / 200 + 11 * (k3m + k4m) / 40 + 4 * (k5m + k6m) / 25);
	h[itime + 1] = h[itime] + timeStep * (13 * (k1h + k7h) / 200 + 11 * (k3h + k4h) / 40 + 4 * (k5h + k6h) / 25);
	w[itime + 1] = w[itime] + timeStep * (13 * (k1w + k7w) / 200 + 11 * (k3w + k4w) / 40 + 4 * (k5w + k6w) / 25);
	time[itime + 1] = time[itime] + timeStep;
	Gexc[itime + 1] = Ge(time[itime + 1]);
	Ginh[itime + 1] = Gi(time[itime + 1]);
	I[itime + 1] = I_total(time[itime + 1]);
	
	//	update internal time

	itime = itime + 1;
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
	I[itime + 1] = I_total(time[itime + 1]);
	//	update internal time

	itime = itime + 1;
}

void HHI_final::DP5_step()
{
	double m1, n1, h1, w1, m3, n4;
	double v, t;
	double k1V, k2V, k3V, k4V, k5V, k6V, k7V;
	double k1n, k2n, k3n, k4n, k5n, k6n, k7n;
	double k1m, k2m, k3m, k4m, k5m, k6m, k7m;
	double k1h, k2h, k3h, k4h, k5h, k6h, k7h;
	double k1w, k2w, k3w, k4w, k5w, k6w, k7w;

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

	m1 = m[itime] + timeStep * k1m / 5;
	n1 = n[itime] + timeStep * k1n / 5;
	h1 = h[itime] + timeStep * k1h / 5;
	w1 = w[itime] + timeStep * k1w / 5;
	v = voltage[itime] + timeStep * k1V / 5;
	t = time[itime] + timeStep / 5;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k2V = kV(v, t, h1, w1, m3, n4);
	k2n = kn(v, n1);
	k2m = km(v, m1);
	k2h = kh(v, h1);
	k2w = kw(v, w1);

	m1 = m[itime] + timeStep * (3 * k1m + 9 * k2m) / 40;
	n1 = n[itime] + timeStep * (3 * k1n + 9 * k2n) / 40;
	h1 = h[itime] + timeStep * (3 * k1h + 9 * k2h) / 40;
	w1 = w[itime] + timeStep * (3 * k1w + 9 * k2w) / 40;
	v = voltage[itime] + timeStep * (3 * k1V + 9 * k2V) / 40;
	t = time[itime] + 3 * timeStep / 10;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k3V = kV(v, t, h1, w1, m3, n4);
	k3n = kn(v, n1);
	k3m = km(v, m1);
	k3h = kh(v, h1);
	k3w = kw(v, w1);

	m1 = m[itime] + timeStep * (44 * k1m / 45 - 56 * k2m / 15 + 32 * k3m / 9);
	n1 = n[itime] + timeStep * (44 * k1n / 45 - 56 * k2n / 15 + 32 * k3n / 9);
	h1 = h[itime] + timeStep * (44 * k1h / 45 - 56 * k2h / 15 + 32 * k3h / 9);
	w1 = w[itime] + timeStep * (44 * k1w / 45 - 56 * k2w / 15 + 32 * k3w / 9);
	v = voltage[itime] + timeStep * (44 * k1V / 45 - 56 * k2V / 15 + 32 * k3V / 9);
	t = time[itime] + 4 * timeStep / 5;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k4V = kV(v, t, h1, w1, m3, n4);
	k4n = kn(v, n1);
	k4m = km(v, m1);
	k4h = kh(v, h1);
	k4w = kw(v, w1);

	m1 = m[itime] + timeStep * (19372 * k1m / 6561 - 25360 * k2m / 2187 + 64448 * k3m / 6561 - 212 * k4m / 729);
	n1 = n[itime] + timeStep * (19372 * k1n / 6561 - 25360 * k2n / 2187 + 64448 * k3n / 6561 - 212 * k4n / 729);
	h1 = h[itime] + timeStep * (19372 * k1h / 6561 - 25360 * k2h / 2187 + 64448 * k3h / 6561 - 212 * k4h / 729);
	w1 = w[itime] + timeStep * (19372 * k1w / 6561 - 25360 * k2w / 2187 + 64448 * k3w / 6561 - 212 * k4w / 729);
	v = voltage[itime] + timeStep * (19372 * k1V / 6561 - 25360 * k2V / 2187 + 64448 * k3V / 6561 - 212 * k4V / 729);
	t = time[itime] + 8 * timeStep / 9;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k5V = kV(v, t, h1, w1, m3, n4);
	k5n = kn(v, n1);
	k5m = km(v, m1);
	k5h = kh(v, h1);
	k5w = kw(v, w1);
	
	m1 = m[itime] + timeStep * (9017 * k1m / 3168 - 355 * k2m / 33 + 46732 * k3m / 5247 + 49 * k4m / 176 - 5103 * k5m / 18656);
	n1 = n[itime] + timeStep * (9017 * k1n / 3168 - 355 * k2n / 33 + 46732 * k3n / 5247 + 49 * k4n / 176 - 5103 * k5n / 18656);
	h1 = h[itime] + timeStep * (9017 * k1h / 3168 - 355 * k2h / 33 + 46732 * k3h / 5247 + 49 * k4h / 176 - 5103 * k5h / 18656);
	w1 = w[itime] + timeStep * (9017 * k1w / 3168 - 355 * k2w / 33 + 46732 * k3w / 5247 + 49 * k4w / 176 - 5103 * k5w / 18656);
	v = voltage[itime] + timeStep * (9017 * k1V / 3168 - 355 * k2V / 33 + 46732 * k3V / 5247 + 49 * k4V / 176 - 5103 * k5V / 18656);
	t = time[itime] + timeStep;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k6V = kV(v, t, h1, w1, m3, n4);
	k6n = kn(v, n1);
	k6m = km(v, m1);
	k6h = kh(v, h1);
	k6w = kw(v, w1);
	
	/*
	m1 = m[itime] + timeStep * (35 * k1m / 384 + 500 * k3m / 1113 + 125 * k4m / 192 - 2187 * k5m / 6784 + 11 * k6m / 84);
	n1 = n[itime] + timeStep * (35 * k1n / 384 + 500 * k3n / 1113 + 125 * k4n / 192 - 2187 * k5n / 6784 + 11 * k6n / 84);
	h1 = h[itime] + timeStep * (35 * k1h / 384 + 500 * k3h / 1113 + 125 * k4h / 192 - 2187 * k5h / 6784 + 11 * k6h / 84);
	w1 = w[itime] + timeStep * (35 * k1w / 384 + 500 * k3w / 1113 + 125 * k4w / 192 - 2187 * k5w / 6784 + 11 * k6w / 84);
	v = voltage[itime] + timeStep * (35 * k1V / 384 + 500 * k3V / 1113 + 125 * k4V / 192 - 2187 * k5V / 6784 + 11 * k6V / 84);
	t = time[itime] + timeStep;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k7V = kV(v, t, h1, w1, m3, n4);
	k7n = kn(v, n1);
	k7m = km(v, m1);
	k7h = kh(v, h1);
	k7w = kw(v, w1);
	*/

	//	update all values for next time point

	voltage[itime + 1] = voltage[itime] + timeStep * (35 * k1V / 384 + 500 * k3V / 1113 + 125 * k4V / 192 - 2187 * k5V / 6784 + 11 * k6V / 84);
	n[itime + 1] = n[itime] + timeStep * (35 * k1n / 384 + 500 * k3n / 1113 + 125 * k4n / 192 - 2187 * k5n / 6784 + 11 * k6n / 84);
	m[itime + 1] = m[itime] + timeStep * (35 * k1m / 384 + 500 * k3m / 1113 + 125 * k4m / 192 - 2187 * k5m / 6784 + 11 * k6m / 84);
	h[itime + 1] = h[itime] + timeStep * (35 * k1h / 384 + 500 * k3h / 1113 + 125 * k4h / 192 - 2187 * k5h / 6784 + 11 * k6h / 84);
	w[itime + 1] = w[itime] + timeStep * (35 * k1w / 384 + 500 * k3w / 1113 + 125 * k4w / 192 - 2187 * k5w / 6784 + 11 * k6w / 84);
	time[itime + 1] = time[itime] + timeStep;
	Gexc[itime + 1] = Ge(time[itime + 1]);
	Ginh[itime + 1] = Gi(time[itime + 1]);
	I[itime + 1] = I_total(time[itime + 1]);

	// calculate errors
	/*
	double vBetter, vError;

	vBetter = voltage[itime] + timeStep * (5179 * k1V / 57600 + 7571 * k3V / 16695 + 393 * k4V / 640 - 92097 * k5V / 339200 + 187 * k6V / 2100 + k7V / 40);

	vError = fabs(voltage[itime + 1] - vBetter);

	std::cout << "Error is " << vError << std::endl;
	*/
	//	update internal time

	itime = itime + 1;
}

void HHI_final::DP8_step()
{

	double m1, n1, h1, w1, m3, n4;
	double v, t;
	double k1V, k2V, k3V, k4V, k5V, k6V, k7V, k8V, k9V, k10V, k11V, k12V, k13V;
	double k1n, k2n, k3n, k4n, k5n, k6n, k7n, k8n, k9n, k10n, k11n, k12n, k13n;
	double k1m, k2m, k3m, k4m, k5m, k6m, k7m, k8m, k9m, k10m, k11m, k12m, k13m;
	double k1h, k2h, k3h, k4h, k5h, k6h, k7h, k8h, k9h, k10h, k11h, k12h, k13h;
	double k1w, k2w, k3w, k4w, k5w, k6w, k7w, k8w, k9w, k10w, k11w, k12w, k13w;

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

	m1 = m[itime] + timeStep * k1m / 18;
	n1 = n[itime] + timeStep * k1n / 18;
	h1 = h[itime] + timeStep * k1h / 18;
	w1 = w[itime] + timeStep * k1w / 18;
	v = voltage[itime] + timeStep * k1V / 18;
	t = time[itime] + timeStep / 18;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k2V = kV(v, t, h1, w1, m3, n4);
	k2n = kn(v, n1);
	k2m = km(v, m1);
	k2h = kh(v, h1);
	k2w = kw(v, w1);

	m1 = m[itime] + timeStep * (k1m + 3 * k2m) / 48;
	n1 = n[itime] + timeStep * (k1n + 3 * k2n) / 48;
	h1 = h[itime] + timeStep * (k1h + 3 * k2h) / 48;
	w1 = w[itime] + timeStep * (k1w + 3 * k2w) / 48;
	v = voltage[itime] + timeStep * (k1V + 3 * k2V) / 48;
	t = time[itime] + timeStep / 12;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k3V = kV(v, t, h1, w1, m3, n4);
	k3n = kn(v, n1);
	k3m = km(v, m1);
	k3h = kh(v, h1);
	k3w = kw(v, w1);

	m1 = m[itime] + timeStep * (k1m + 3 * k3m) / 32;
	n1 = n[itime] + timeStep * (k1n + 3 * k3n) / 32;
	h1 = h[itime] + timeStep * (k1h + 3 * k3h) / 32;
	w1 = w[itime] + timeStep * (k1w + 3 * k3w) / 32;
	v = voltage[itime] + timeStep * (k1V + 3 * k3V) / 32;
	t = time[itime] + timeStep / 8;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k4V = kV(v, t, h1, w1, m3, n4);
	k4n = kn(v, n1);
	k4m = km(v, m1);
	k4h = kh(v, h1);
	k4w = kw(v, w1);

	m1 = m[itime] + timeStep * (20 * k1m - 75 * k3m + 75 * k4m) / 64;
	n1 = n[itime] + timeStep * (20 * k1n - 75 * k3n + 75 * k4n) / 64;
	h1 = h[itime] + timeStep * (20 * k1h - 75 * k3h + 75 * k4h) / 64;
	w1 = w[itime] + timeStep * (20 * k1w - 75 * k3w + 75 * k4w) / 64;
	v = voltage[itime] + timeStep * (20 * k1V - 75 * k3V + 75 * k4V) / 64;
	t = time[itime] + 5 * timeStep / 16;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k5V = kV(v, t, h1, w1, m3, n4);
	k5n = kn(v, n1);
	k5m = km(v, m1);
	k5h = kh(v, h1);
	k5w = kw(v, w1);
	
	m1 = m[itime] + timeStep * (3 * k1m + 15 * k4m + 12 * k5m) / 80;
	n1 = n[itime] + timeStep * (3 * k1n + 15 * k4n + 12 * k5n) / 80;
	h1 = h[itime] + timeStep * (3 * k1h + 15 * k4h + 12 * k5h) / 80;
	w1 = w[itime] + timeStep * (3 * k1w + 15 * k4w + 12 * k5w) / 80;
	v = voltage[itime] + timeStep * (3 * k1V + 15 * k4V + 12 * k5V) / 80;
	t = time[itime] + 3 * timeStep / 8;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k6V = kV(v, t, h1, w1, m3, n4);
	k6n = kn(v, n1);
	k6m = km(v, m1);
	k6h = kh(v, h1);
	k6w = kw(v, w1);
	
	m1 = m[itime] + timeStep * (29443841 * k1m / 614563906 + 77736538 * k4m / 692538347 - 28693883 * k5m / 1125000000 + 23124283 * k6m / 1800000000);
	n1 = n[itime] + timeStep * (29443841 * k1n / 614563906 + 77736538 * k4n / 692538347 - 28693883 * k5n / 1125000000 + 23124283 * k6n / 1800000000);
	h1 = h[itime] + timeStep * (29443841 * k1h / 614563906 + 77736538 * k4h / 692538347 - 28693883 * k5h / 1125000000 + 23124283 * k6h / 1800000000);
	w1 = w[itime] + timeStep * (29443841 * k1w / 614563906 + 77736538 * k4w / 692538347 - 28693883 * k5w / 1125000000 + 23124283 * k6w / 1800000000);
	v = voltage[itime] + timeStep * (29443841 * k1V / 614563906 + 77736538 * k4V / 692538347 - 28693883 * k5V / 1125000000 + 23124283 * k6V / 1800000000);
	t = time[itime] + 59 * timeStep / 400;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k7V = kV(v, t, h1, w1, m3, n4);
	k7n = kn(v, n1);
	k7m = km(v, m1);
	k7h = kh(v, h1);
	k7w = kw(v, w1);
	
	m1 = m[itime] + timeStep * (16016141 * k1m / 946692911 + 61564180 * k4m / 158732637 + 22789713 * k5m / 633445777 + 545815736 * k6m / 2771057229 - 180193667 * k7m / 1043307555);
	
	n1 = n[itime] + timeStep * (16016141 * k1n / 946692911 + 61564180 * k4n / 158732637 + 22789713 * k5n / 633445777 + 545815736 * k6n / 2771057229 - 180193667 * k7n / 1043307555);
	
	h1 = h[itime] + timeStep * (16016141 * k1h / 946692911 + 61564180 * k4h / 158732637 + 22789713 * k5h / 633445777 + 545815736 * k6h / 2771057229 - 180193667 * k7h / 1043307555);
	
	w1 = w[itime] + timeStep * (16016141 * k1w / 946692911 + 61564180 * k4w / 158732637 + 22789713 * k5w / 633445777 + 545815736 * k6w / 2771057229 - 180193667 * k7w / 1043307555);
	
	v = voltage[itime] + timeStep * (16016141 * k1V / 946692911 + 61564180 * k4V / 158732637 + 22789713 * k5V / 633445777 + 545815736 * k6V / 2771057229 - 180193667 * k7V / 1043307555);
	
	t = time[itime] + 93 * timeStep / 200;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k8V = kV(v, t, h1, w1, m3, n4);
	k8n = kn(v, n1);
	k8m = km(v, m1);
	k8h = kh(v, h1);
	k8w = kw(v, w1);
	
	m1 = m[itime] + timeStep * (39632708 * k1m / 573591083 - 433636366 * k4m / 683701615 - 421739975 * k5m / 2616292301 + 100302831 * k6m / 723423059 + 790204164 * k7m / 839813087 + 800635310 * k8m / 3783071287);
	
	n1 = n[itime] + timeStep * (39632708 * k1n / 573591083 - 433636366 * k4n / 683701615 - 421739975 * k5n / 2616292301 + 100302831 * k6n / 723423059 + 790204164 * k7n / 839813087 + 800635310 * k8n / 3783071287);
	
	h1 = h[itime] + timeStep * (39632708 * k1h / 573591083 - 433636366 * k4h / 683701615 - 421739975 * k5h / 2616292301 + 100302831 * k6h / 723423059 + 790204164 * k7h / 839813087 + 800635310 * k8h / 3783071287);
	
	w1 = w[itime] + timeStep * (39632708 * k1w / 573591083 - 433636366 * k4w / 683701615 - 421739975 * k5w / 2616292301 + 100302831 * k6w / 723423059 + 790204164 * k7w / 839813087 + 800635310 * k8w / 3783071287);
	
	v = voltage[itime] + timeStep * (39632708 * k1V / 573591083 - 433636366 * k4V / 683701615 - 421739975 * k5V / 2616292301 + 100302831 * k6V / 723423059 + 790204164 * k7V / 839813087 + 800635310 * k8V / 3783071287);
	
	t = time[itime] + 5490023248 * timeStep / 9719169821;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k9V = kV(v, t, h1, w1, m3, n4);
	k9n = kn(v, n1);
	k9m = km(v, m1);
	k9h = kh(v, h1);
	k9w = kw(v, w1);
	
	m1 = m[itime] + timeStep * (246121993 * k1m / 1340847787 - 37695042795 * k4m / 15268766246 - 309121744 * k5m / 1061227803 - 12992083 * k6m / 490766935 + 6005943493 * k7m / 2108947869 + 393006217 * k8m / 1396673457 + 123872331 * k9m / 1001029789);
	
	n1 = n[itime] + timeStep * (246121993 * k1n / 1340847787 - 37695042795 * k4n / 15268766246 - 309121744 * k5n / 1061227803 - 12992083 * k6n / 490766935 + 6005943493 * k7n / 2108947869 + 393006217 * k8n / 1396673457 + 123872331 * k9n / 1001029789);
	
	h1 = h[itime] + timeStep * (246121993 * k1h / 1340847787 - 37695042795 * k4h / 15268766246 - 309121744 * k5h / 1061227803 - 12992083 * k6h / 490766935 + 6005943493 * k7h / 2108947869 + 393006217 * k8h / 1396673457 + 123872331 * k9h / 1001029789);
	
	w1 = w[itime] + timeStep * (246121993 * k1w / 1340847787 - 37695042795 * k4w / 15268766246 - 309121744 * k5w / 1061227803 - 12992083 * k6w / 490766935 + 6005943493 * k7w / 2108947869 + 393006217 * k8w / 1396673457 + 123872331 * k9w / 1001029789);
	
	v = voltage[itime] + timeStep * (246121993 * k1V / 1340847787 - 37695042795 * k4V / 15268766246 - 309121744 * k5V / 1061227803 - 12992083 * k6V / 490766935 + 6005943493 * k7V / 2108947869 + 393006217 * k8V / 1396673457 + 123872331 * k9V / 1001029789);
	
	t = time[itime] + 13 * timeStep / 20;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k10V = kV(v, t, h1, w1, m3, n4);
	k10n = kn(v, n1);
	k10m = km(v, m1);
	k10h = kh(v, h1);
	k10w = kw(v, w1);
	
	m1 = m[itime] + timeStep * (- 1028468189 * k1m / 846180014 + 8478235783 * k4m / 508512852 + 1311729495 * k5m / 1432422823 - 10304129995 * k6m / 1701304382 - 48777925059 * k7m / 3047939560 + 15336726248 * k8m / 1032824649 - 45442868181 * k9m / 3398467696 + 3065993473 * k10m / 597172653);
	
	n1 = n[itime] + timeStep * (- 1028468189 * k1n / 846180014 + 8478235783 * k4n / 508512852 + 1311729495 * k5n / 1432422823 - 10304129995 * k6n / 1701304382 - 48777925059 * k7n / 3047939560 + 15336726248 * k8n / 1032824649 - 45442868181 * k9n / 3398467696 + 3065993473 * k10n / 597172653);
	
	h1 = h[itime] + timeStep * (- 1028468189 * k1h / 846180014 + 8478235783 * k4h / 508512852 + 1311729495 * k5h / 1432422823 - 10304129995 * k6h / 1701304382 - 48777925059 * k7h / 3047939560 + 15336726248 * k8h / 1032824649 - 45442868181 * k9h / 3398467696 + 3065993473 * k10h / 597172653);
	
	w1 = w[itime] + timeStep * (- 1028468189 * k1w / 846180014 + 8478235783 * k4w / 508512852 + 1311729495 * k5w / 1432422823 - 10304129995 * k6w / 1701304382 - 48777925059 * k7w / 3047939560 + 15336726248 * k8w / 1032824649 - 45442868181 * k9w / 3398467696 + 3065993473 * k10w / 597172653);
	
	v = voltage[itime] + timeStep * (- 1028468189 * k1V / 846180014 + 8478235783 * k4V / 508512852 + 1311729495 * k5V / 1432422823 - 10304129995 * k6V / 1701304382 - 48777925059 * k7V / 3047939560 + 15336726248 * k8V / 1032824649 - 45442868181 * k9V / 3398467696 + 3065993473 * k10V / 597172653);
	
	t = time[itime] + 1201146811 * timeStep / 12990119798;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k11V = kV(v, t, h1, w1, m3, n4);
	k11n = kn(v, n1);
	k11m = km(v, m1);
	k11h = kh(v, h1);
	k11w = kw(v, w1);
	
	m1 = m[itime] + timeStep * (185892177 * k1m / 718116043 - 3185094517 * k4m / 667107341 - 477755414 * k5m / 1098053517 - 703635378 * k6m / 230739211 + 5731566787 * k7m / 1027545527 + 5232866602 * k8m / 850066563 - 4093664535 * k9m / 808688257 + 3962137247 * k10m / 1805957418 + 65686358 * k11m / 487910083);
	
	n1 = n[itime] + timeStep * (185892177 * k1n / 718116043 - 3185094517 * k4n / 667107341 - 477755414 * k5n / 1098053517 - 703635378 * k6n / 230739211 + 5731566787 * k7n / 1027545527 + 5232866602 * k8n / 850066563 - 4093664535 * k9n / 808688257 + 3962137247 * k10n / 1805957418 + 65686358 * k11n / 487910083);
	
	h1 = h[itime] + timeStep * (185892177 * k1h / 718116043 - 3185094517 * k4h / 667107341 - 477755414 * k5h / 1098053517 - 703635378 * k6h / 230739211 + 5731566787 * k7h / 1027545527 + 5232866602 * k8h / 850066563 - 4093664535 * k9h / 808688257 + 3962137247 * k10h / 1805957418 + 65686358 * k11h / 487910083);
	
	w1 = w[itime] + timeStep * (185892177 * k1w / 718116043 - 3185094517 * k4w / 667107341 - 477755414 * k5w / 1098053517 - 703635378 * k6w / 230739211 + 5731566787 * k7w / 1027545527 + 5232866602 * k8w / 850066563 - 4093664535 * k9w / 808688257 + 3962137247 * k10w / 1805957418 + 65686358 * k11w / 487910083);
	
	v = voltage[itime] + timeStep * (185892177 * k1V / 718116043 - 3185094517 * k4V / 667107341 - 477755414 * k5V / 1098053517 - 703635378 * k6V / 230739211 + 5731566787 * k7V / 1027545527 + 5232866602 * k8V / 850066563 - 4093664535 * k9V / 808688257 + 3962137247 * k10V / 1805957418 + 65686358 * k11V / 487910083);
	
	t = time[itime] + timeStep;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k12V = kV(v, t, h1, w1, m3, n4);
	k12n = kn(v, n1);
	k12m = km(v, m1);
	k12h = kh(v, h1);
	k12w = kw(v, w1);

	/*
	m1 = m[itime] + timeStep * (403863854 * k1m / 491063109 - 5068492393 * k4m / 434740067 - 411421997 * k5m / 543043805 + 652783627 * k6m / 914296604 + 11173962825 * k7m / 925320556 - 13158990841 * k8m / 6184727034 + 3936647629 * k9m / 1978049680 - 160528059 * k10m / 685178525 + 248638103 * k11m / 1413531060;
	
	n1 = n[itime] + timeStep * (403863854 * k1n / 491063109 - 5068492393 * k4n / 434740067 - 411421997 * k5n / 543043805 + 652783627 * k6n / 914296604 + 11173962825 * k7n / 925320556 - 13158990841 * k8n / 6184727034 + 3936647629 * k9n / 1978049680 - 160528059 * k10n / 685178525 + 248638103 * k11n / 1413531060;
	
	h1 = h[itime] + timeStep * (403863854 * k1h / 491063109 - 5068492393 * k4h / 434740067 - 411421997 * k5h / 543043805 + 652783627 * k6h / 914296604 + 11173962825 * k7h / 925320556 - 13158990841 * k8h / 6184727034 + 3936647629 * k9h / 1978049680 - 160528059 * k10h / 685178525 + 248638103 * k11h / 1413531060;
	
	w1 = w[itime] + timeStep * (403863854 * k1w / 491063109 - 5068492393 * k4w / 434740067 - 411421997 * k5w / 543043805 + 652783627 * k6w / 914296604 + 11173962825 * k7w / 925320556 - 13158990841 * k8w / 6184727034 + 3936647629 * k9w / 1978049680 - 160528059 * k10w / 685178525 + 248638103 * k11w / 1413531060;

	v = voltage[itime] + timeStep * (403863854 * k1V / 491063109 - 5068492393 * k4V / 434740067 - 411421997 * k5V / 543043805 + 652783627 * k6V / 914296604 + 11173962825 * k7V / 925320556 - 13158990841 * k8V / 6184727034 + 3936647629 * k9V / 1978049680 - 160528059 * k10V / 685178525 + 248638103 * k11V / 1413531060;
	
	t = time[itime] + timeStep;
	n4 = n1 * n1 * n1 * n1;
	m3 = m1 * m1 * m1;

	k13V = kV(v, t, h1, w1, m3, n4);
	k13n = kn(v, n1);
	k13m = km(v, m1);
	k13h = kh(v, h1);
	k13w = kw(v, w1);
	*/
	//	update all values for next time point

	voltage[itime + 1] = voltage[itime] + timeStep * (13451932 * k1V / 455176623 - 808719846 * k6V / 976000145 + 1757004468 * k7V / 5645159321 + 656045339 * k8V / 265891186 - 3867574721 * k9V / 1518517206 + 465885868 * k10V / 322736535 + 53011238 * k11V / 667516719 + 2 * k12V / 45);
	n[itime + 1] = n[itime] + timeStep * (13451932 * k1n / 455176623 - 808719846 * k6n / 976000145 + 1757004468 * k7n / 5645159321 + 656045339 * k8n / 265891186 - 3867574721 * k9n / 1518517206 + 465885868 * k10n / 322736535 + 53011238 * k11n / 667516719 + 2 * k12n / 45);
	m[itime + 1] = m[itime] + timeStep * (13451932 * k1m / 455176623 - 808719846 * k6m / 976000145 + 1757004468 * k7m / 5645159321 + 656045339 * k8m / 265891186 - 3867574721 * k9m / 1518517206 + 465885868 * k10m / 322736535 + 53011238 * k11m / 667516719 + 2 * k12m / 45);
	h[itime + 1] = h[itime] + timeStep * (13451932 * k1h / 455176623 - 808719846 * k6h / 976000145 + 1757004468 * k7h / 5645159321 + 656045339 * k8h / 265891186 - 3867574721 * k9h / 1518517206 + 465885868 * k10h / 322736535 + 53011238 * k11h / 667516719 + 2 * k12h / 45);
	w[itime + 1] = w[itime] + timeStep * (13451932 * k1w / 455176623 - 808719846 * k6w / 976000145 + 1757004468 * k7w / 5645159321 + 656045339 * k8w / 265891186 - 3867574721 * k9w / 1518517206 + 465885868 * k10w / 322736535 + 53011238 * k11w / 667516719 + 2 * k12w / 45);
	time[itime + 1] = time[itime] + timeStep;
	Gexc[itime + 1] = Ge(time[itime + 1]);
	Ginh[itime + 1] = Gi(time[itime + 1]);
	I[itime + 1] = I_total(time[itime + 1]);

	// calculate errors
	/*
	double vBetter, vError;


	vBetter = voltage[itime] + timeStep * (14005451 * k1V / 335480064 - 59238493 * k6V / 1068277825 + 181606767 * k7V / 758867731 + 561292985 * k8V / 797845732 - 1041891430 * k9V / 1371343529 + 760417239 * k10V / 1151165299 + 118820643 * k11V / 751138087 - 528747749 * k12V / 2220607170 + k13V / 4);

	vError = fabs(voltage[itime + 1] - vBetter);

	std::cout << "Error is " << vError << std::endl;
	*/
	//	update internal time

	itime = itime + 1;
}

double HHI_final::bisection(DDfunction f, double x, double xmax, double eps, int Nmax)
{
	if ( (fabs(f(x)) <= eps) )
	{
		//std::cout << "Precision is fulfilled for the value on the previous step" << std::endl;

		return x;
	}
	
	double a; // left point of the interval
	double b; // right point of the interval
	
	// first let's find values for x where f(x) has different signs

	double step = 0.1; // step
	int i = 1;

	double f_a; // value of f at x=a
	double f_b; // value of f at x=b
	
	do
	{
		a = x - step * i;
		b = x + step * i;

		f_a = f(a);
		f_b = f(b);
		
		if (f_a * f_b < 0)
		{	
		//	std::cout << "a = " << a << " f(a) = " << f(a) << std::endl;
		//	std::cout << "b = " << b << " f(b) = " << f(b) << std::endl;
		//	std::cout << "x = " << x << " f(x) = " << f(x) << std::endl;
		//	std::cout << "number of iterations used to find interval = " << i << std::endl;
			break;
		}

		i ++;
	} while (b <= xmax);

	if (b > xmax)
		std::cerr << "bisection method did not find interval where function takes opposite signs" << std::endl;
	
	
	
	int N = 0;// number of iterations performed
	double c = (a + b) / 2; // middle of the interval

	double f_c; // value of f at x=c


	while (N < Nmax)
	{
		N++;

		f_c = f(c);

		if ( (fabs(f_c) <= eps) || (b - c) <= eps)
		{
			//std::cout << "Number of iterations performed for eps = " << eps << " is " << N << std::endl;
			return c;
		}
		

		//std::cout << "Before updating c" << std::endl;
		//std::cout << "a = " << a << " f(a) = " << f(a) << std::endl;
		//std::cout << "b = " << b << " f(b) = " << f(b) << std::endl;
		//std::cout << "c = " << c << " f(c) = " << f(c) << std::endl;

		if ( (f_c * f_b < 0) || (f_a * f_c > 0) )
		{
			a = c;
			f_a = f_c;
		}
		else if ( (f_a * f_c < 0) || (f_b * f_c > 0) )
		{
			b = c;
			f_b = f_c;
		}
		
		c = (a + b) / 2;
		
		//std::cout << "After updating c" << std::endl;
		//std::cout << "a = " << a << " f(a) = " << f(a) << std::endl;
		//std::cout << "b = " << b << " f(b) = " << f(b) << std::endl;
		//std::cout << "c = " << c << " f(c) = " << f(c) << std::endl;
	}
	//std::cout << "All " << N << " iterations were made in bisection method! fabs(f(c)) - eps = " << fabs(f_c) - eps << std::endl;
	
	return c;
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
		targets[i]->raiseI(targetsG[i]);
	}
}

double HHI_final::I_total(double t)
{
	if (injected_current)
		return I_injected(t) + mu;
	else
		return I_default(t) + mu;

}

double HHI_final::f_iTrapezoid(double fv, double fn, double fm, double fh, double fw,  double v)
{
	double m_ = m_iTrapezoid(fm, v);
	double n_ = n_iTrapezoid(fn, v);

	double m3 = m_ * m_ * m_;
	double n4 = n_ * n_ * n_ * n_;

	return v - voltage[itime] - 0.5 * timeStep * ( fv + (-gL * (v - El) - 
				gNa * h_iTrapezoid(fh, v) * m3 * (v - Ena) - 
				gKdr * n4 * (v - Ek) - 
				gKHT * w_iTrapezoid(fw, v) * (v - Ek) - Gexc[itime + 1] * v - Ginh[itime + 1] * (v - Ei) + I[itime + 1]) / cm) ;
}


double HHI_final::m_iTrapezoid(double fm, double v)
{
	return (m[itime] + 0.5 * timeStep * (fm + am(v))) / (1 + 0.5 * timeStep * (am(v) + bm(v)));
}

double HHI_final::n_iTrapezoid(double fn, double v)
{
	return (n[itime] + 0.5 * timeStep * (fn + an(v))) / (1 + 0.5 * timeStep * (an(v) + bn(v)));
}

double HHI_final::h_iTrapezoid(double fh, double v)
{
	return (h[itime] + 0.5 * timeStep * (fh + ah(v))) / (1 + 0.5 * timeStep * (ah(v) + bh(v)));
}

double HHI_final::w_iTrapezoid(double fw, double v)
{
	return (w[itime] + 0.5 * timeStep * (fw + wInf(v) / tauW(v))) / (1 + 0.5 * timeStep / tauW(v));
}

double HHI_final::f_iEuler(double v)
{
	double m_ = m_iEuler(v);
	double n_ = n_iEuler(v);

	double m3 = m_ * m_ * m_;
	double n4 = n_ * n_ * n_ * n_;

	return v - voltage[itime] - timeStep * ( (-gL * (v - El) - 
				gNa * h_iEuler(v) * m3 * (v - Ena) - 
				gKdr * n4 * (v - Ek) - 
				gKHT * w_iEuler(v) * (v - Ek) - Gexc[itime + 1] * v - Ginh[itime + 1] * (v - Ei) + I[itime + 1]) / cm) ;
}


double HHI_final::m_iEuler(double v)
{
	return (m[itime] + timeStep * am(v)) / (1 + timeStep * (am(v) + bm(v)));
}

double HHI_final::n_iEuler(double v)
{
	return (n[itime] + timeStep * an(v)) / (1 + timeStep * (an(v) + bn(v)));
}

double HHI_final::h_iEuler(double v)
{
	return (h[itime] + timeStep * ah(v)) / (1 + timeStep * (ah(v) + bh(v)));
}

double HHI_final::w_iEuler(double v)
{
	return (w[itime] + timeStep * wInf(v) / tauW(v)) / (1 + timeStep / tauW(v));
}

double HHI_final::f_bdf4(double v)
{
	double m_ = m_bdf4(v);
	double n_ = n_bdf4(v);

	double m3 = m_ * m_ * m_;
	double n4 = n_ * n_ * n_ * n_;
	
	return v - (48 * voltage[itime] - 36 * voltage[itime - 1] + 16 * voltage[itime - 2] - 3 * voltage[itime - 3] + 12 * timeStep * ( (-gL * (v - El) - 
				gNa * h_bdf4(v) * m3 * (v - Ena) - gKdr * n4 * (v - Ek) - 
				gKHT * w_bdf4(v) * (v - Ek) - Gexc[itime + 1] * v - Ginh[itime + 1] * (v - Ei) + I[itime + 1]) / cm) ) / 25;
}

double HHI_final::m_bdf4(double v)
{
	return (48 * m[itime] - 36 * m[itime - 1] + 16 * m[itime - 2] - 3 * m[itime - 3] + 12 * timeStep * am(v)) / (25 + 12 * timeStep * (am(v) + bm(v)));
}

double HHI_final::n_bdf4(double v)
{
	return (48 * n[itime] - 36 * n[itime - 1] + 16 * n[itime - 2] - 3 * n[itime - 3] + 12 * timeStep * an(v)) / (25 + 12 * timeStep * (an(v) + bn(v)));
}

double HHI_final::h_bdf4(double v)
{
	return (48 * h[itime] - 36 * h[itime - 1] + 16 * h[itime - 2] - 3 * h[itime - 3] + 12 * timeStep * ah(v)) / (25 + 12 * timeStep * (ah(v) + bh(v)));
}


double HHI_final::w_bdf4(double v)
{
	return (48 * w[itime] - 36 * w[itime - 1] + 16 * w[itime - 2] - 3 * w[itime - 3] + 12 * timeStep * wInf(v) / tauW(v)) / (25 + 12 * timeStep / tauW(v));
}

double HHI_final::f_bdf6(double v)
{
	double m_ = m_bdf6(v);
	double n_ = n_bdf6(v);

	double m3 = m_ * m_ * m_;
	double n4 = n_ * n_ * n_ * n_;
	
	return v - (360 * voltage[itime] - 450 * voltage[itime - 1] + 400 * voltage[itime - 2] - 225 * voltage[itime - 3] + 72 * voltage[itime - 4] - 
			10 * voltage[itime - 5] + 60 * timeStep * ( (-gL * (v - El) - gNa * h_bdf4(v) * m3 * (v - Ena) - gKdr * n4 * (v - Ek) - 
				gKHT * w_bdf4(v) * (v - Ek) - Gexc[itime + 1] * v - Ginh[itime + 1] * (v - Ei) + I[itime + 1]) / cm) ) / 147;
}

double HHI_final::m_bdf6(double v)
{
	return (360 * m[itime] - 450 * m[itime - 1] + 400 * m[itime - 2] - 225 * m[itime - 3] + 72 * m[itime - 4] - 10 * m[itime - 5] + 60 * timeStep * am(v)) /
			(147 + 60 * timeStep * (am(v) + bm(v)));
}

double HHI_final::n_bdf6(double v)
{
	return (360 * n[itime] - 450 * n[itime - 1] + 400 * n[itime - 2] - 225 * n[itime - 3] + 72 * n[itime - 4] - 10 * n[itime - 5] + 60 * timeStep * an(v)) /
			(147 + 60 * timeStep * (an(v) + bn(v)));
}

double HHI_final::h_bdf6(double v)
{
	return (360 * h[itime] - 450 * h[itime - 1] + 400 * h[itime - 2] - 225 * h[itime - 3] + 72 * h[itime - 4] - 10 * h[itime - 5] + 60 * timeStep * ah(v)) /
			(147 + 60 * timeStep * (ah(v) + bh(v)));
}
double HHI_final::w_bdf6(double v)
{
	return (360 * w[itime] - 450 * w[itime - 1] + 400 * w[itime - 2] - 225 * w[itime - 3] + 72 * w[itime - 4] - 10 * w[itime - 5] +
	60 * timeStep * wInf(v) / tauW(v)) / (147 + 60 * timeStep / tauW(v));
}

double HHI_final::Ge(double t){return Gexc[itime] * exp(-(t - time[itime]) / tExc);}
double HHI_final::Gi(double t){return Ginh[itime] * exp(-(t - time[itime]) / tInh);}
double HHI_final::kV(double v, double t, double h, double w, double m3, double n4){
	return (-gL * (v - El) - gNa * h * m3 * (v - Ena) - gKdr * n4 * (v - Ek)
		- gKHT * w * (v - Ek) - Ge(t) * v - Gi(t) * (v - Ei) + I_total(t)) / cm;}
double HHI_final::kn(double v, double n){return HHI_final::an(v)*(1 - n) - HHI_final::bn(v)*n;}
double HHI_final::km(double v, double m){return HHI_final::am(v)*(1 - m) - HHI_final::bm(v)*m;}
double HHI_final::kh(double v, double h){return HHI_final::ah(v)*(1 - h) - HHI_final::bh(v)*h;}
double HHI_final::kw(double v, double w){return (HHI_final::wInf(v) - w) / HHI_final::tauW(v);}

