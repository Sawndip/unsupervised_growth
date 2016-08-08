#include "HHI_final_pool.h"

void HHI_final_pool::set_coordinates(double xx, double yy)
{
	x = xx;
	y = yy;
}

void HHI_final_pool::get_coordinates(double& xx, double& yy)
{
	xx = x;
	yy = y;
}

void HHI_final_pool::reset()
{
     if ((std::isnan(voltage.back()))||(std::isinf(voltage.back())))
        printf("Dynamics of HHI neuron blew up\n");

	noise_exc -= itime;
	noise_inh -= itime;

	itime = 0;

	time[0] = 0;
	voltage[0] = voltage.back();
	n[0] = n.back();
	m[0] = m.back();
	h[0] = h.back();
	w[0] = w.back();
	Gexc[0] = Gexc.back();
	Ginh[0] = Ginh.back();
	flag[0] = flag.back();
	I[0] = Iext(time[0]);

	//this->initialize_noise(noise_exc);
	//this->initialize_noise(noise_inh);
}
