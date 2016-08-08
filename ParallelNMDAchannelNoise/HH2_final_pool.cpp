#include "HH2_final_pool.h"
#include <cmath>

void HH2_final_pool::set_coordinates(double xx, double yy)
{
	x = xx;
	y = yy;
}

void HH2_final_pool::get_coordinates(double& xx, double& yy)
{
	xx = x;
	yy = y;
}

void HH2_final_pool::reset()
{
    if ((std::isnan(Vs.back()))||(std::isinf(Vs.back()))||
        (std::isnan(Vd.back()))||(std::isinf(Vd.back())))
        printf("Dynamics of RA neuron blew up\n");

	noise_exc_soma -= itime;
	noise_inh_soma -= itime;
	noise_exc_dend -= itime;
	noise_inh_dend -= itime;


	itime = 0;



	time[0] = 0;
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

	flag_soma[0] = flag_soma.back();
    flag_dend[0] = flag_dend.back();

	Is[0] = IsExt(time[0]);
	Id[0] = IdExt(time[0]);

	//this->initialize_noise(noise_exc_soma);
	//this->initialize_noise(noise_inh_soma);
	//this->initialize_noise(noise_exc_dend);
	//this->initialize_noise(noise_inh_dend);
}

void HH2_final_pool::set_Ei(double E){Ei = E;}

void HH2_final_pool::get_state(iter_v_doubles& it_begin, iter_v_doubles& it_end)
{
    std::vector<double> state;

    state.push_back(time[itime]);
    state.push_back(Vs[itime]);

    it_begin = state.begin();
    it_end = state.end();
}
