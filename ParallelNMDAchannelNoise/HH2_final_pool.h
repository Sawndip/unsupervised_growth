#pragma once
#ifndef HH2_FINAL_POOL
#define HH2_FINAL_POOL

#include "HH2_final.h"
#include <vector>

typedef std::vector<double>::iterator iter_v_doubles;

class HH2_final_pool : public HH2_final
{
public:
	HH2_final_pool() : HH2_final(){};
	
	void set_coordinates(double xx, double yy); // set coordinates of a neuron
	void get_coordinates(double& xx, double& yy); // get coordinates of a neuron
	void set_Ei(double E); // set reverse GABA potential
	void reset(); // last values in all vectors become first values
    void get_state(iter_v_doubles& it_begin, iter_v_doubles& it_end); // return state of the neuron

protected:
	double x; // x - coordinate
	double y; // y - coordinate
	
};


#endif
