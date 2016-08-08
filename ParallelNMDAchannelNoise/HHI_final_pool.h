#pragma once
#ifndef HHI_FINAL_POOL
#define HHI_FINAL_POOL

#include "HHI_final.h"

class HHI_final_pool : public HHI_final
{
public:
	HHI_final_pool() : HHI_final(){};
	
	void set_coordinates(double xx, double yy); // set coordinates of a neuron
	void get_coordinates(double& xx, double& yy); // get coordinates of a neuron
	void reset(); // reset dynamics (last values in all arrays become first values)
protected:
	double x; // x - coordinate of a neuron
	double y; // y - coordinate of a neuron
};

#endif