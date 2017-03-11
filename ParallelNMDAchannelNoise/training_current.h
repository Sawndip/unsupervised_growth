#pragma once
#ifndef TRAINING_CURRENT
#define TRAINING_CURRENT

double training_current(double start, double t)
{
	//double start = 0;
	double duration = 10;
	double I_ampl = 1500;

	if ((t >= start)&&(t <= (start + duration)))
		return I_ampl;
	else
		return 0;
}

#endif
