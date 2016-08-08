#pragma once
#ifndef POISSON_NOISE.H
#define POISSON_NOISE.H
#endif

#include <random>
#include <tgmath.h>

class Poisson_noise
{
public:
	Poisson_noise() : d(0.0,1.0), mean(0), sd(1){};
	void set_seed(unsigned s);	//	function to set seed to randon number generator
	double get_spike_time(double lambda);	//	get time for the next noisy spike
	double random(double G);	//	get random number in range (0; G)
	void set_normal_distribution(double mu, double sigma); // set parameters of normal distribution
	double normal_distribution(); // get number sampled from normal distribution
private:
	unsigned seed;	//	seed for random number generator
	std::mt19937 generator;	//	marsene-twister generator
	std::normal_distribution<double> d; // normal distribution with zero mean and unit variance
	double mean; // mean for normal distribution
	double sd; // standard deviation for noraml distribution
};