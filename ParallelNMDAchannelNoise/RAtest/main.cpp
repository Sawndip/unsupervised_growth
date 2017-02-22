#include "../HH2_final.h"
#include "../HHI_final.h"
#include "../poisson_noise.h"
#include <string>

#define start 250
#define duration 5
#define ampl 20.0

double I(double t)
{
    if ( (t >= start) && (t <= duration + start) )
        return ampl;
    else
        return 0.0;
}

using namespace std;

int main(int argc, char** argv)
{
    std::string outFileRA = "/home/eugene/forFigures/RA_inh_input.bin";
    std::string outFileI = "/home/eugene/forFigures/I.bin";

    HH2_final n;
    double trial_duration = 500; // trial duration in ms
    double timestep = 0.02; // dynamics timestep
    
    // noise
    double mean_soma = 175.0;
    double std_soma = 50.0;
    double mean_dend = 175.0;
    double std_dend = 100.0;

    // noise generator
    Poisson_noise noise_generator;

    n.set_noise_generator(&noise_generator);
    n.set_white_noise(mean_soma, std_soma, mean_dend, std_dend);
    n.set_dynamics(trial_duration, timestep);

    // set GABA potential

    double Ei = -80.0; // reverse GABA potential

    n.set_Ei(Ei);

    // create inhibitory neuron
    HHI_final nInh;

    nInh.set_dynamics(trial_duration, timestep);

    // connect interneuron to HVC(RA) neuron
    double Gie = 1.5;

    nInh.set_target(&n, 0, Gie);

    // set injection current to interneuron
    nInh.set_injected_current(&I);


    int size = static_cast<int>(trial_duration / timestep);

    for (int i = 0; i < size; i++)
    {
        n.Debraband_step_no_target_update();
        nInh.DP8_step_with_target_update();
    }
    n.writeToFile(outFileRA.c_str());
    nInh.writeToFile(outFileI.c_str());

	return 0;

}
