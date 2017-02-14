#include "../Configuration.h"
#include <string>

int main()
{
    std::string filename = "/home/eugene/projects/chaingrowth/ParallelNMDAchannelNoise/config/parameters_written.cfg";
    std::string outfile = "/home/eugene/projects/chaingrowth/ParallelNMDAchannelNoise/config/parameters_rewritten.cfg";

    struct NetworkParameters network_params;
    struct SynapticParameters synaptic_params;
    struct SpatialParameters spatial_params;
    struct GabaParameters gaba_params;
    struct TimeParameters time_params;
    struct NoiseParameters noise_params;

    Configuration config;

    config.read_configuration(filename.c_str());
    config.print_configuration();

    config.write_configuration(outfile.c_str());

    //read_parameters(filename.c_str(), &network_params, &synaptic_params, &spatial_params, &gaba_params, &time_params, &noise_params);



    return 0;
}
