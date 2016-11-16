#include "PoolParallel.h"
#include <mpi.h>

using namespace std;

int main(int argc, char** argv)
{
    double interval = 500;
    double timeStep = 0.01;
    int num_trials = 0;

    double network_update = 0.01;
	double Ei = -55;
	double beta = 0;
	double beta_s = 0;
	double Ap = 0;
	double Ad = 0;
	double Ap_super = 0;
	double Ad_super = 0;
	double f0 = 0;
	double activation = 0;
	double super_threshold = 0;
	double Gmax = 0;
	double Gie_mean = 0;
	double Tp = 0;
	double Td = 0;
	double tauP = 0;
	double tauD = 0;

    int N_RA = 0;

	int num_inh_clusters_in_row = 0;
	int num_inh_in_cluster = 0;
	int N_ss = 0;
	int N_TR = 0;;
	double a = 0;
	double b = 0;
	double s_rai = 0;
	double s_ira = 0;;

	double maturation_threshold = 0; // maturation threshold for a neuron

	double sigma_soma; // white noise amplitude in soma compartment
	double sigma_dend; // white noise amplitude in dendritic compartment
	double mu_soma; // white noise mean in soma
	double mu_dend; // white noise mean in dend

    string outputDirectory;
	string filenumber; // string containing number of the files from which to read network connectivity
	

    int rank; // MPI process rank
	int count; // file number from which to start writing output data

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // parse command line parameters
    //
    if (argc > 1)
    {
		mu_soma = atof(argv[1]);
		sigma_soma = atof(argv[2]);
		mu_dend = atof(argv[3]);
		sigma_dend = atof(argv[4]);
        N_RA = atoi(argv[22]);
        num_inh_clusters_in_row = atoi(argv[23]);
        num_inh_in_cluster = atoi(argv[24]);
        N_ss = atoi(argv[25]);
        N_TR = atoi(argv[26]);
        outputDirectory = argv[27];
		filenumber = argv[29];

        if (rank == 0)
            printf("Output directory is %s\n", outputDirectory.c_str());
    }
    else
    {
        printf("Command line parameters are needed!\n");
        return -1;
    }

    string fileRA = outputDirectory + "RA.bin";
    string fileI = outputDirectory + "I.bin";
    string fileIxy = outputDirectory + "I_xy.bin";
    string fileRAxy = outputDirectory + "RA_xy.bin";
    string fileRA2I = outputDirectory + "RA_I_connections.bin";
    string fileI2RA = outputDirectory + "I_RA_connections.bin";
    string fileTimeInfo = outputDirectory + "timeInfo.bin";
    
	string fileSimInfo = outputDirectory + "sim_info.bin";
	string fileSynapticInfo = outputDirectory + "synaptic_info.bin";
	string fileMaturePerm = outputDirectory;
	
	string fileMatureInfo = outputDirectory + "mature" + filenumber + ".bin"; // file from which to read mature information
	string fileAllInfo = outputDirectory + "weights" + filenumber + ".bin"; // file from which to read all RA-RA connections
	string fileActiveInfo = outputDirectory + "active" + filenumber + ".net"; // file from which to read all active RA-RA connections
	string fileSuperInfo = outputDirectory + "super" + filenumber + ".net"; // file from which to read all super RA-RA connections

    string fileActive = outputDirectory + "RA_RA_connections.bin";
    string fileSuper = outputDirectory + "RA_RA_super_connections.bin";
    string fileTimeSoma = outputDirectory + "time_info_soma.bin";
    string fileTimeDend = outputDirectory + "time_info_dend.bin";
    string fileTimeInterneuron = outputDirectory + "time_info_interneuron.bin";

    string fileChainTest = outputDirectory + "chain_test.bin";
    
	string RAdir = outputDirectory + "RAneurons/";
    string Idir = outputDirectory + "Ineurons/";

    int N_I = num_inh_clusters_in_row * num_inh_clusters_in_row * num_inh_in_cluster;
	
	PoolParallel pool(a, s_rai, b, s_ira, network_update, Ei, beta, beta_s, Tp, Td, tauP, tauD, Ap, Ad, Ap_super, Ad_super, 
					  f0, activation, super_threshold, maturation_threshold, Gmax, N_RA, num_inh_clusters_in_row, num_inh_in_cluster, N_ss, N_TR);

	pool.initialize_generator();


		pool.read_from_file(fileRAxy.c_str(), fileIxy.c_str(), fileAllInfo.c_str(), fileActiveInfo.c_str(), fileSuperInfo.c_str(), fileRA2I.c_str(), fileI2RA.c_str(), fileMatureInfo.c_str(), fileTimeInfo.c_str());

		pool.send_connections();
		pool.send_simulation_parameters();
		
		pool.update_synaptic_info();
	
    //pool.print_invariable_connections();

    pool.set_generator4neurons();
    pool.set_dynamics(interval, timeStep);
    pool.set_white_noise_RA(mu_soma, sigma_soma, mu_dend, sigma_dend);

	pool.print_simulation_parameters();
    
    double start_time = MPI_Wtime();
	pool.mature_chain_test(num_trials, fileTimeSoma.c_str(), fileTimeDend.c_str(), fileChainTest.c_str());
	
    double end_time = MPI_Wtime();

    if (rank == 0)
        printf("Execution time = %f\n", end_time - start_time);

	MPI_Finalize();


	return 0;

}
