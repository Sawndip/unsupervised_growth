#include "PoolParallel.h"
#include <mpi.h>

using namespace std;

int main(int argc, char** argv)
{
    bool training = true;

	double Gei_mean = 0.15;
	double Gei_var = 0.0;
	double Gie_var = 0.0;

    double interval = 1000;
    double timeStep = 0.01;
    int trials = 100000;

    double beta, beta_s, Ap, Ad, Ap_super, Ad_super, activation, super_threshold, Gmax, Gie_mean, Tp, Td, tauP, tauD;
    int N_RA, num_inh_clusters_in_row, num_inh_in_cluster, N_ss, N_TR;

    string outputDirectory;
    string workDirectory = "/home/eugene/";

    int rank;

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // parse command line parameters
    //
    if (argc > 1)
    {
        Gie_mean = atof(argv[1]);
        beta = atof(argv[2]);
        beta_s = atof(argv[3]);
		Tp = atof(argv[4]);
		Td = atof(argv[5]);
		tauP = atof(argv[6]);
		tauD = atof(argv[7]);
        Ap = atof(argv[8]);
        Ad = atof(argv[9]);
        Ap_super = atof(argv[10]);
        Ad_super = atof(argv[11]);
        activation = atof(argv[12]);
        super_threshold = atof(argv[13]);
        Gmax = atof(argv[14]);
        N_RA = atoi(argv[15]);
        num_inh_clusters_in_row = atoi(argv[16]);
        num_inh_in_cluster = atoi(argv[17]);
        N_ss = atoi(argv[18]);
        N_TR = atoi(argv[19]);
        outputDirectory = argv[20];

        
        if (rank == 0)
            printf("Output directory is %s\n", outputDirectory.c_str());
    }
    else
    {
        printf("Command line parameters are needed!\n");
        return -1;
    }

    string fileRA = workDirectory + outputDirectory + "RA.bin";
    string fileI = workDirectory + outputDirectory + "I.bin";
    string filePajekSuper = workDirectory + outputDirectory + "super.net";
    string filePajekSuperPerm = workDirectory + outputDirectory;
    string filePajekActive = workDirectory + outputDirectory + "active.net";
    string filePajekActivePerm = workDirectory + outputDirectory;
    string filePajekAll = workDirectory + outputDirectory + "all.net";
    string filePajekAllPerm = workDirectory + outputDirectory;
    string filePajekFixed = workDirectory + outputDirectory + "fixed.net";
    string fileIxy = workDirectory + outputDirectory + "I_xy.bin";
    string fileRAxy = workDirectory + outputDirectory + "RA_xy.bin";
    string fileRA2I = workDirectory + outputDirectory + "RA_I_connections.bin";
    string fileI2RA = workDirectory + outputDirectory + "I_RA_connections.bin";
    
	string fileSimInfo = workDirectory + outputDirectory + "sim_info.bin";
	string fileSynapticInfo = workDirectory + outputDirectory + "synaptic_info.bin";
	string fileMaturePerm = workDirectory + outputDirectory;
	
	string fileMatureInfo = workDirectory + outputDirectory + "mature23.bin"; // file from which to read mature information
	string fileAllInfo = workDirectory + outputDirectory + "all23.net"; // file from which to read all RA-RA connections
	string fileActiveInfo = workDirectory + outputDirectory + "active23.net"; // file from which to read all active RA-RA connections
	string fileSuperInfo = workDirectory + outputDirectory + "super23.net"; // file from which to read all super RA-RA connections

    string fileActive = workDirectory + outputDirectory + "RA_RA_connections.bin";
    string fileSuper = workDirectory + outputDirectory + "RA_RA_super_connections.bin";
    string fileTimeSoma = workDirectory + outputDirectory + "time_info_soma.bin";
    string fileTimeDend = workDirectory + outputDirectory + "time_info_dend.bin";

    string fileWeights = workDirectory + outputDirectory + "weights.bin";
    string fileWeightsPerm = workDirectory + outputDirectory;
    string RAdir = workDirectory + outputDirectory + "RAneurons/";
    string Idir = workDirectory + outputDirectory + "Ineurons/";


    //printf("My rank is %d\n", rank);


	PoolParallel pool(beta, beta_s, Tp, Td, tauP, tauD, Ap, Ad, Ap_super, Ad_super, activation, super_threshold, Gmax, N_RA, num_inh_clusters_in_row, num_inh_in_cluster, N_ss, N_TR);

	pool.print_simulation_parameters();
	pool.initialize_generator();
	pool.initialize_inhibitory_clusters();
	pool.initialize_RA_for_inh_clusters();
	//pool.initialize_coordinates();
	//pool.initialize_equal_clusters();
	//pool.initialize_connections_for_clusters(Gei_mean, Gei_var, Gie_mean, Gie_var);
	pool.initialize_connections_for_inhibitory_clusters(Gei_mean, Gei_var, Gie_mean, Gie_var);
	//pool.read_from_file(fileRAxy.c_str(), fileIxy.c_str(), fileAllInfo.c_str(), fileActiveInfo.c_str(), fileSuperInfo.c_str(), fileRA2I.c_str(), fileI2RA.c_str(), fileMatureInfo.c_str());
	
	pool.send_connections();
	//pool.send_RARA_connections();
	
	//pool.initialize_test_allI2RA_connections(0.2);
    //pool.initialize_test_allRA2I_connections(0.1);

    double mu_soma = 30;
    double sigma_soma = 100;
    double mu_dend = 50;
    double sigma_dend = 185;


    //pool.print_invariable_connections();

    pool.set_generator4neurons();
    pool.set_dynamics(interval, timeStep);
    pool.set_no_noise_RA();
    pool.set_no_noise_I();
    pool.set_white_noise_distribution_soma(mu_soma, sigma_soma);
    pool.set_white_noise_distribution_dend(mu_dend, sigma_dend);
    pool.set_white_noise_soma();
    pool.set_white_noise_dend();

	//pool.update_synaptic_info();
    pool.write_coordinates(fileRAxy.c_str(), fileIxy.c_str());
    pool.write_invariable_synapses(fileRA2I.c_str(), fileI2RA.c_str());
    pool.write_pajek_fixed(filePajekFixed.c_str());
    
   // pool.write_pajek_active("/home/eugene/Output/active_test.net");
	//pool.write_pajek_super("/home/eugene/Output/super_test.net");

    //pool.write_pajek_all("/home/eugene/Output/all_test.net");
	//pool.write_mature("/home/eugene/Output/mature_test.net");

    //double start_time = MPI_Wtime();
    string weightsFilename, pajekSuperFilename, pajekActiveFilename, pajekAllFilename, fileAllRAneurons, fileAllIneurons, fileMature;
    int count = 0;
   
   	int synapses_trials_update = 15;
	int weights_trials_update = 50;

	pool.write_sim_info(fileSimInfo.c_str(), synapses_trials_update, weights_trials_update);
	
	
    for (int i = 0; i < trials; i++)
    {
        if (rank == 0)
            printf("Trial %d\n", i);
        pool.trial(training);
        pool.gather_data();
	//pool.statistics();

        if (i % synapses_trials_update == 0)
        {
           	pool.write_num_synapses(fileSynapticInfo.c_str());
		    pool.write_soma_time_info(fileTimeSoma.c_str());
            pool.write_dend_time_info(fileTimeDend.c_str());
            pool.write_weights(fileWeights.c_str());
            pool.write_active_synapses(fileActive.c_str());
	    	pool.write_supersynapses(fileSuper.c_str());
            pool.write_RA(fileRA.c_str(), 3);
            pool.write_I(fileI.c_str(), 1);
            pool.write_pajek_super(filePajekSuper.c_str());
            pool.write_pajek_active(filePajekActive.c_str());
            //pool.write_pajek_all(filePajekAll.c_str());
        }
	
			

		if (i % weights_trials_update == 0)
		{
	 		weightsFilename = fileWeightsPerm + "weights" + std::to_string(count) + ".bin";
			pool.write_weights(weightsFilename.c_str());
			//pool.write_weights(fileWeights.c_str());

			pajekSuperFilename = filePajekSuperPerm + "super" + std::to_string(count) + ".net";
	    	pool.write_pajek_super(pajekSuperFilename.c_str());
	    
	   		pajekActiveFilename = filePajekActivePerm + "active" + std::to_string(count) + ".net";
	    	pool.write_pajek_active(pajekActiveFilename.c_str());

	    	pajekAllFilename = filePajekAllPerm + "all" + std::to_string(count) + ".net";
	    	pool.write_pajek_all(pajekAllFilename.c_str());

			fileMature = fileMaturePerm + "mature" + std::to_string(count) + ".bin";
			pool.write_mature(fileMature.c_str());

		
			//for (int i = 0; i < N_RA; i++)
	    	//{
	        //	fileAllRAneurons = RAdir + "RA" + std::to_string(i) + ".bin";
			//	pool.write_RA(fileAllRAneurons.c_str(), i);
	    	//}

	    	//for (int i = 0; i < N_I; i++)
	    	//{
	        //	fileAllIneurons = Idir + "I" + std::to_string(i) + ".bin";
			//	pool.write_I(fileAllIneurons.c_str(), i);
	    	//}


	    	count++;
		}

        pool.randomize_after_trial();
    }
	
	
    //double end_time = MPI_Wtime();

    //if (rank == 0)
      //  printf("Execution time = %f\n", end_time - start_time);

	MPI_Finalize();


	return 0;

}
