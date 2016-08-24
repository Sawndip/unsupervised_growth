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

    int N_I = num_inh_clusters_in_row * num_inh_clusters_in_row * num_inh_in_cluster;

	double sigma_soma; // white noise amplitude in soma compartment
	double sigma_dend; // white noise amplitude in dendritic compartment

    string outputDirectory;
	string filenumber; // string containing number of the files from which to read network connectivity
	

    int rank; // MPI process rank
	int reading; // indicator whether we read connections from the files or not
	int testing; // indicator if it is a test run
	int count; // file number from which to start writing output data

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // parse command line parameters
    //
    if (argc > 1)
    {
		sigma_soma = atof(argv[1]);
		sigma_dend = atof(argv[2]);
        Gie_mean = atof(argv[3]);
        beta = atof(argv[4]);
        beta_s = atof(argv[5]);
		Tp = atof(argv[6]);
		Td = atof(argv[7]);
		tauP = atof(argv[8]);
		tauD = atof(argv[9]);
        Ap = atof(argv[10]);
        Ad = atof(argv[11]);
        Ap_super = atof(argv[12]);
        Ad_super = atof(argv[13]);
        activation = atof(argv[14]);
        super_threshold = atof(argv[15]);
        Gmax = atof(argv[16]);
        N_RA = atoi(argv[17]);
        num_inh_clusters_in_row = atoi(argv[18]);
        num_inh_in_cluster = atoi(argv[19]);
        N_ss = atoi(argv[20]);
        N_TR = atoi(argv[21]);
        outputDirectory = argv[22];
		reading = atoi(argv[23]);
		filenumber = argv[24];
		testing = atoi(argv[25]);

        
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
    string filePajekSuper = outputDirectory + "super.net";
    string filePajekSuperPerm = outputDirectory;
    string filePajekActive = outputDirectory + "active.net";
    string filePajekActivePerm = outputDirectory;
    string filePajekAll = outputDirectory + "all.net";
    string filePajekAllPerm = outputDirectory;
    string filePajekFixed = outputDirectory + "fixed.net";
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

    string fileWeights = outputDirectory + "weights.bin";
    string fileWeightsPerm = outputDirectory;
    string RAdir = outputDirectory + "RAneurons/";
    string Idir = outputDirectory + "Ineurons/";

	PoolParallel pool(beta, beta_s, Tp, Td, tauP, tauD, Ap, Ad, Ap_super, Ad_super, activation, super_threshold, Gmax, N_RA, num_inh_clusters_in_row, num_inh_in_cluster, N_ss, N_TR);

	pool.print_simulation_parameters();
	pool.initialize_generator();

	// if we need to read connections from network
	if (reading > 0)
	{
		count = atoi(filenumber.c_str()) + 1;

		pool.read_from_file(fileRAxy.c_str(), fileIxy.c_str(), fileAllInfo.c_str(), fileActiveInfo.c_str(), fileSuperInfo.c_str(), fileRA2I.c_str(), fileI2RA.c_str(), fileMatureInfo.c_str(), fileTimeInfo.c_str());

		pool.send_connections();
		pool.send_simulation_parameters();
		
		pool.update_synaptic_info();
	}
	else
	{
		count = 1;

		// if testing run assemble neurons into clusters
		if (testing > 0)
		{

			pool.initialize_inhibitory_clusters();
			pool.initialize_RA_for_inh_clusters();
			pool.initialize_connections_for_inhibitory_clusters(Gei_mean, Gei_var, Gie_mean, Gie_var);

		}
		else
		// if real run set distribution for connections between neurons
		{
			pool.initialize_coordinates();
			pool.initialize_connections(Gei_mean, Gei_var, Gie_mean, Gie_var);
		}

		pool.send_connections();
    	pool.write_coordinates(fileRAxy.c_str(), fileIxy.c_str());
    	pool.write_invariable_synapses(fileRA2I.c_str(), fileI2RA.c_str());
    	pool.write_pajek_fixed(filePajekFixed.c_str());
	}
	
    double mu_soma = 30;
    //double sigma_soma = 100;
    double mu_dend = 50;
    //double sigma_dend = 185;


    //pool.print_invariable_connections();

    pool.set_generator4neurons();
    pool.set_dynamics(interval, timeStep);
    pool.set_no_noise_RA();
    pool.set_no_noise_I();
    pool.set_white_noise_distribution_soma(mu_soma, sigma_soma);
    pool.set_white_noise_distribution_dend(mu_dend, sigma_dend);
    pool.set_white_noise_soma();
    pool.set_white_noise_dend();

    //double start_time = MPI_Wtime();
    string weightsFilename, pajekSuperFilename, pajekActiveFilename, pajekAllFilename, fileAllRAneurons, fileAllIneurons, fileMature;
   
   	int synapses_trials_update = 5;
	int weights_trials_update = 10;

	pool.write_sim_info(fileSimInfo.c_str(), synapses_trials_update, weights_trials_update);
	
	
    while (true)
    {
        pool.trial(training);
        pool.gather_data();
	    
        int trial_number = pool.get_trial_number();

        if (rank == 0)
            printf("Trial %d\n", trial_number);
        
        //pool.statistics();

        if (trial_number % synapses_trials_update == 0)
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
	
			

		if (trial_number % weights_trials_update == 0)
		{
	 		weightsFilename = fileWeightsPerm + "weights" + std::to_string(count) + ".bin";
			pool.write_weights(weightsFilename.c_str());
			//pool.write_weights(fileWeights.c_str());

			pajekSuperFilename = filePajekSuperPerm + "super" + std::to_string(count) + ".net";
	    	pool.write_pajek_super(pajekSuperFilename.c_str());
	    
	   		pajekActiveFilename = filePajekActivePerm + "active" + std::to_string(count) + ".net";
	    	pool.write_pajek_active(pajekActiveFilename.c_str());

	    	//pajekAllFilename = filePajekAllPerm + "all" + std::to_string(count) + ".net";
	    	//pool.write_pajek_all(pajekAllFilename.c_str());

			fileMature = fileMaturePerm + "mature" + std::to_string(count) + ".bin";
			pool.write_mature(fileMature.c_str());

            pool.write_time_info(fileTimeInfo.c_str());

		
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
