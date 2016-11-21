#include "PoolParallel.h"
#include <mpi.h>

using namespace std;

int main(int argc, char** argv)
{
	double Gei_mean = 0.15;
	double Gei_var = 0.0;
	double Gie_var = 0.0;

    double interval = 1000;
    double timeStep = 0.01;
    int trials = 100000;

    double network_update, Ei, beta, beta_s, Ap, Ad, Ap_super, Ad_super, f0, activation, super_threshold, Gmax, Gie_mean, Tp, Td, tauP, tauD;
    int N_RA, num_inh_clusters_in_row, num_inh_in_cluster, N_ss, N_TR;
	double a, b, s_rai, s_ira;

	double maturation_threshold; // maturation threshold for a neuron

	double sigma_soma; // white noise amplitude in soma compartment
	double sigma_dend; // white noise amplitude in dendritic compartment
	double mu_soma; // white noise mean in soma
	double mu_dend; // white noise mean in dend

    string outputDirectory;
	string filenumber; // string containing number of the files from which to read network connectivity
	

    int rank; // MPI process rank
	int reading; // indicator whether we read connections from the files or not
	int testing; // indicator if it is a test run
	int training; // indcator if we innervate training neurons
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
        Gie_mean = atof(argv[5]);
		Ei = atof(argv[6]);
        beta = atof(argv[7]);
        beta_s = atof(argv[8]);
		Tp = atof(argv[9]);
		Td = atof(argv[10]);
		tauP = atof(argv[11]);
		tauD = atof(argv[12]);
        Ap = atof(argv[13]);
        Ad = atof(argv[14]);
        Ap_super = atof(argv[15]);
        Ad_super = atof(argv[16]);
		f0 = atof(argv[17]);
        activation = atof(argv[18]);
        super_threshold = atof(argv[19]);
        maturation_threshold = atof(argv[20]);
        Gmax = atof(argv[21]);
        N_RA = atoi(argv[22]);
        num_inh_clusters_in_row = atoi(argv[23]);
        num_inh_in_cluster = atoi(argv[24]);
        N_ss = atoi(argv[25]);
        N_TR = atoi(argv[26]);
        outputDirectory = argv[27];
		reading = atoi(argv[28]);
		filenumber = argv[29];
		testing = atoi(argv[30]);
		training = atoi(argv[31]);
        network_update = atof(argv[32]);
		a = atof(argv[33]);
		s_rai = atof(argv[34]);
		b = atof(argv[35]);
		s_ira = atof(argv[36]);

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
    string fileTimeInterneuron = outputDirectory + "time_info_interneuron.bin";

    string fileWeights = outputDirectory + "weights.bin";
    string fileWeightsPerm = outputDirectory;
    string RAdir = outputDirectory + "RAneurons/";
    string Idir = outputDirectory + "Ineurons/";

    int N_I = num_inh_clusters_in_row * num_inh_clusters_in_row * num_inh_in_cluster;
	
	PoolParallel pool(a, s_rai, b, s_ira, network_update, Ei, beta, beta_s, Tp, Td, tauP, tauD, Ap, Ad, Ap_super, Ad_super, 
					  f0, activation, super_threshold, maturation_threshold, Gmax, N_RA, num_inh_clusters_in_row, num_inh_in_cluster, N_ss, N_TR);

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
	
    //double mu_soma = 30;
    //double sigma_soma = 100;
    //double mu_dend = 50;
    //double sigma_dend = 185;


    //pool.print_invariable_connections();

    pool.set_generator4neurons();
    pool.set_dynamics(interval, timeStep);
    pool.set_white_noise_RA(mu_soma, sigma_soma, mu_dend, sigma_dend);

	pool.print_simulation_parameters();
    
	string weightsFilename, pajekSuperFilename, pajekActiveFilename, pajekAllFilename, fileAllRAneurons, fileAllIneurons, fileMature;
   
   	int synapses_trials_update = 20;
	int weights_trials_update = 70;

	pool.write_sim_info(fileSimInfo.c_str(), synapses_trials_update, weights_trials_update);
	
	
    double start_time = MPI_Wtime();
    bool data_gathered; // indicator if data was already gathered

	while (true)
    {
		//break;
		data_gathered = false;
		//break;
        pool.trial(training);
        //pool.gather_data();
	    
        int trial_number = pool.get_trial_number();

        if (rank == 0)
            printf("Trial %d\n", trial_number);
        
        //pool.statistics();

		std::vector<int> RAtoWrite{0, 1, 2, 3, 65, 275, 107, 67, 189, 16, 77};

        if (trial_number % synapses_trials_update == 0)
        {
			pool.gather_data();
			data_gathered = true;

           	pool.write_num_synapses(fileSynapticInfo.c_str());
		    pool.write_soma_time_info(fileTimeSoma.c_str());
            pool.write_dend_time_info(fileTimeDend.c_str());
            pool.write_interneuron_time_info(fileTimeInterneuron.c_str());
            pool.write_weights(fileWeights.c_str());
            pool.write_active_synapses(fileActive.c_str());
	    	pool.write_supersynapses(fileSuper.c_str());
            pool.write_RA(fileRA.c_str(), 3);
            pool.write_I(fileI.c_str(), 1);
            pool.write_pajek_super(filePajekSuper.c_str());
            pool.write_pajek_active(filePajekActive.c_str());
           
			//for (int i = 0; i < (int) RAtoWrite.size(); i++)
	    	//{
	  		//	fileAllRAneurons = RAdir + "RA" + std::to_string(RAtoWrite[i]) + ".bin";
			//	pool.write_RA(fileAllRAneurons.c_str(), RAtoWrite[i]);
	    	//}

			//fileMature = fileMaturePerm + "mature" + std::to_string(count) + ".bin";
			//pool.write_mature(fileMature.c_str());
			//pool.write_pajek_all(filePajekAll.c_str());
        }
	
			

		if (trial_number % weights_trials_update == 0)
		{
			if (!data_gathered)
				pool.gather_data();

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

		
			
	    	count++;
		}
		
		//break;    
		//for (int i = 0; i < N_I; i++)
	    //{
	      //  fileAllIneurons = Idir + "I" + std::to_string(i) + ".bin";
		//	pool.write_I(fileAllIneurons.c_str(), i);
	    //}


        pool.randomize_after_trial();
		
		//break;
    }
	
	
    double end_time = MPI_Wtime();

    if (rank == 0)
        printf("Execution time = %f\n", end_time - start_time);

	MPI_Finalize();


	return 0;

}
