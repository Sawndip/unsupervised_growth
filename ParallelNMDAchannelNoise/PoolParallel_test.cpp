#include "PoolParallel.h"
#include <mpi.h>

using namespace std;

int main(int argc, char** argv)
{
    bool training = true;

	double Gei_mean = 0.15;
	double Gei_var = 0.0;
	double Gie_mean = 0.60;
	double Gie_var = 0.0;

    double interval = 1000;
    double timeStep = 0.01;
    int trials = 100000;

    string fileRA = "/home/eugene/Output/RA.bin";
    string fileI = "/home/eugene/Output/I.bin";
    string filePajekSuper = "/home/eugene/Output/super.net";
    string filePajekSuperPerm = "/home/eugene/Output/";
    string filePajekActive = "/home/eugene/Output/active.net";
    string filePajekActivePerm = "/home/eugene/Output/";
    string filePajekAll = "/home/eugene/Output/all.net";
    string filePajekAllPerm = "/home/eugene/Output/";
    string filePajekFixed = "/home/eugene/Output/fixed.net";
    string fileIxy = "/home/eugene/Output/I_xy.bin";
    string fileRAxy = "/home/eugene/Output/RA_xy.bin";
    string fileRA2I = "/home/eugene/Output/RA_I_connections.bin";
    string fileI2RA = "/home/eugene/Output/I_RA_connections.bin";
    
	string fileSimInfo = "/home/eugene/Output/sim_info.bin";
	string fileSynapticInfo = "/home/eugene/Output/synaptic_info.bin";
	string fileMaturePerm = "/home/eugene/Output/";
	
	string fileMatureInfo = "/home/eugene/Output/mature393.bin"; // file from which to read mature information
	string fileAllInfo = "/home/eugene/Output/all393.net"; // file from which to read all RA-RA connections
	string fileActiveInfo = "/home/eugene/Output/active393.net"; // file from which to read all active RA-RA connections
	string fileSuperInfo = "/home/eugene/Output/super393.net"; // file from which to read all super RA-RA connections

    string fileActive = "/home/eugene/Output/RA_RA_connections.bin";
    string fileSuper = "/home/eugene/Output/RA_RA_super_connections.bin";
    string fileTimeSoma = "/home/eugene/Output/time_info_soma.bin";
    string fileTimeDend = "/home/eugene/Output/time_info_dend.bin";

    string fileWeights = "/home/eugene/Output/weights.bin";
    string fileWeightsPerm = "/home/eugene/Output/";
    string RAdir = "/home/eugene/Output/RAneurons/";
    string Idir = "/home/eugene/Output/Ineurons/";

    int rank;
    int N_RA = 100;
    //int N_I = 36;
	int num_inh_clusters_in_row = 4;
	int num_inh_in_cluster = 2;
    int N_TR = 2;
    int N_ss = 2;

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	PoolParallel pool(N_TR, N_RA, num_inh_clusters_in_row, num_inh_in_cluster, N_ss);;


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
           // pool.write_weights(fileWeights.c_str());
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
			pool.write_weights(fileWeights.c_str());

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
