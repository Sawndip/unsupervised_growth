# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 19:27:45 2018

@author: jingroup

Script determines inhibition received before the neuron burst and compares it to inhibition received
by other neurons at the same time
"""
import reading
import utils
import numpy as np
import os
import matplotlib.pyplot as plt
import sklearn.preprocessing
import Tkinter as tk

def maximize_figure(figure=None):
    root = tk.Tk()
    width = root.winfo_screenmmwidth() / 25.4
    height = root.winfo_screenmmheight() / 25.4
    if figure is not None:
        plt.figure(figure).set_size_inches(width, height)
    return width, height

def compare_networkAndPoolConductance(dataDir, testDataDir, outFigureDir, simName, trial):
    """
    Function reads conductances for all neurons in the network, calculates average total conductances
    and conductances alinged to average bursting times
    """
    N_RA, N_I = reading.read_num_neurons(os.path.join(dataDir, "num_neurons.bin"))
    training_neurons = reading.read_training_neurons(os.path.join(dataDir, "training_neurons.bin"))
    N_TR = len(training_neurons)
    
    _, numTestTrials, \
        probability_soma_spike, average_num_soma_spikes_in_trial, mean_first_soma_spike_time, std_first_soma_spike_time,\
        probability_dend_spike, average_num_dend_spikes_in_trial, mean_first_dend_spike_time, std_first_dend_spike_time = reading.read_jitter(os.path.join(testDataDir, "jitter.bin")) 
    
    neuronsWithRobustDendriticSpike = np.where(probability_dend_spike >= 0.75)[0]
    meanFirstSomaSpikeOfNeuronsWithRoburstDendriticSpike = mean_first_soma_spike_time[neuronsWithRobustDendriticSpike]
    
    t, Vs, Vd, Gexc_d, Ginh_d, n, h, r, c, Ca = reading.read_hh2_buffer_full(os.path.join(testDataDir, "RA/RA0_trial0.bin"))
    
    Ginh = np.zeros((N_RA, len(t)), np.float32)
    Gexc = np.zeros((N_RA, len(t)), np.float32)
    
    numOtherNeurons = N_RA - 1
    
    GinhSumAll = np.zeros(len(t), np.float32)
    GexcSumAll = np.zeros(len(t), np.float32)
    
    for testTrial in range(numTestTrials):
        print "Test trial: ",testTrial
        #if testTrial == 1:
          #  break
        for neuronId in range(N_RA):
            t, Vs, Vd, Gexc_d, Ginh_d, n, h, r, c, Ca = reading.read_hh2_buffer_full(os.path.join(testDataDir, "RA/RA"+str(neuronId)+"_trial" + str(testTrial) + ".bin"))
            
            Ginh[neuronId] += Ginh_d
            Gexc[neuronId] += Gexc_d
        
            # sum of conductances for all neurons excluding training
            if neuronId not in training_neurons:    
                GinhSumAll += Ginh_d 
                GexcSumAll += Gexc_d 
                
    for neuronId in range(N_RA):
        Ginh[neuronId] /= float(numTestTrials)
        
    GinhSumAll = GinhSumAll / (float(numTestTrials) * float(N_RA - N_TR))
    GexcSumAll = GexcSumAll / (float(numTestTrials) * float(N_RA - N_TR))
           
    #print np.max(spike_times_d)
    #print np.max(t)
    
    #window = 100.0 # window size in ms
    window = 50.0
    
    Gbursted = None # conductance aligned to bursting time
    Gother = None # conductance of neurons that did npt burst
     
    
    dt = t[1] - t[0]
    window_t = [float(i)*dt-window/2. for i in range(int(window / dt)-1)]
    
    GburstedInh_window = np.empty((len(neuronsWithRobustDendriticSpike) - N_TR, int(window / dt)-1), np.float32) # conductances of all burst neurons aligned to burst time
    GburstedExc_window = np.empty((len(neuronsWithRobustDendriticSpike) - N_TR, int(window / dt)-1), np.float32) # conductances of all burst neurons aligned to burst time
    
    # plot conductances for several random bursted neurons
    np.random.seed(1991)
    
    nid_toPlot = np.random.choice(neuronsWithRobustDendriticSpike, 16, replace=False)
    
    nrows = 4
    ncols = 4
    
    fInh, axarrInh = plt.subplots(nrows=nrows, ncols=ncols)
    fExc, axarrExc = plt.subplots(nrows=nrows, ncols=ncols)
    
    neuronPlotCounter = 0
    neuronSavedCounter = 0
    
    for nid, meanFirstSpikeTime in zip(neuronsWithRobustDendriticSpike, meanFirstSomaSpikeOfNeuronsWithRoburstDendriticSpike):
        if nid in training_neurons:
            continue
    
        meanFirstSpikeTime =  round(int(meanFirstSpikeTime / dt) *dt, 2)
        
        GburstedInh_window[neuronSavedCounter] = Ginh[nid][(t >meanFirstSpikeTime - window/2.)&(t < meanFirstSpikeTime + window/2.)]    
        GburstedExc_window[neuronSavedCounter] = Gexc[nid][(t >meanFirstSpikeTime - window/2.)&(t < meanFirstSpikeTime + window/2.)]    
    
        # normalize conductance by max value    
        #Gbursted_window[neuronSavedCounter] /= np.max(Gbursted_window[neuronSavedCounter])    
        
        # normalize to o mean and unit variance
        #Gbursted_window[neuronSavedCounter] = sklearn.preprocessing.scale(Gbursted_window[neuronSavedCounter], axis=0, with_mean=True, with_std=True, copy=True)    
        
        if nid in nid_toPlot:
            row = neuronPlotCounter // 4
            col = neuronPlotCounter % 4 
            axarrInh[row, col].plot(window_t, GburstedInh_window[neuronSavedCounter])
            axarrExc[row, col].plot(window_t, GburstedExc_window[neuronSavedCounter])
            #axarr[row, col].vlines(meanFirstSpikeTime, 0, np.max(Ginh[nid]))
            
            if row == 3:
                axarrInh[row, col].set_xlabel('Time (ms)')
                axarrExc[row, col].set_xlabel('Time (ms)')
    
            if col == 0:
                axarrInh[row, col].set_ylabel('Ginh (mS/cm^2)')
                axarrExc[row, col].set_ylabel('Gexc (mS/cm^2)')
    
            neuronPlotCounter += 1
        
        neuronSavedCounter += 1
    
        #if Gbursted == None:
          #  Gbursted = Ginh[nid[0]][(t > dend_spike_time[0] - window/2.)&(t < dend_spike_time[0] + window/2.)]    
        #else:
          #  Gbursted += Ginh[nid[0]][(t > dend_spike_time[0] - window/2.)&(t < dend_spike_time[0] + window/2.)]    
        
         
        #if Gother == None:
          #  Gother = (GsumAll-Ginh[nid[0]])[(t > dend_spike_time[0] - window/2.)&(t < dend_spike_time[0] + window/2.)]    
        #else:
          #  Gother += (GsumAll-Ginh[nid[0]])[(t > dend_spike_time[0] - window/2.)&(t < dend_spike_time[0] + window/2.)]    
        
        #plt.figure()
        #plt.plot(t, Gexc_d)
        #plt.vlines(dend_spike_time[0], 0, np.max(Gexc_d))
        
        
        
        #plt.figure()
        #plt.plot(t, G)
        #plt.vlines(dend_spike_time[0], 0, np.max(G))
    
    w, h = maximize_figure(fInh.number)
    fInh.savefig(outFigureDir + simName + "_trial" + str(trial) + '_Ginh_examples.png', bbox_inches='tight')
    plt.close(fInh)
    
    w, h = maximize_figure(fExc.number)
    fExc.savefig(outFigureDir + simName + "_trial" + str(trial) + '_Gexc_examples.png', bbox_inches='tight')
    plt.close(fExc)
    
    GburstedInh = np.sum(GburstedInh_window, axis=0) / float(len(neuronsWithRobustDendriticSpike)-N_TR)
    std_GburstedInh = np.std(GburstedInh_window, axis=0)
    
    GburstedExc = np.sum(GburstedExc_window, axis=0) / float(len(neuronsWithRobustDendriticSpike)-N_TR)
    std_GburstedExc = np.std(GburstedExc_window, axis=0)
    
    
    f = plt.figure()
    plt.plot(window_t, GburstedInh , label='bursted neurons')
    plt.plot(window_t, GburstedInh +std_GburstedInh / np.sqrt(len(neuronsWithRobustDendriticSpike) - N_TR))
    plt.plot(window_t, GburstedInh - std_GburstedInh / np.sqrt(len(neuronsWithRobustDendriticSpike) - N_TR))
    plt.xlabel('Time (ms)')
    plt.ylabel('average Ginh (mS/cm^2)')
    
    f.savefig(outFigureDir + simName + "_trial" + str(trial) + '_average_Ginh.png', bbox_inches='tight')
    plt.close(f)
    
    #plt.plot(window_t, Gother / float((N_RA - N_TR - 1)*numNeuronsWithDendSpike), label='other neurons')
    #plt.legend()
       
    f = plt.figure()
    plt.plot(t, GinhSumAll)
    plt.xlabel('Time (ms)')
    plt.ylabel('total Ginh (mS/cm^2)')
    f.savefig(outFigureDir + simName + "_trial" + str(trial) + '_total_Ginh.png', bbox_inches='tight')
    plt.close(f)
    
    f = plt.figure()
    plt.plot(window_t, GburstedExc , label='bursted neurons')
    plt.plot(window_t, GburstedExc +std_GburstedExc / np.sqrt(len(neuronsWithRobustDendriticSpike) - N_TR))
    plt.plot(window_t, GburstedExc - std_GburstedExc / np.sqrt(len(neuronsWithRobustDendriticSpike) - N_TR))
    plt.xlabel('Time (ms)')
    plt.ylabel('average Gexc (mS/cm^2)')
    f.savefig(outFigureDir + simName + "_trial" + str(trial) + '_average_Gexc.png', bbox_inches='tight')
    plt.close(f)
    
    #plt.plot(window_t, Gother / float((N_RA - N_TR - 1)*numNeuronsWithDendSpike), label='other neurons')
    #plt.legend()
    
    f = plt.figure()
    plt.plot(t, GexcSumAll)
    plt.xlabel('Time (ms)')
    plt.ylabel('total Gexc (mS/cm^2)')
    f.savefig(outFigureDir + simName + "_trial" + str(trial) + '_total_Gexc.png', bbox_inches='tight')
    plt.close(f)
    
    
    
def compare_inhibitory_weights(dataDir, testDataDir, outFigureDir, simName, trial):
    """
    Compare inhibitory inputs to network and pool neurons
    """
    N_RA, N_I = reading.read_num_neurons(os.path.join(dataDir, "num_neurons.bin"))
    training_neurons = reading.read_training_neurons(os.path.join(dataDir, "training_neurons.bin"))
    N_TR = len(training_neurons)
    
    _, numTestTrials, \
        probability_soma_spike, average_num_soma_spikes_in_trial, mean_first_soma_spike_time, std_first_soma_spike_time,\
        probability_dend_spike, average_num_dend_spikes_in_trial, mean_first_dend_spike_time, std_first_dend_spike_time = reading.read_jitter(os.path.join(testDataDir, "jitter.bin")) 
    
    neuronsWithRobustDendriticSpike = np.where(probability_dend_spike >= 0.75)[0]
    
    # analyze inhibitory weights to neurons
    (_, targets_ID, weights_I2RA, syn_lengths, axonal_delays) = reading.read_connections(os.path.join(dataDir, "I_RA_connections_"+str(trial)+".bin"))  
    
     
    inhibition_weights_on_network_neurons = [] 
    inhibition_weights_on_pool_neurons = [] 
     
    total_inhibition_on_network_neurons = {}
    total_inhibition_on_pool_neurons = {} 
    
     
    set_neuronsWithRobustDendriticSpike = set(neuronsWithRobustDendriticSpike) 
     
    for i in range(N_I):
        for j, target in enumerate(targets_ID[i]):
            if target not in training_neurons:
                if target in set_neuronsWithRobustDendriticSpike:
                    if target in total_inhibition_on_network_neurons:
                        total_inhibition_on_network_neurons[target] += weights_I2RA[i][j]
                    else:
                        total_inhibition_on_network_neurons[target] = weights_I2RA[i][j]
    
                    inhibition_weights_on_network_neurons.append(weights_I2RA[i][j])
                else:
                    if target in total_inhibition_on_pool_neurons:
                        total_inhibition_on_pool_neurons[target] += weights_I2RA[i][j]
                    else:
                        total_inhibition_on_pool_neurons[target] = weights_I2RA[i][j]
                    inhibition_weights_on_pool_neurons.append(weights_I2RA[i][j])
    
    totalInhId, totalInhW =   total_inhibition_on_network_neurons.keys(), total_inhibition_on_network_neurons.values()  
    
    f = plt.figure()
    plt.scatter(mean_first_soma_spike_time[totalInhId], totalInhW)
    plt.xlabel('Mean first spike time (ms)')
    plt.ylabel('Total inhibitory input (mS/cm^2)')
    f.savefig(outFigureDir + simName + "_trial" + str(trial) + '_total_Ginh_vs_burstTime.png', bbox_inches='tight')
    plt.close(f)
    
    
    nbin = 20
    hist_on_network, bin_edges_on_network = np.histogram(inhibition_weights_on_network_neurons, bins=nbin)
    hist_on_network = hist_on_network.astype(float) / float(len(neuronsWithRobustDendriticSpike)-N_TR)
    bin_centers_on_network = bin_edges_on_network[:-1:1] + bin_edges_on_network[1] - bin_edges_on_network[0]
    
    hist_on_pool, bin_edges_on_pool = np.histogram(inhibition_weights_on_pool_neurons, bins=nbin)
    hist_on_pool = hist_on_pool.astype(float) / float(N_RA - len(neuronsWithRobustDendriticSpike))
    bin_centers_on_pool = bin_edges_on_pool[:-1:1] + bin_edges_on_pool[1] - bin_edges_on_pool[0]
    
    f = plt.figure()
    plt.xlabel('Inhibitory input weight (mS/cm^2)')
    plt.ylabel('# of inputs per neuron')
    #plt.hist( , fill=False, label='network neurons', edgecolor='r')
    #plt.hist(np.array(inhibition_weights_on_pool_neurons) / float(N_RA - len(neuronsWithRobustDendriticSpike) - N_TR), fill=False, label='pool neurons', edgecolor='b')
    plt.bar(bin_centers_on_network, hist_on_network, align='center',  fill=False, edgecolor='b', width=bin_edges_on_network[1] - bin_edges_on_network[0], label='network neurons')
    plt.bar(bin_centers_on_pool, hist_on_pool, align='center',  fill=False, edgecolor='r', width=bin_edges_on_pool[1] - bin_edges_on_pool[0], label='pool neurons')
    plt.legend(loc=4)
    f.savefig(outFigureDir + simName + "_trial" + str(trial) + '_Ginh_weights_comparison.png', bbox_inches='tight')
    plt.close(f)
    
    
    nbin = 20
    hist_on_network_total, bin_edges_on_network_total = np.histogram(total_inhibition_on_network_neurons.values(), bins=nbin)
    hist_on_network_total = hist_on_network_total.astype(float) / float(len(neuronsWithRobustDendriticSpike)-N_TR)
    bin_centers_on_network_total = bin_edges_on_network_total[:-1:1] + bin_edges_on_network_total[1] - bin_edges_on_network_total[0]
    
    hist_on_pool_total, bin_edges_on_pool_total = np.histogram(total_inhibition_on_pool_neurons.values(), bins=nbin)
    hist_on_pool_total = hist_on_pool_total.astype(float) / float(N_RA - len(neuronsWithRobustDendriticSpike))
    bin_centers_on_pool_total = bin_edges_on_pool_total[:-1:1] + bin_edges_on_pool_total[1] - bin_edges_on_pool_total[0]
    
    f = plt.figure()
    plt.xlabel('Total inhibitory input weight (mS/cm^2)')
    plt.ylabel('Norm. # of neurons')
    #plt.hist( , fill=False, label='network neurons', edgecolor='r')
    #plt.hist(np.array(inhibition_weights_on_pool_neurons) / float(N_RA - len(neuronsWithRobustDendriticSpike) - N_TR), fill=False, label='pool neurons', edgecolor='b')
    plt.bar(bin_centers_on_network_total, hist_on_network_total, align='center',  fill=False, edgecolor='b', width=bin_edges_on_network_total[1] - bin_edges_on_network_total[0], label='network neurons')
    plt.bar(bin_centers_on_pool_total, hist_on_pool_total, align='center',  fill=False, edgecolor='r', width=bin_edges_on_pool_total[1] - bin_edges_on_pool_total[0], label='pool neurons')
    plt.legend(loc=1)
    f.savefig(outFigureDir + simName + "_trial" + str(trial) + '_total_Ginh_comparison.png', bbox_inches='tight')
    plt.close(f)
    
    f = plt.figure()
    plt.bar([1, 2], [np.mean(total_inhibition_on_pool_neurons.values()), np.mean(total_inhibition_on_network_neurons.values())], align='center', width=0.1,
            yerr=[np.std(total_inhibition_on_pool_neurons.values()) / np.sqrt(float(N_RA - len(neuronsWithRobustDendriticSpike))), np.std(total_inhibition_on_network_neurons.values()) / np.sqrt(float(len(neuronsWithRobustDendriticSpike)-N_TR))])
    plt.xticks([1,2], ['pool', 'network'])
    plt.ylabel('Mean inhibitory input (mS/cm^2)')
    f.savefig(outFigureDir + simName + "_trial" + str(trial) + '_mean_Ginh_comparison.png', bbox_inches='tight')
    plt.close(f)
    
    from scipy.stats import ranksums
    print ranksums(total_inhibition_on_pool_neurons.values(), total_inhibition_on_network_neurons.values())
    
    plt.show()


def integral(a, dx):
    """
    Function computes integral under the curve for array a with data point resolution dx
    """
    return (np.sum(a) - (a[0]+a[-1])/2. ) * dx
    
def compare_conductanceOfNetworkNeurons(dataDir, testDataDir, outFigureDir, simName, trial):
    """
    Function reads conductances for network  neurons and compares inhibitory conductance of neurons with 
    conductances of the neurons that burst later
    """
    N_RA, N_I = reading.read_num_neurons(os.path.join(dataDir, "num_neurons.bin"))
    training_neurons = reading.read_training_neurons(os.path.join(dataDir, "training_neurons.bin"))
    N_TR = len(training_neurons)
    
    _, numTestTrials, \
        probability_soma_spike, average_num_soma_spikes_in_trial, mean_first_soma_spike_time, std_first_soma_spike_time,\
        probability_dend_spike, average_num_dend_spikes_in_trial, mean_first_dend_spike_time, std_first_dend_spike_time = reading.read_jitter(os.path.join(testDataDir, "jitter.bin")) 
    
    neuronsWithRobustDendriticSpike = np.where(probability_dend_spike >= 0.75)[0]
    meanFirstSomaSpikeOfNeuronsWithRoburstDendriticSpike = mean_first_soma_spike_time[neuronsWithRobustDendriticSpike]
    
    t, Vs, Vd, Gexc_d, Ginh_d, n, h, r, c, Ca = reading.read_hh2_buffer_full(os.path.join(testDataDir, "RA/RA0_trial0.bin"))
    
    Ginh = np.zeros((len(neuronsWithRobustDendriticSpike), len(t)), np.float32)

    
    print "Robustly firing neurons: ",neuronsWithRobustDendriticSpike    
    print "First soma spikes of robust neurons: ", meanFirstSomaSpikeOfNeuronsWithRoburstDendriticSpike
    
    for testTrial in range(numTestTrials):
        print "Test trial: ",testTrial
        #f testTrial == 1:
          #  break
        for neuronCounter, neuronId in enumerate(neuronsWithRobustDendriticSpike):
            t, Vs, Vd, Gexc_d, Ginh_d, n, h, r, c, Ca = reading.read_hh2_buffer_full(os.path.join(testDataDir, "RA/RA"+str(neuronId)+"_trial" + str(testTrial) + ".bin"))    
            Ginh[neuronCounter] += Ginh_d
                
                
    for neuronId in range(len(Ginh)):
        Ginh[neuronId] /= float(numTestTrials)
    
    window = 30.0
    dt = t[1] - t[0]
    window_t = [float(i)*dt-window/2. for i in range(int(window / dt)-1)]
    
    #GburstedInh_window = np.empty((len(neuronsWithRobustDendriticSpike) - N_TR, int(window / dt)-1), np.float32) # conductances of all burst neurons aligned to burst time
    
    # plot conductances for several random bursted neurons
    #np.random.seed(1991)
    
    #nid_toPlot = np.random.choice(neuronsWithRobustDendriticSpike, 16, replace=False)
    
    #nrows = 4
    #ncols = 4
    
    #fInh, axarrInh = plt.subplots(nrows=nrows, ncols=ncols)
    #fExc, axarrExc = plt.subplots(nrows=nrows, ncols=ncols)
    
    #neuronPlotCounter = 0
    #neuronSavedCounter = 0
    
    integralConductanceOfAllSpiked = []
    integralConductanceOfAllOther = []    
    
    for neuronCounter, (nid, meanFirstSpikeTime) in enumerate(zip(neuronsWithRobustDendriticSpike, meanFirstSomaSpikeOfNeuronsWithRoburstDendriticSpike)):
        if nid in training_neurons:
            continue
    
    
        print "counter = {0}; mean first spike time = {1}".format(neuronCounter, meanFirstSpikeTime)
        
        meanFirstSpikeTime =  round(int(meanFirstSpikeTime / dt) *dt, 2)
        
        #GburstedInh_window[neuronSavedCounter] = Ginh[nid][(t >meanFirstSpikeTime - window/2.)&(t < meanFirstSpikeTime + window/2.)]    
        
        #f = plt.figure()
        #ax  = f.add_subplot(511)      
        
        GburstedInh = Ginh[neuronCounter][(t >meanFirstSpikeTime - window/2.)&(t < meanFirstSpikeTime + window/2.)]    
        
        integralConductanceOfSpiked = integral(GburstedInh, dt)
        integralConductanceOfAllSpiked.append(integralConductanceOfSpiked)     
        
        
        #ax.plot(window_t, GburstedInh)
        
        # neurons that spike later:
        integralConductanceOfOther = []
        
        laterSpikedNeurons = np.where(meanFirstSomaSpikeOfNeuronsWithRoburstDendriticSpike > meanFirstSpikeTime + window)[0]
        print "ind of neurons that spiked later: ",laterSpikedNeurons
        #numToPlot = 4
        for nid_later in laterSpikedNeurons:
            #if i >= numToPlot:
               # break
            GInhLater = Ginh[nid_later][(t >meanFirstSpikeTime - window/2.)&(t < meanFirstSpikeTime + window/2.)]    
            integralConductanceOfOther.append(integral(GInhLater, dt))    
        
        integralConductanceOfAllOther.extend(integralConductanceOfOther)     
            
        plt.figure()
        plt.hist(integralConductanceOfOther, fill=False, edgecolor='r', label='spiked later')
        ylim=plt.ylim()
        plt.vlines(integralConductanceOfSpiked, 0, ylim[1])
        plt.legend()

        if neuronCounter == 20:
            break
            #ax = f.add_subplot(510 + i + 2)
            #ax.plot(window_t, GInhLater)
             
            


        #break

    plt.figure()
    plt.hist(integralConductanceOfAllSpiked, fill=False, edgecolor='b', label='spiked neurons')
    plt.hist(integralConductanceOfAllOther, fill=False, edgecolor='r', label='spiked later')
    plt.legend()
    plt.show()


        # normalize conductance by max value    
        #Gbursted_window[neuronSavedCounter] /= np.max(Gbursted_window[neuronSavedCounter])    
        
        # normalize to o mean and unit variance
        #Gbursted_window[neuronSavedCounter] = sklearn.preprocessing.scale(Gbursted_window[neuronSavedCounter], axis=0, with_mean=True, with_std=True, copy=True)    
        
        #if nid in nid_toPlot:
          #  row = neuronPlotCounter // 4
            #col = neuronPlotCounter % 4 
            #axarrInh[row, col].plot(window_t, GburstedInh_window[neuronSavedCounter])
            #axarrExc[row, col].plot(window_t, GburstedExc_window[neuronSavedCounter])
            #axarr[row, col].vlines(meanFirstSpikeTime, 0, np.max(Ginh[nid]))
            
            #if row == 3:
                #axarrInh[row, col].set_xlabel('Time (ms)')
                #axarrExc[row, col].set_xlabel('Time (ms)')
    
            #if col == 0:
               # axarrInh[row, col].set_ylabel('Ginh (mS/cm^2)')
               # axarrExc[row, col].set_ylabel('Gexc (mS/cm^2)')
    
            #neuronPlotCounter += 1
        
        #neuronSavedCounter += 1
    
        
    
if __name__ == "__main__":
    trial = 25200
    simName = "matTrans63"
    
    dataDir = "/mnt/hodgkin/eugene/results/immature/clusters/" + simName + "/"
    testDataDir = "/mnt/hodgkin/eugene/results/immature/clusters/test/" + simName + "/trial" + str(trial) + "/"
    outFigureDir = "/mnt/hodgkin/eugene/figures/chainGrowth/111818/"
    
    
    #compare_conductanceOfNetworkNeurons(dataDir, testDataDir, outFigureDir, simName, trial)
    compare_inhibitory_weights(dataDir, testDataDir, outFigureDir, simName, trial)
    #x = np.linspace(0, 1, 1000)
    #a = x*x
    #print integral(a, x[1]-x[0])