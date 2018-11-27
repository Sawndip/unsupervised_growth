# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 10:14:01 2018

@author: jingroup

Script analyzes burst tighness for syllable locked neurons
"""
import reading
import os
import utils
import numpy as np
import matplotlib.pyplot as plt

def plotSpikes(spikes, axis):
    """
    Function makes a raster plot of spikes
    
    Input: list of list with spikes
           axis - axis where to plot spikes
    """
    if len(spikes) == 0:
        print "Spike array is empty!"
        return
    
    mint = 0.0
    maxt = -1.0    
    
    for i in range(len(spikes)):
        for spike in spikes[i]:
            axis.plot([spike, spike], [i-0.5, i+0.5],c='k')
            if spike > maxt:
                maxt = spike
            
    axis.set_xlim([mint-1.0, maxt+1.0])


def loadLockingInfo(filename):
    """
    Function reads average conductances from a file
    """
    locking = np.load(filename)
    return locking['syllableLockingTIme'], locking['pvalues'], locking['meanFsi'], locking['medianFsi'], locking['meanIsi'], locking['medianIsi']


def calculateAverageFiringRate(spikeTrials, tmin, tmax, binsize):
    """
    Function calculates an average firing rate for several trials. 
    
    Input: spikeTrials - array with spikes in different trials
           tmin - left time boundary
           tmax - right time boundary
           binsize - firing rate resolution
           
    Output: binCenters, smoothed average firing rate
    """
    
    binCenters = [tmin + i*binsize for i in range(int((tmax-tmin)/binsize))]
    firingRate = np.zeros(len(binCenters), np.float32)

    for spikes in spikeTrials:
        for spike in spikes:
            if spike >= tmin and spike <= tmax:
                binInd = int((spike - binCenters[0] - binsize/2.)/binsize)
                #print spike, binInd
                firingRate[binInd] += 1.0
        
    averageFiringRate = firingRate / float(len(spikeTrials))
    
    return binCenters, averageFiringRate

def moving_average(a, n=3):
    """
    Function smoothes array over n points
    """
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return np.r_[a[0:(n-1)/2],ret[n - 1:]/n,a[-(n-1)/2:]] 

def shiftSpikes(spikes, offset, leftT, rightT):
    """
    Shift spikes by adding random offset so that spikes remain between
    left and right boundaries
    """
    if np.random.uniform() >= 0.5:
        sign = 1.0
    else:
        sign = -1.0
    
    shiftedSpikes = np.array(spikes) + sign * np.random.uniform() * offset
    
    #print "lT, rT = ",leftT,rightT
    #print "Spikes: ",spikes
    #print "Shifted spikes: ",shiftedSpikes
    #print "shift = ",sign*offset
    
    smallerLeft = np.where(shiftedSpikes < leftT)
    largerRight = np.where(shiftedSpikes > rightT)
    
    shiftedSpikes[smallerLeft] = rightT - (leftT - shiftedSpikes[smallerLeft])
    shiftedSpikes[largerRight] = leftT + (shiftedSpikes[largerRight] - rightT)
    
    #print "Shifted spikes adjusted: ",shiftedSpikes
    
    return shiftedSpikes

def computeLockingInfo(dataDir, testDir):
    """
    Compute syllable locking significance and spike pattern properties for HVC-RA neurons from testing trials
    """
    N_RA, _ = reading.read_num_neurons(os.path.join(dataDir, "num_neurons.bin"))
    training_neurons = reading.read_training_neurons(os.path.join(dataDir, "training_neurons.bin"))
    
    
    files = os.listdir(testDir)
    
    somatic_spikes = [[] for i in range(N_RA)]
    
    for f in files:
        if "test_spike_times_soma" in f:
            spike_times_s,  neuron_id_s,  ordered_spikes_s, neuron_ordered_id_s = utils.getSpikes(os.path.join(testDir, f))
            
            
            # find earliest first spike time of training neurons to allign spikes to syllable onsets
            minTrainingFirstSpike = 1e6 # earliest first spike time of training neuron
            
            for spikes, neuronId in zip(spike_times_s, neuron_id_s):            
                if neuronId[0] in training_neurons:
                    if spikes[0] < minTrainingFirstSpike:
                        minTrainingFirstSpike = spikes[0]
    
            #print "earliest first spike time of training neurons: ",minTrainingFirstSpike
            print minTrainingFirstSpike
            for spikes, neuronId in zip(spike_times_s, neuron_id_s):        
                somatic_spikes[neuronId[0]].append([spike - minTrainingFirstSpike for spike in spikes])
                
    #print somatic_spikes
    
    earliestLocking = -500 # earliest locking to syllable onset in ms
    latestLocking = 500 # latest locking to syllable onset in ms
    nbootstrap = 1000 # number of bootstrap samples
    interBurstGap = 30.0 # gap between bursts in ms
    
    if not os.path.isdir(os.path.join(testDir, 'figures')):
        os.mkdir(os.path.join(testDir, 'figures'))       
    
    plt.ioff()
    
    
    pvalues = []
    
    meanFsi = []
    medianFsi = []
    meanIsi = []
    medianIsi = []
    
    syllableLockingTIme = []
    
    for nid in range(N_RA):
        if len(somatic_spikes[nid]) == 0:
            pvalues.append(1.0)
            syllableLockingTIme.append(np.nan)
            meanFsi.append(np.nan)
            medianFsi.append(np.nan)
            meanIsi.append(np.nan)
            medianIsi.append(np.nan)
    
            continue
        
        #print somatic_spikes[nid]
        binCenters, averageFiringRate = calculateAverageFiringRate(somatic_spikes[nid], earliestLocking, latestLocking, 1.0)    
    
        smoothWindow = 21    
        
        smoothedFiringRate = moving_average(averageFiringRate, smoothWindow)  
        
       
           
        # estimate significance of peak in firing rate
        peakFiringRate = np.max(smoothedFiringRate[smoothWindow:-smoothWindow])
        peakFiringRateTime = np.array(binCenters)[smoothWindow:-smoothWindow][np.argmax(smoothedFiringRate[smoothWindow:-smoothWindow])]
          
        
        bootstrapOffset = 500.0 # offset of spikes in bootstrap samples in ms    
        
        bootstrapPeakFiringRates = np.empty(nbootstrap, np.float32)
        
        for i in range(nbootstrap):
            spikesShifted = []
            
            for spikes in somatic_spikes[nid]:
                spikesShifted.append(shiftSpikes(spikes, bootstrapOffset, earliestLocking, latestLocking))
        
            binCentersBootstrap, averageFiringRateBootstrap = calculateAverageFiringRate(spikesShifted, earliestLocking, latestLocking, 1.0)    
        
            smoothedFiringRateBootstrap = moving_average(averageFiringRateBootstrap, smoothWindow)  
            
            bootstrapPeakFiringRates[i] = np.max(smoothedFiringRateBootstrap[smoothWindow:-smoothWindow])
        
        isi = []
        fsi = []        
        
        for spikes in somatic_spikes[nid]:
             #print spikes
             if len(spikes) > 1:
                 bursts = utils.getBurstsForNeuron(spikes, interBurstGap)
                 #print bursts
                 for burst in bursts:
                     if len(burst) > 1:
                         isi.extend(np.diff(burst))
                         fsi.append(burst[1]-burst[0])
           

        print len(fsi)
        print len(somatic_spikes[nid])              
        #print isi
        #print fsi
                        
            # plot example of bootstrap spikes
    #==============================================================================
    #         if i == 0:
    #             f2 = plt.figure()
    #             plt.suptitle("Bootstrap example")
    #             ax1_2 = f2.add_subplot(211)
    #             
    #             plotSpikes(spikesShifted, ax1_2)
    #             
    #             ax1_2.set_xlim([earliestLocking - 100, latestLocking + 100])
    #             
    #             ax2_2 = f2.add_subplot(212)
    #             ax2_2.step(binCenters, averageFiringRate, where='mid')
    #             ax2_2.step(binCenters, smoothedFiringRate, where='mid', linewidth=3.0)
    #             ax2_2.set_xlim([earliestLocking - 100, latestLocking + 100])
    #             ax2_2.set_xlabel('Time relative to syllable onset (ms)')
    #             ax2_2.set_ylabel('Firing rate (1/ms)')
    #             
    #             plt.savefig(neuronDir + "/" + 'bootstrapEx.png', bbox_inches='tight')
    #             plt.close(f2)
    #               
    #==============================================================================
        
        pvalue = float(len(np.where(bootstrapPeakFiringRates >= peakFiringRate)[0])) / float(nbootstrap)
        print "p-value = ",pvalue
        print "Syllable locking time = ",peakFiringRateTime
        
        f1 = plt.figure()
        plt.suptitle('Neuron {0} with p = {1}'.format(nid, pvalue))
        ax1 = f1.add_subplot(211)
        
        plotSpikes(somatic_spikes[nid], ax1)
        
        ax1.set_xlim([earliestLocking - 100, latestLocking + 100])
        
        ax2 = f1.add_subplot(212)
        ax2.step(binCenters, averageFiringRate, where='mid', linewidth=3.0, zorder=2, color='b', label='firing rate')
    
        ax2.step(binCenters, smoothedFiringRate, where='mid', linewidth=3.0, zorder=2, color='r', label='smooth firing rate')
        ax2.set_xlim([earliestLocking - 100, latestLocking + 100])
        
        #ax1.set_xlim([-100, 300])     
        #ax2.set_xlim([-100, 300])
        
        ax2.set_xlabel('Time relative to syllable onset (ms)')
        ax2.set_ylabel('Firing rate (1/ms)')
        plt.legend()
        
         
        plt.savefig(os.path.join(testDir, 'figures/firingRate{0}.png'.format(nid)), bbox_inches='tight')
        plt.close(f1)
        
        pvalues.append(pvalue)
        syllableLockingTIme.append(peakFiringRateTime)
    
        if len(fsi) > 0:
            meanFsi.append(np.mean(fsi))
            medianFsi.append(np.median(fsi))
            meanIsi.append(np.mean(isi))
            medianIsi.append(np.median(isi))
        else:
            meanFsi.append(np.nan)
            medianFsi.append(np.nan)
            meanIsi.append(np.nan)
            medianIsi.append(np.nan)

        np.savez(os.path.join(testDir, "lockingInfo.npz"), syllableLockingTIme=syllableLockingTIme, pvalues=pvalues, 
                 meanFsi=meanFsi, medianFsi=medianFsi, meanIsi=meanIsi, medianIsi=medianIsi)
                 
        
        #if nid == 5:
          #  break
    
if __name__ == "__main__":

    trial = 23800
    simName = "matTrans62"
    dataDir = "/mnt/hodgkin/eugene/results/immature/clusters/" + simName
    testDir = "/mnt/hodgkin/eugene/results/immature/clusters/test/" + simName + "/trial" + str(trial)
    
    #computeLockingInfo(dataDir, testDir)
    # read from file
    training_neurons = reading.read_training_neurons(os.path.join(dataDir, "training_neurons.bin"))

    print "Training neurons:", training_neurons

    syllableLockingTIme, pvalues, meanFsi, medianFsi, meanIsi, medianIsi = loadLockingInfo(os.path.join(testDir, "lockingInfo.npz"))

    print "nan p values:", np.any(np.isnan(pvalues))
    print "nan syllable locking time:", np.any(np.isnan(syllableLockingTIme))
    
    print "Locked neurons with negative syllable locking time: ",np.where((syllableLockingTIme < 0)&(pvalues < 0.05))
    #print pvalues[syllableLockingTIme < 0][np.where(pvalues[syllableLockingTIme < 0] <0.05)]
    
    print "Locked neurons with negative syllable locking time and non-nan median fsi:", np.where((~np.isnan(syllableLockingTIme))&(syllableLockingTIme < 0)&(pvalues < 0.05)&(~np.isnan(medianFsi)))
    
    ind_neurons = np.where((~np.isnan(syllableLockingTIme))&(syllableLockingTIme < 0)&(pvalues < 0.05)&(~np.isnan(medianFsi)))
    print syllableLockingTIme[ind_neurons]
    print medianFsi[ind_neurons]
    print pvalues[ind_neurons]
    #print "Locked neurons with negative syllable locking time and non-nan median fsi:", np.where((~np.isnan(syllableLockingTIme))&(syllableLockingTIme < 0)&(pvalues < 0.05)&(~np.isnan(medianFsi)))

    print syllableLockingTIme
    print pvalues
    print meanFsi
    print medianFsi
    print meanIsi
    print medianIsi
    
 
    plt.figure()
    plt.scatter(syllableLockingTIme[(pvalues < 0.05) & (~np.isnan(medianFsi))], medianFsi[(pvalues < 0.05) & (~np.isnan(medianFsi))])
    
    plt.show()
