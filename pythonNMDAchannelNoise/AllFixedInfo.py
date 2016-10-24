# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 17:33:44 2016

@author: jingroup

Script calculates how number of input interneurons from the cell decays with distance
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 09:39:36 2016

@author: jingroup
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May 11 12:41:50 2016

@author: eugene

Script plots spatial distributions of fixed synapses
"""

filename = "/home/eugene/Output/super.net"

RA_xy = "/home/eugene/Output/RA_xy.bin"
I_xy = "/home/eugene/Output/I_xy.bin"
RA2I = "/home/eugene/Output/RA_I_connections.bin"
I2RA = "/home/eugene/Output/I_RA_connections.bin"
fileSilent = "/home/eugene/Output/silent.bin"

import reading
import math
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import random

SQUARE_SIDE = 100
SIDE = SQUARE_SIDE * math.sqrt(2)

def get_distribution_for_combined(n, synapses):
    """
    Function calculates distribution of I inputs due to training group for neurons
    with number of convergent silent connections = n
    """
    numI = []
    
    for keys, values in synapses.iteritems():
        if values[0] == n:
            numI.append(values[1])
            
    return Counter(numI)

def combine_silent_and_I_inputs_from_group(targetedSilent, targetedThroughI):
    """
    Funcrion estimates how many silent synapses and I inputs neuron reveive from training group
    """
    synapses = {}
    
    for i in xrange(len(targetedSilent)):
        for j in xrange(len(targetedSilent[i])):
            for k in xrange(len(targetedThroughI)):
                if targetedSilent[i][j] in targetedThroughI[k]:
                    synapses[targetedSilent[i][j]] = (i, k)
    return synapses

def get_number_I_inputs_from_group(group, RAtargets, Itargets):
    """
    Function calculates distribution of number of I inputs to RA neurons due 
    to RA neurons in the group
    """
    targetedI = set() # interneurons targeted by RA neurons in the group
    
    for n in group:
        targetedI.update(RAtargets[n])
    
    targetedRA = [] # RA neurons targeted by RA neurons in the group through interneurons
    
    for n in targetedI:
        targetedRA.extend(Itargets[n])
        
    numInputsFromGroup = Counter(targetedRA)
    distOfInputsFromGroup  = Counter(numInputsFromGroup.values())# distribution of number of inputs from neuron in the group
     
    maxNumOfInputs = max(numInputsFromGroup.values())

    targeted = [[] for i in  xrange(maxNumOfInputs+1)]
        
    for x,y in numInputsFromGroup.iteritems():
        targeted[y].append(x)
        
    return distOfInputsFromGroup, targeted
    


def get_average_number_I_inputs(groupSize, RAtargets, Itargets, nRepeats):
    """ 
    Function finds average number of connections pool RA neurons rerecive due to indirect 
    connections through interneurons from the group of randomly selected neurons
    Here we ignore multiple connections from neurons in the group to the same interneuron.    
    """
    sumDist = Counter() # dictionary which stores the sum of all distributions    
       
    groups = set() # groups that were already counted
    
    for j in xrange(nRepeats):
        # create random group of neuron of size groupSize
        group = set()    
        
        stop = False      
        
        while not stop:
            for i in xrange(groupSize):
                rand = random.randint(0, N_RA - 1)
                while rand in group:
                    rand = random.randint(0, N_RA - 1)
                group.add(rand)
            
            # check if group was already counted before
            if group not in groups:
                groups.add(frozenset(group))
                stop = True
        
        distOfInputsFromGroup, unused = get_number_I_inputs_from_group(group, RAtargets, Itargets)
        sumDist += distOfInputsFromGroup                
    
    #print "All groups considered: ",groups
    sumDist = {x: y/float(nRepeats) for x,y in sumDist.items()}    
    
    return sumDist        

def get_number_converging_inputs_from_group(group, RAtargets):
    """
    Function calculates distribution of number of converging inputs that interneurons
    receive due to RA neurons in the group
    """
    Itargets = []    
        
    for n in group:
        Itargets.extend(RAtargets[n])
        
    numInputsFromGroup = Counter(Itargets)
    distOfInputsFromGroup  = Counter(numInputsFromGroup.values())# distribution of number of inputs from neuron in the group
    
    return distOfInputsFromGroup

def get_average_number_converging_inputs(groupSize, RAtargets, nRepeats):
    """
    Function finds average number of inhibitory targets with converging inputs from randomly selected RA group
    of size groupSize. Average is calculated by picking random neurons nRepeats times.    
    """
    
    sumDist = Counter() # dictionary which stores the sum of all distributions    

    groups = set() # groups that were already counted
    
    for j in xrange(nRepeats):
        # create random group of neuron of size groupSize
        group = set()    
        
        stop = False      
        
        while not stop:
            for i in xrange(groupSize):
                rand = random.randint(0, N_RA - 1)
                while rand in group:
                    rand = random.randint(0, N_RA - 1)
                group.add(rand)
            
            # check if group was already counted before
            if group not in groups:
                groups.add(frozenset(group))
                stop = True
        
        
        distOfInputsFromGroup  = get_number_converging_inputs_from_group(group, RAtargets)# distribution of number of inputs from neuron in the group
        
        sumDist += distOfInputsFromGroup
        
        #print "sumDist:", sumDist    
    
    sumDist = {x: y/float(nRepeats) for x,y in sumDist.items()}    
    
    return sumDist

def get_convergent_silent_inputs_from_group(group, RAtargets):
    """
    Function return list of arrays when each list contains neurons with certain amount of convergent
    inputs from the group
    """
    
    targetedRA = [] # RA neurons targeted by neurons in the group
    
    for n in group:
        targetedRA.extend(RAtargets[n])
        
    dist = Counter(targetedRA)
    maxNumOfInputs = max(dist.values())

    targetedSilent = [[] for i in  xrange(maxNumOfInputs+1)]
        
    for x,y in dist.iteritems():
        targetedSilent[y].append(x)
        
    return targetedSilent
        
        
        

def d(x1, y1, x2, y2):
    """ Computes distance between two points in space"""
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

def std(x):
    """Computes standard deviation"""
    mean = np.mean(x)
    
    sum = 0    
    for e in x:
        sum += (e - mean)**2
        
    return math.sqrt(sum/len(x))
        



(x_I, y_I) = reading.read_coordinates(I_xy)
(x_RA, y_RA) = reading.read_coordinates(RA_xy)



#print "(x_I, y_I):", zip(x_I, y_I)
#print "(x_RA, y_RA):", zip(x_RA, y_RA)


(N_RA, RA_targets, RA_targets_G) = reading.read_connections(RA2I)
(N_I, I_targets, I_targets_G) = reading.read_connections(I2RA)
(N_RA, RA_silent_targets, RA_silent_targets_G) = reading.read_connections(fileSilent)

num_targetsRAI = [len(t) for t in RA_targets]
num_targetsIRA = [len(t) for t in I_targets]
num_silent_targets = [len(t) for t in RA_silent_targets]

targetsRAI = [t for sublist in RA_targets for t in sublist]
targetsIRA = [t for sublist in I_targets for t in sublist]
silentTargets = [t for sublist in RA_silent_targets for t in sublist]

print RA_silent_targets[3]

# for each neuron calculate how number of its indirect links via interneurons decays with distance

#print RA_targets
#print I_targets

distances = []
nIinputs = []

for j in xrange(N_RA):
    for i in xrange(N_RA):
        if i!=j:
            dist = d(x_RA[i], y_RA[i], x_RA[j], y_RA[j])
            num_I_inputs = 0
            # loop through all inhibitory targets and if a neuron is targeted by interneuron
            # increase count
            
            for Itarget in RA_targets[j]:
                #@print Itarget
                if i in I_targets[Itarget]:
                    num_I_inputs += 1
            
            distances.append(dist/SIDE)
            nIinputs.append(num_I_inputs)

# sort lists by distance
distances, nIinputs = zip(*sorted(zip(distances, nIinputs)))

#print distances

# calculate mean distribution
numBins = 400

spaceBins = np.linspace(0, 1, numBins)
averageIinputs = []
#print spaceBins
#print len(spaceBins)
ind = 0

for i in xrange(numBins - 1):

    numPoints = 0
    inputs = 0

    while (distances[ind] >= spaceBins[i]) and (distances[ind] < spaceBins[i+1]) and (ind < len(distances) -1):
        numPoints += 1
        inputs += nIinputs[ind]
        ind += 1

    if (numPoints > 0):
        averageIinputs.append(float(inputs/float(numPoints)))
    else:
        averageIinputs.append(0)

averageIinputs.append(0)
    
#print averageIinputs
#print len(averageIinputs)

# calculate distributions for each number of inputs

# sort array by nIinputs
nIinputs, distances = zip(*sorted(zip(nIinputs, distances)))

minIinputs = min(nIinputs)
maxIinputs = max(nIinputs)

# estimate number of points with corresponding number of I inputs
numOfSamples = Counter(nIinputs).values()
distDiffNumInputs = []

counter = 0
for i in range(minIinputs, maxIinputs+1):
    dist = []
    for j in xrange(numOfSamples[i - minIinputs]):
        dist.append(distances[counter])
        counter += 1
    distDiffNumInputs.append(dist)
    

# calculate number of converging inputs to interneurons
group = [0, 1, 2, 3]

IgroupTargets = []

for n in group:
    IgroupTargets.extend(RA_targets[n])

numInputsFromGroup = Counter(IgroupTargets)
#print numInputsFromGroup

repeat = 500
groupSize = 4
trainingGroup = [0, 1, 2, 3]

# calculate number of converging inputs interneurons receive from a random group 
# of training neurons
distIconvergingForTrainingGroup = get_number_converging_inputs_from_group(trainingGroup, RA_targets)
print "distribution of number of converging inputs from training group: ",distIconvergingForTrainingGroup

averageConvergingInputs = get_average_number_converging_inputs(groupSize, RA_targets, repeat)
print "Average distribution of converging inputs to interneurons: ",averageConvergingInputs

# calculate number of inputs RA neurons receive from a random group 
# of training neurons through inhibitory connections
distForTrainingGroup, targetedThroughIFromTraining = get_number_I_inputs_from_group(trainingGroup, RA_targets, I_targets)
print "distribution of number of I inputs pool neurons receive through indirect connections from training group: ",distForTrainingGroup
#print targetedThroughIFromTraining

averageIinputsFromGroup = get_average_number_I_inputs(groupSize, RA_targets, I_targets, repeat)
print "Average distribution of number of I inputs pool neurons receive through interneurons: ",averageIinputsFromGroup

# calculate convergence of silent inputs
targetedSilentFromTraining = get_convergent_silent_inputs_from_group(trainingGroup, RA_silent_targets)
#print targetedSilentFromTraining

# check how many inputs through I do targeted silent have
synapses = combine_silent_and_I_inputs_from_group(targetedSilentFromTraining, targetedThroughIFromTraining)
print synapses

# get distribution for all numbers of convergent silent inputs
distCombined = []

for i in xrange(groupSize):    
    distCombined.append(get_distribution_for_combined(i+1, synapses))
print distCombined
#print distances
#print nIinputs
# Here we calculate mean distance between all RA and I neurons

all_dist_RAandI = []

for i in xrange(N_RA):
    for j in xrange(N_I):
        dist = d(x_RA[i], y_RA[i], x_I[j], y_I[j])
        all_dist_RAandI.append(dist)

mean_distance_RAandI = np.mean(all_dist_RAandI)



# Here we calculate mean distance between all RA and RA neurons
        
all_dist_RAandRA = []


for i in xrange(N_RA):
    for j in xrange(i, N_RA):
        if j!=i:
            dist = d(x_RA[i], y_RA[i], y_RA[j], y_RA[j])
        all_dist_RAandRA.append(dist)

mean_distance_RAandRA = np.mean(all_dist_RAandRA)



# Here we calculate mean distance between RA and their I targets

dist_RA2I = []

for i in xrange(N_RA):
    if len(RA_targets[i]) > 1:
        for target_ind in RA_targets[i]:
            dist_RA2I.append(d(x_RA[i], y_RA[i], x_I[target_ind], y_I[target_ind])/SIDE)


# Here we calculate mean distance between I and their RA targets

dist_I2RA = []

for i in xrange(N_I):
    if len(I_targets[i]) > 1:
        for target_ind in I_targets[i]:
            dist_I2RA.append(d(x_I[i], y_I[i], x_RA[target_ind], y_RA[target_ind])/SIDE)
            

# let's count number of input I connections

numInputI = Counter(targetsIRA).values()
numWithIinputs = len(numInputI)
numInputI.extend([0]*(N_RA - numWithIinputs))

# let's count number of input RA connections
numInputRA = Counter(targetsRAI).values()
numWithRAinputs = len(numInputRA)
numInputRA.extend([0]*(N_I - numWithRAinputs))

#print dist_RA2I
#print "test:", d(17,31,33,45)
print
print "Mean distance between all RA and I neurons normalized by side of square area: ", mean_distance_RAandI/SIDE
print "Standard deviation for all distances between RA and I: ", std(all_dist_RAandI)/SIDE
print "Mean distance between all RA and RA neurons normalized by side of square area: ", mean_distance_RAandRA/SIDE
print "Standard deviation for all distances between RA and RA: ", std(all_dist_RAandRA)/SIDE

print "Mean distance for RA to I connections normalized by side of square area: ", np.mean(dist_RA2I)
print "Standard deviation for distances from RA to I: ", std(dist_RA2I)
print "Mean distance for I to RA connections normalized by side of square area: ", np.mean(dist_I2RA)
print "Standard deviation for distances from I to RA: ", std(dist_I2RA)

# plot spatial distribution of RA-I and I-RA connections
numBins = 30

f = plt.figure()
ax = f.add_subplot(211)
ax.hist(dist_RA2I, numBins)
ax.set_xlabel("Normalized distance")
ax.set_ylabel("Number of connections")
ax.set_title("Spatial distribution of connections from RA onto I neurons")
ax.set_xlim([0,1])

numBins = 25

ax = f.add_subplot(212)
ax.hist(dist_I2RA, numBins)
ax.set_xlabel("Normalized distance")
ax.set_ylabel("Number of connections")
ax.set_title("Spatial distribution of connections from I onto RA neurons")
ax.set_xlim([0,1])

# plot distribution of in and out connections for I and RA neurons
numBins = 25

f = plt.figure()
ax = f.add_subplot(221)
ax.hist(num_targetsIRA, numBins)
ax.set_xlabel("# of output")
ax.set_title("Distribution of outputs from I")

ax = f.add_subplot(222)
ax.hist(numInputRA, numBins)
ax.set_xlabel("# of inputs")
ax.set_title("Distribution of inputs to I")

ax = f.add_subplot(223)
ax.hist(num_targetsRAI, numBins)
ax.set_xlabel("# of output")
ax.set_title("Distribution of outputs from RA")

ax = f.add_subplot(224)
ax.hist(numInputI, numBins)
ax.set_xlabel("# of inputs")
ax.set_title("Distribution of inputs to RA")



# plot spatial distributions of I inputs due to central RA neuron
f = plt.figure()
ax = f.add_subplot(211)
ax.scatter(distances, nIinputs)
#ax.hist(numInputRA, numBins)
ax.set_xlim([0, max(distances)])
ax.set_ylim([0, max(nIinputs) +1])
ax.set_ylabel("# of inhibitory inputs due to central RA neuron")
ax.set_xlabel("distance from RA neuron")
ax.set_title("Spatial distribution of number of inhibitory inputs to RA neuron due to central RA neuron")

ax = f.add_subplot(212)
ax.scatter(spaceBins, averageIinputs)
#ax.hist(numInputRA, numBins)
ax.set_xlim([0, max(spaceBins)])
ax.set_ylim([0, max(averageIinputs) +1])
ax.set_ylabel("# of inhibitory inputs due to central RA neuron")
ax.set_xlabel("distance from RA neuron")
ax.set_title("Average spatial distribution of number of inhibitory inputs to RA neuron due to central RA neuron")


# plot spatial distribution for ach number of I inputs
numBins = 50
numPlots = maxIinputs - minIinputs + 1
nrows = int(math.sqrt(numPlots))
ncol = numPlots/nrows

#print minIinputs
#print maxIinputs

if numPlots % nrows != 0:
    ncol = ncol + 1

#print nrows
#print ncol

f, axarr = plt.subplots(nrows, ncol)

for i in xrange(len(distDiffNumInputs)):
    
    #ax.scatter(spaceBins, averageIinputs)
    axarr[i/ncol, i%ncol].hist(distDiffNumInputs[i], numBins, normed=True)
    axarr[i/ncol, i%ncol].set_xlim([0, 1])
    #ax.set_ylim([0, max(averageIinputs) +1])
    axarr[i/ncol, i%ncol].set_ylabel("pdf")
    axarr[i/ncol, i%ncol].set_xlabel("distance from RA neuron")
    axarr[i/ncol, i%ncol].set_title("{0} I inputs".format(minIinputs+i))

# plot distribution of inputs from training group
width = 0.25

f = plt.figure()
ax1 = f.add_subplot(211) # converging inputs to interneurons
ax1.bar(averageConvergingInputs.keys(), averageConvergingInputs.values(), width)
#ax.hist(numInputRA, numBins)
ax1.set_xticks(averageConvergingInputs.keys())
ax1.set_xlim([0, max(averageConvergingInputs.keys())+1])
ax1.set_ylim([0, max(averageConvergingInputs.values()) +1])
ax1.set_ylabel("# of occurences")
ax1.set_xlabel("# of converging inputs")
ax1.set_title("Average distribution of number of inputs to interneurons from group")

ax2 = f.add_subplot(212)
ax2.bar(averageIinputsFromGroup.keys(), averageIinputsFromGroup.values(), width)
#ax.hist(numInputRA, numBins)
ax2.set_xticks(averageIinputsFromGroup.keys())
ax2.set_xlim([0, max(averageIinputsFromGroup.keys())+1])
ax2.set_ylim([0, max(averageIinputsFromGroup.values()) +1])
ax2.set_ylabel("# of occurences")
ax2.set_xlabel("# of I inputs from group")
ax2.set_title("Average distribution of number of inputs to RA neurons through interneurons from group")


# plot distribution of combined distribution
width = 0.25
numBins = 50
numPlots = groupSize
nrows = int(math.sqrt(numPlots))
ncol = numPlots/nrows

#print minIinputs
#print maxIinputs

if numPlots % nrows != 0:
    ncol = ncol + 1

#print nrows
#print ncol

f, axarr = plt.subplots(nrows, ncol)

for i in xrange(len(distCombined)):
    
    #ax.scatter(spaceBins, averageIinputs)
    axarr[i/ncol, i%ncol].bar(distCombined[i].keys(), distCombined[i].values(), width)
    axarr[i/ncol, i%ncol].set_xlim([0, max(distCombined[i].keys())+1])
    axarr[i/ncol, i%ncol].set_ylim([0, max(distCombined[i].values()) +1])
    axarr[i/ncol, i%ncol].set_ylabel("# of occurences")
    axarr[i/ncol, i%ncol].set_xlabel("# of I inputs")
    #ax.set_ylim([0, max(averageIinputs) +1])
    axarr[i/ncol, i%ncol].set_title("{0} convergent silent inputs".format(i+1))


#numBins = 25

#f = plt.figure()
#ax = f.add_subplot(111)
#ax.hist(num_silent_targets, numBins)
#ax.set_xlabel("# of output")
#ax.set_title("Distribution of silent outputs from RA")



plt.show()

