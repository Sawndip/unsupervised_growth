# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 17:10:58 2016

@author: jingroup

Script reads gaba potentials and maturation info
"""

import reading

fileMature = "/home/eugene/Output/mature.bin"

gaba_potential, maturation_triggered, mature = reading.read_maturation_info(fileMature)

mature_id = [i for i,e in enumerate(mature) if e == 1]
maturation_triggered_id = [i for i,e in enumerate(maturation_triggered) if e == 1]

print "gaba_potential = ", gaba_potential
print "Mature neurons: ", mature_id
print "Neurons with maturation triggered: ",maturation_triggered_id

immature_id = [i for i,e in enumerate(mature) if e == 0]

#print "Immature neurons: ", immature_id