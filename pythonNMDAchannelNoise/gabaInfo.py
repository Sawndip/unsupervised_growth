# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 17:10:58 2016

@author: jingroup

Script reads gaba potentials and maturation info
"""

import reading

fileGaba = "/home/eugene/Output/gaba_potential.bin"
fileMature = "/home/eugene/Output/mature.bin"

gaba_potential = reading.read_gaba_potential(fileGaba)
mature = reading.read_mature(fileMature)

mature_id = [i for i,e in enumerate(mature) if e == 1]

print "gaba_potential = ", gaba_potential
print "Mature neurons: ", mature_id


immature_id = [i for i,e in enumerate(mature) if e == 0]

#print "Immature neurons: ", immature_id