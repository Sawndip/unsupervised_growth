# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 10:39:21 2017

@author: jingroup

script checks if network was saved correctly
"""

import reading

num_file = 0 # file number

filename_RAxy_before = "/home/eugene/Output/networks/networkTest/RA_xy_initial.bin"
filename_RAxy_after = "/home/eugene/Output/networks/networkTest/RA_xy_" + str(num_file) + "_.bin"

filename_RAI_before = "/home/eugene/Output/networks/networkTest/RA_I_connections_initial.bin"
filename_RAI_after = "/home/eugene/Output/networks/networkTest/RA_I_connections_" + str(num_file) + "_.bin"

filename_IRA_before = "/home/eugene/Output/networks/networkTest/I_RA_connections_initial.bin"
filename_IRA_after = "/home/eugene/Output/networks/networkTest/I_RA_connections_" + str(num_file) + "_.bin"


# check RA coordinates
xx_before, yy_before = reading.read_coordinates(filename_RAxy_before)
xx_after, yy_after = reading.read_coordinates(filename_RAxy_after)

print xx_before == xx_after
print yy_before == yy_after

# check RA -> I 
(N_RA, targets_RAI_ID_before, targets_RAI_G_before) = reading.read_connections(filename_RAI_before)
(_, targets_RAI_ID_after, targets_RAI_G_after) = reading.read_connections(filename_RAI_after)


print targets_RAI_ID_before == targets_RAI_ID_after
print targets_RAI_G_before == targets_RAI_G_after


# check I -> RA connections
(_, targets_IRA_ID_before, targets_IRA_G_before) = reading.read_connections(filename_IRA_before)
(_, targets_IRA_ID_after, targets_IRA_G_after) = reading.read_connections(filename_IRA_after)


print targets_IRA_ID_before == targets_IRA_ID_after
print targets_IRA_G_before == targets_IRA_G_after
