# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 09:24:03 2018

@author: jingroup
"""
import os
import reading

dirname = "/mnt/hodgkin/eugene/results/immature/clusters/matTrans69"
trial_number = 27800
fileActivityHistory = os.path.join(dirname, "activity_history_" + str(trial_number) + ".bin")

(_, _, activity_history) = reading.read_activity_history(fileActivityHistory) 


print list(activity_history[55])