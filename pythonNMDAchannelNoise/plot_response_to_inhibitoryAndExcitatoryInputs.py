# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 13:27:30 2018

@author: jingroup

Script plots neuron responses to inhibitory and excitatory inputs at different
ages
"""
import matplotlib.pyplot as plt

WEIGHT_UNITS = 10 # pS

################### Age 0 Different inhibitory kicks #########################

weights_kick0_age0 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]
p_fire_soma_kick0_age0 = [0, 0, 0, 0, 0, 0.09, 0.26, 0.72, 0.91, 0.99, 1, 1, 1, 1, 1]
p_fire_dend_kick0_age0 = [0, 0, 0, 0, 0, 0.09, 0.26, 0.72, 0.91, 0.99, 1, 1, 1, 1, 1]

weights_kick1_age0 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]
p_fire_soma_kick1_age0 = [0, 0, 0, 0.014, 0.076, 0.368, 0.686, 0.914, 0.99, 1, 1, 1, 1, 1, 1]
p_fire_dend_kick1_age0 = [0, 0, 0, 0.014, 0.076, 0.368, 0.686, 0.914, 0.99, 1, 1, 1, 1, 1, 1]

weights_kick2_age0 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]
p_fire_soma_kick2_age0 = [0, 0.002, 0.028, 0.118, 0.312, 0.66, 0.886, 0.984, 1, 1, 1, 1, 1, 1, 1]
p_fire_dend_kick2_age0 = [0, 0.002, 0.028, 0.118, 0.312, 0.66, 0.886, 0.984, 1, 1, 1, 1, 1, 1, 1]


weights_kick3_age0 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]
p_fire_soma_kick3_age0 = [0.008, 0.026, 0.096, 0.282, 0.512, 0.81, 0.952, 0.99, 1, 1, 1, 1, 1, 1, 1]
p_fire_dend_kick3_age0 = [0.008, 0.026, 0.096, 0.282, 0.512, 0.81, 0.952, 0.99, 1, 1, 1, 1, 1, 1, 1]


f = plt.figure()

plt.plot([w*WEIGHT_UNITS for w in weights_kick0_age0], p_fire_soma_kick0_age0, label="G_kick = 0.0")
plt.plot([w*WEIGHT_UNITS for w in weights_kick1_age0], p_fire_soma_kick1_age0, label="G_kick = 1.0")
plt.plot([w*WEIGHT_UNITS for w in weights_kick2_age0], p_fire_soma_kick2_age0, label="G_kick = 2.0")
plt.plot([w*WEIGHT_UNITS for w in weights_kick3_age0], p_fire_soma_kick3_age0, label="G_kick = 3.0")

plt.xlabel('Excitatory weight [pS]')
plt.ylabel('probability to fire')
plt.legend()

################### Kick 0 Different neuron ages #########################

weights_kick0_age1 = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950]
p_fire_soma_kick0_age1 = [0, 0, 0, 0, 0, 0, 0, 0, 0.824, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
p_fire_dend_kick0_age1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.77, 0.996, 1, 1, 1, 1, 1, 1, 1, 1]

weights_kick0_age05 = [0, 35, 70, 105, 140, 175, 210, 245, 280, 315, 350, 385, 420, 455, 490]
p_fire_soma_kick0_age05 = [0, 0, 0, 0, 0, 0, 0, 0.598, 1, 1, 1, 1, 1, 1, 1]
p_fire_dend_kick0_age05 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.004, 0.482, 0.89, 1, 1, 1]

weights_kick0_age025 = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240, 255, 270, 285, 300, 315, 330, 345, 360]
p_fire_soma_kick0_age025 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.032, 0.296, 0.826, 0.996, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
p_fire_dend_kick0_age025 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.032, 0.296, 0.826, 0.996, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]


weights_kick0_age100 = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950]
p_fire_soma_kick0_age100 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.714, 1, 1, 1, 1, 1]
p_fire_dend_kick0_age100 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.104, 0.952, 1, 1, 1]

f = plt.figure()

plt.plot([w*WEIGHT_UNITS for w in weights_kick0_age0], p_fire_soma_kick0_age0, label="age = 0.0")
plt.plot([w*WEIGHT_UNITS for w in weights_kick0_age025], p_fire_soma_kick0_age025, label="age = 0.25")
plt.plot([w*WEIGHT_UNITS for w in weights_kick0_age05], p_fire_soma_kick0_age05, label="age = 0.5")
plt.plot([w*WEIGHT_UNITS for w in weights_kick0_age1], p_fire_soma_kick0_age1, label="age = 1.0")
plt.plot([w*WEIGHT_UNITS for w in weights_kick0_age100], p_fire_soma_kick0_age100, label="age = 100.0")

plt.xlabel('Excitatory weight [pS]')
plt.ylabel('probability to fire')
plt.legend()


################### Kick 3 Different neuron ages #########################
weights_kick3_age100 = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000]
p_fire_soma_kick3_age100 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]
p_fire_dend_kick3_age100 = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1]

weights_kick3_age1 = [0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250, 3500]
p_fire_soma_kick3_age1 = [0, 0, 0, 0, 0, 0, 0, 0, 0.012, 1, 1, 1, 1, 1, 1]
p_fire_dend_kick3_age1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.254, 1, 1, 1, 1, 1]

weights_kick3_age05 = [0, 150, 300, 450, 600, 750, 900, 1050, 1200, 1350, 1500, 1650, 1800, 1950, 2100]
p_fire_soma_kick3_age05 = [0, 0, 0, 0, 0, 0, 0, 0.996, 1, 1, 1, 1, 1, 1, 1]
p_fire_dend_kick3_age05 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.138, 0.998, 1, 1]

f = plt.figure()

plt.plot([w*WEIGHT_UNITS for w in weights_kick3_age0], p_fire_soma_kick3_age0, label="age = 0.0")
#plt.plot([w*WEIGHT_UNITS for w in weights_kick0_age025], p_fire_soma_kick0_age025, label="age = 0.25")
plt.plot([w*WEIGHT_UNITS for w in weights_kick3_age05], p_fire_soma_kick3_age05, label="age = 0.5")
plt.plot([w*WEIGHT_UNITS for w in weights_kick3_age1], p_fire_soma_kick3_age1, label="age = 1.0")
plt.plot([w*WEIGHT_UNITS for w in weights_kick3_age100], p_fire_soma_kick3_age100, label="age = 100.0")

plt.xlabel('Excitatory weight [pS]')
plt.ylabel('probability to fire')
plt.legend()
plt.show()