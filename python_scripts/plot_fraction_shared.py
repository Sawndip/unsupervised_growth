#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 21 13:44:42 2018

@author: jingroup

Script plots mean fraction of shared inputs for different
network topologies
"""
import matplotlib.pyplot as plt

windows_1 = [0.1, 0.5, 1.0, 2.0, 5.0]
#mean_fraction_shared_1 = [0.73850846244604818, 0.23041884868946555, 0.057017380766853178, 0.0070876728742229984, 0.0]
# std/mean = 0.49; 
# from multiple trials
mean_fraction_shared_1 = [0.73401441405783641, 0.24045572410735583, 0.057790439780326608, 0.0051146477828070839, 0.0]
mean_fraction_1 = [0.86412690138529946, 0.57543516241722403, 0.42670363988249238, 0.33323610366334944, 0.24340463821969377]
label_1 = "matTrans31; <std/mean> = 0.40"

windows_2 = [0.1, 0.5, 1.0, 2.0, 5.0]
#mean_fraction_shared_2 = [0.77662995383308664, 0.32248550317831109, 0.12803461764400112, 0.021615020323320692, 0.0]
#std/mean = 0.66;

# from multiple trials
mean_fraction_shared_2 = [0.83860862422555393, 0.32060199350666724, 0.13493493298607456, 0.024803040548539674, 0.0]
mean_fraction_2 = [0.91891816754417133, 0.64407537028018802, 0.52594930686952135, 0.41964751123635563, 0.31017871219916326]
label_2 = "matTrans36; <std/mean> = 0.57"

windows_3 = [0.1, 0.5, 1.0, 2.0, 5.0]
#mean_fraction_shared_3 = [0.8747257848360791, 0.53657856629196099, 0.32880152250940808, 0.12844324773481575, 0.0038643816859289339]
# std/mean = 0.75; 
# from multiple trials
mean_fraction_shared_3 = [0.88529596564847424, 0.57406941632549446, 0.33752653509016195, 0.11300436761987483, 0.0035876361963318483]
mean_fraction_3 = [0.94237865532199439, 0.79080646089515716, 0.6767502741287732, 0.55710180744242421, 0.43094816732609403]
label_3 = "matTrans35; <std/mean> = 0.56"


windows_4 = [0.1, 0.5, 1.0, 2.0, 5.0]
#mean_fraction_shared_4 = [1.0, 0.99831536388140174, 0.99194820509679005, 0.98109526087120436, 0.95482583454281578]
# std/mean = 1.55; 
# from multiple trials
mean_fraction_shared_4 = [1.0, 0.99831536388140174, 0.99288123433170605, 0.98033069524814809, 0.95276215529753272]
mean_fraction_4 = [1.0, 0.99900606469002695, 0.99665161741896913, 0.99255554374215926, 0.97485548060430138]
label_4 = "matTrans40; <std/mean> = 1.36"


windows_5 = [0.1, 0.5, 1.0, 2.0, 5.0]

# from multiple trials
mean_fraction_shared_5 = [0.91846901812179593, 0.67728989243831295, 0.45682772620216383, 0.20021464368375932, 0.0024449197985059699]
mean_fraction_5 = [0.95937488651609026, 0.84762109359791138, 0.75239959162154402, 0.63444876481690593, 0.42533899272903397]
label_5 = "matTrans44; <std/mean> = 0.40"

plt.figure()

#plt.plot(windows_1, mean_fraction_shared_1, '-o', label=label_1)
#plt.plot(windows_2, mean_fraction_shared_2, '-o', label=label_2)
#plt.plot(windows_3, mean_fraction_shared_3, '-o', label=label_3)
#plt.plot(windows_4, mean_fraction_shared_4, '-o', label=label_4)
#plt.plot(windows_5, mean_fraction_shared_5, '-o', label=label_5)

plt.plot(windows_1, mean_fraction_1, '-o', label=label_1)
plt.plot(windows_2, mean_fraction_2, '-o', label=label_2)
plt.plot(windows_3, mean_fraction_3, '-o', label=label_3)
plt.plot(windows_4, mean_fraction_4, '-o', label=label_4)
plt.plot(windows_5, mean_fraction_5, '-o', label=label_5)


plt.xlabel('Window size (ms)')
plt.ylabel('Mean fraction of inputs to synchronous')

plt.legend()


#==============================================================================
# # comparison of different pruned chains. No delay
# windows = [0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
# 
# mean_fraction_shared_1 = [1.0, 1.0, 1.0, 1.0, 0.98571428571428577, 0.97959183673469385, 0.97959183673469385, 0.97959183673469385, 0.97959183673469385, 0.96122448979591835]
# mean_fraction_1 = [1.0, 1.0, 1.0, 1.0, 0.99761904761904774, 0.99659863945578242, 0.99659863945578242, 0.99659863945578242, 0.99594453165881747, 0.98565821933168873]
# label_1 = "p=1.0"
# label_shared_1 = "p=1.0 shared"
# 
# mean_fraction_shared_2 = [0.83191124068675093, 0.4326911240686751, 0.22744897959183674, 0.10768545513443473, 0.079608843537414958, 0.064489795918367343, 0.058367346938775509, 0.051836734693877548, 0.034693877551020408, 0.0051020408163265302]
# mean_fraction_2 = [0.92477459237663318, 0.80183632410673233, 0.77471317590195132, 0.75834914990527236, 0.75146952065319406, 0.7417761863886998, 0.72732395608926215, 0.70400377880162812, 0.59894407773639224, 0.41933050711031589]
# label_2 = "p=0.75"
# label_shared_2 = "p=0.75 shared"
# 
# mean_fraction_shared_3 = [0.81236935789404929, 0.40096593507087336, 0.20743108628911097, 0.04478819648572735, 0.010210497093213143, 0.0022633744855967081, 0.00041152263374485601, 0.0, 0.0, 0.0]
# mean_fraction_3 = [0.90576553552376604, 0.71111357987386792, 0.61935528710734467, 0.53710062410753179, 0.50855156764575027, 0.48082218269857258, 0.44557566511592023, 0.38464629678220152, 0.30832605994778345, 0.24280955638017424]
# label_3 = "p=0.5"
# label_shared_3 = "p=0.5 shared"
# 
# mean_fraction_shared_4 = [0.87208767627222572, 0.41286957558416781, 0.17448906601267117, 0.033096600585870974, 0.0042918454935622317, 0.0042918454935622317, 0.0010729613733905579, 0.0, 0.0, 0.0]
# mean_fraction_4 = [0.93523911710606988, 0.68553034947884717, 0.52135591558338334, 0.36724097465935296, 0.29330361623950679, 0.25145331976854957, 0.21704466135535525, 0.18630186990462327, 0.16266403442600963, 0.14112016549531545]
# label_4 = "p=0.25"
# label_shared_4 = "p=0.25 shared"
# 
# plt.plot(windows, mean_fraction_1, '-ro', label=label_1)
# plt.plot(windows, mean_fraction_2, '-go', label=label_2)
# plt.plot(windows, mean_fraction_3, '-bo', label=label_3)
# plt.plot(windows, mean_fraction_4, '-mo', label=label_4)
# 
# plt.plot(windows, mean_fraction_shared_1, '-rx', label=label_shared_1)
# plt.plot(windows, mean_fraction_shared_2, '-gx', label=label_shared_2)
# plt.plot(windows, mean_fraction_shared_3, '-bx', label=label_shared_3)
# plt.plot(windows, mean_fraction_shared_4, '-mx', label=label_shared_4)
# 
# plt.legend()
# plt.xlabel('Window size (ms)')
# plt.ylabel('Mean fraction')
#==============================================================================

#==============================================================================
# # comparison of same pruned chain with different delays
# windows = [0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
# 
# mean_fraction_shared_1 = [0.81236935789404929, 0.40096593507087336, 0.20743108628911097, 0.04478819648572735, 0.010210497093213143, 0.0022633744855967081, 0.00041152263374485601, 0.0, 0.0, 0.0]
# mean_fraction_1 = [0.90576553552376604, 0.71111357987386792, 0.61935528710734467, 0.53710062410753179, 0.50855156764575027, 0.48082218269857258, 0.44557566511592023, 0.38464629678220152, 0.30832605994778345, 0.24280955638017424]
# label_1 = "delay=0.0"
# label_shared_1 = "delay=0.0 shared"
# 
# mean_fraction_shared_2 = [0.81401192920700105, 0.40244450963136802, 0.20340031289723284, 0.063871777321469314, 0.027783481633584302, 0.001642710472279261, 0.00061601642710472284, 0.0, 0.0, 0.0]
# mean_fraction_2 = [0.90752501548189435, 0.71348355338858416, 0.62071307333627657, 0.55721432457572084, 0.53561879057515605, 0.50728004871749821, 0.49356413479136724, 0.47565717755882486, 0.44132836391827845, 0.38272474981278731]
# label_2 = "delay=1.0ms"
# label_shared_2 = "delay=1.0ms shared"
# 
# mean_fraction_shared_3 = [0.81186735623211048, 0.4009782396565183, 0.20926440931563883, 0.066915495706479311, 0.031528916211293265, 0.013627049180327869, 0.0018442622950819673, 0.0014344262295081969, 0.0, 0.0]
# mean_fraction_3 = [0.90682775500910739, 0.71049065129239308, 0.62346885866865387, 0.55718480299679074, 0.53645315979867958, 0.52192201366024638, 0.51140279522436127, 0.49960931863170471, 0.48466451638746694, 0.47179215363312738]
# label_3 = "delay=2.0ms"
# label_shared_3 = "delay=2.0ms shared"
# 
# mean_fraction_shared_4 = [0.79481840747272847, 0.39469429747207524, 0.20350202495264225, 0.057868737344045997, 0.033990626428898038, 0.0092592592592592587, 0.0059670781893004111, 0.0016460905349794241, 0.00041152263374485601, 0.0]
# mean_fraction_4 = [0.89961705532693181, 0.70859276183350262, 0.62380292984716856, 0.55606403789016945, 0.53966373577430149, 0.52396312043228743, 0.51927184210947774, 0.51450870897887357, 0.50933343052896429, 0.50003430184539965]
# label_4 = "delay=3.0ms"
# label_shared_4 = "delay=3.0ms shared"
# 
# plt.plot(windows, mean_fraction_1, '-ro', label=label_1)
# plt.plot(windows, mean_fraction_2, '-go', label=label_2)
# plt.plot(windows, mean_fraction_3, '-bo', label=label_3)
# plt.plot(windows, mean_fraction_4, '-mo', label=label_4)
# 
# plt.plot(windows, mean_fraction_shared_1, '-rx', label=label_shared_1)
# plt.plot(windows, mean_fraction_shared_2, '-gx', label=label_shared_2)
# plt.plot(windows, mean_fraction_shared_3, '-bx', label=label_shared_3)
# plt.plot(windows, mean_fraction_shared_4, '-mx', label=label_shared_4)
# 
# plt.legend()
# plt.xlabel('Window size (ms)')
# plt.ylabel('Mean fraction')
# 
#==============================================================================

plt.show()
