#!/usr/bash/env python

'''
Calculating 1/k_B * (<E_R> - <E_TS>) as a function of T for
an ensemble of isomers to see if the ensemble shows non-Arrhenius
behavior.
'''

import numpy as np

k = 8.617333262145e-5
E_R  = [0.0, 0.01, 0.1, 0.15, 0.2]
E_TS = [1.91, 0.06, 1.91, 1.76, 1.91]
T = [ i for i in range(300, 1001, 1)]
T_inv  = [1/i for i in T]
exp   = []
prob  = []
data  = []
output = []
exp_R  = []
exp_TS = []
prob_R  = [] 
prob_TS = []
slope = []

for i in range(len(T)):
    exp_R.append( [np.exp(-j/(k*T[i])) for j in E_R] )
    exp_TS.append( [np.exp(-j/(k*T[i])) for j in E_TS] )
for i in range(len(T)):
    prob_R.append( [icount/sum(exp_R[i]) for icount in exp_R[i]] )
    prob_TS.append( [icount/sum(exp_TS[i]) for icount in exp_TS[i]] )
for i in range(len(T)):
    temp  =      ( prob_R[i][0]*E_R[0]  + 
                   prob_R[i][1]*E_R[1]  +
                   prob_R[i][2]*E_R[2]  +
                   prob_R[i][3]*E_R[3]  +
                   prob_R[i][4]*E_R[4]  ) - ( prob_TS[i][0]*E_TS[0] +
                                              prob_TS[i][1]*E_TS[1] +
                                              prob_TS[i][2]*E_TS[2] +
                                              prob_TS[i][3]*E_TS[3] +
                                              prob_TS[i][4]*E_TS[4] )
    temp = temp/k
    slope.append([T[i], temp])

with open('slope.dat', 'w') as f:
    for i in range(len(T)):
        f.write(('%12.8f %12.8f \n') % (float(slope[i][0]), float(slope[i][1])))
