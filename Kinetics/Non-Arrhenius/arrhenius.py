#!/usr/bash/env python

'''
Exploring all Eai combinations for ensemble of 5 isomers to find the combinations which show 
non-Arrhenius behavior. The relative energies can be set initially but they are fixed
during the search.
'''

import numpy as np
import timeit
from scipy import stats
from scipy.interpolate import spline
from itertools import chain

#|===========|
#|   INPUT   |
#|===========|

k = 8.617333262145e-5
E = [0.0, 0.01, 0.1, 0.15, 0.2]
dx1 = 0.01
dx2 = 0.1
Ea1 = 0.1
Ea2 = 2.01

T = [300, 400, 500, 600, 700, 800, 900, 1000]
T_inv  = [1/i for i in T]
exp   = []
prob  = []
data  = []
output = []

start_time = timeit.default_timer()

#define the initial and final Eai and the step in each part to be explored
loop1 = np.array(list(np.arange(0, Ea1, dx1)) + list(np.arange(Ea1, Ea2, dx2)))

for i in range(len(T)):
    exp.append( [np.exp(-j/(k*T[i])) for j in E] )
for i in range(len(T)):
    prob.append( [icount/sum(exp[i]) for icount in exp[i]] )

for i1 in loop1:
    for i2 in loop1:
        for i3 in loop1:
            for i4 in loop1:
                for i5 in loop1:
                    logk = []
                    for i in range(len(T)):
                        kens =     (  prob[i][0]*np.exp(-i1/(k*T[i])) + 
                                      prob[i][1]*np.exp(-i2/(k*T[i])) +
                                      prob[i][2]*np.exp(-i3/(k*T[i])) +
                                      prob[i][3]*np.exp(-i4/(k*T[i])) +
                                      prob[i][4]*np.exp(-i5/(k*T[i])) )
                        data.append([T[i], kens, i1, i2, i3, i4, i5])
                        logk.append(np.log(kens))
                    coeffs1 = np.polyfit(T_inv, logk, 1)
                    poly1   = np.poly1d(coeffs1)
                    slope, intercept, r_value, p_value, std_err = stats.linregress(T_inv, logk)
                    output.append([i1, i2, i3, i4, i5, slope, intercept, r_value**2])

output.sort(key = lambda x: x[7])

with open('output', 'w') as f:
    f.write('%19s %4.2f %4.2f %4.2f %4.2f %4.2f\n' %
            ('Isomer Eenergies = ', E[0], E[1], E[2], E[3], E[4]))
    f.write('%14s %4.2f %2s %4.2f %2s %4.2f %7s %4.2f %2s %4.2f \n' %
            ('Ea explored = ', '0.0', '--', Ea1, '--', Ea2, 'Steps = ', dx1,', ', dx2))
    f.write('%4s %4s %4s %4s %4s %6s %9s %3s\n' %
            ('Ea1', 'Ea2', 'Ea3', 'Ea4', 'Ea5', 'slope', 'intercept', 'R2'))
    for i in range(len(output)):
        f.write('%4.2f %4.2f %4.2f %4.2f %4.2f %8.4f %8.4f %8.4f\n' %
               (output[i][0], output[i][1], output[i][2], output[i][3], output[i][4], output[i][5], output[i][6], output[i][7]))

print('Total time in seconds =', timeit.default_timer() - start_time)
