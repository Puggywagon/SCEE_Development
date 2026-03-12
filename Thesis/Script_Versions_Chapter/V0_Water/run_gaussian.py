#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pylab as plt

import SCEE

water = SCEE.SCEE()

water.g09root='/home/zoe/g09_pgi/'
water.GAUSS_SCRDIR = '/home/zoe/Research/Trial/TIP4P/Gaussian/scratch'

#water.rundir = './opt-orig/'
water.rundir = './opt-test/'

water.init_v0()
water.run_gaussian(step=0)

water.init_v1()
water.run_gaussian(step=1)
#
#multipole = water.read_multipoles('./opt-test/TIP4P2005_c1_q1_v1.out')
#print(multipole)

df = water.get_multipole_statistics()
print(df.describe())
print(df.head())

muL_list = []
for index, row in df.iterrows():
    x = np.array([2.0732078, 2.79883053, 3.52445326])
    y = np.array([row['dipole_l'], row['dipole_m'], row['dipole_h']])
    coeff = np.polyfit(x, y, 2)
    b = coeff[1] - 1
    muL = (-b - np.sqrt(b*b-4*coeff[0]*coeff[2])) / (2*coeff[0])
    muL_list.append(muL)
    
df['muL'] = muL_list
df.to_csv('junk_dipoles.csv', index=False)

Q1 = df['muL'].quantile(0.25)
Q3 =df['muL'].quantile(0.75)
IQR = Q3 - Q1
lower = Q1 - 1.5*IQR
upper = Q3 + 1.5*IQR
print(Q1, Q3, IQR)
print(lower, upper)



plt.hist(df['dipole_l'], 20)
plt.hist(df['dipole_m'], 20)
plt.hist(df['dipole_h'], 20)
plt.show()
