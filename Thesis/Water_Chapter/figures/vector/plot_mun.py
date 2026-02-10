#!/usr/bin/python3

import argparse
import subprocess
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize


################################################################################
################################################################################
################################################################################
parser = argparse.ArgumentParser()

parser.add_argument('--show', type=str, default='True', help='plot figures')

args = parser.parse_args()

show = (args.show == 'True')

fontsize = 12
ticksize = 18
markersize = 4
legendsize = 6

# required data files
csvfile = '../../TIP4P_Results/pcm_dipoles.csv'
#csvfile = '../../data/pcm_dipoles.csv'
#csvfile = './pcm_dipoles.csv'


################################################################################
################################################################################
################################################################################
df = pd.read_csv(csvfile, header=None)
df.columns = ['Temperature', 'Pressure', 'Density', 'Volume','Internal Energy','Enthalpy', 'Entropy','Cv','Cp','Sound','Joule','Viscosity','Therm','Phase','epsilon','g','PMC Dipole','My Dipole','Nearest','diff']


T_list = [298, 400, 500, 600, 700, 800, 900, 1000]


plt.figure(figsize=(3.5, 3))

cmap = plt.cm.rainbow
color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]


for T, color in zip(T_list, color_list):
    tmp = df[df['Temperature'] == T]
    plt.plot(tmp['Nearest'], tmp['My Dipole'], ls='',
             color=color, marker='o', markersize=markersize, label=r'$T='+f'{T}$ K')
    plt.plot(tmp['Nearest'], tmp['PMC Dipole'], ls='',
             color=color, marker='o', markersize=markersize,
             markerfacecolor='none', markeredgewidth=1.0)


mu0 = 1.856

def res(slope):
    err = 0.0
    for index, row in df.iterrows():
        mu = row['My Dipole']
        n = row['Nearest']
        mu_fit = mu0 + slope*n
        err += (mu_fit-mu)**2
    return err
sol = minimize(res, 1.0)
slope = sol['x']

def res(slope):
    err = 0.0
    for index, row in df.iterrows():
        n = row['Nearest']
        if (n <= 1.0):
            mu = row['My Dipole']
            mu_fit = mu0 + slope*n
            err += (mu_fit-mu)**2
    return err
sol = minimize(res, 1.0)
slope1 = sol['x']
print(sol)

plt.plot([0, 4], [mu0+slope*0, mu0+slope*4], color='black', ls='solid')
#plt.plot([0, 4], [mu0+slope1*0, mu0+slope1*4], color='black', ls='solid')


plt.xlabel(r'$\langle N_{\rm neig}\rangle$', fontsize=fontsize)    
plt.ylabel(r'$\bar{\mu}$ / D', fontsize=fontsize)
plt.xlim([0, 4])
plt.ylim([1.6, 3.0])
#plt.xlim([0, 1200])
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)  

plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})   

plt.tight_layout()

plt.savefig('../mun.png')
plt.savefig('../mun.pdf')
if (show):
    plt.show()
else:    
    plt.clf()


################################################################################
################################################################################
################################################################################
