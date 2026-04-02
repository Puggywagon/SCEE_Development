#!/usr/bin/python3

import argparse
import subprocess
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt


################################################################################
################################################################################
################################################################################
parser = argparse.ArgumentParser()

parser.add_argument('--show', type=str, default='True', help='plot figures')

args = parser.parse_args()

show = (args.show == 'True')

fontsize = 10
legendsize = 6
markersize = 4

#nist_csvfile = './NIST_isotherms.csv'
#pcm_csvfile = '../../data/pcm_dipoles.csv'
nist_csvfile = './NIST_isotherms.csv'
pcm_csvfile = '../../TIP4P_Results/5ns/Analysis/data.csv'


################################################################################
################################################################################
################################################################################
T_list = [298, 400, 500, 600, 700, 800, 900, 1000]

df_nist = pd.read_csv(nist_csvfile)
df = pd.read_csv(pcm_csvfile, header=None)
#df.columns = ['Temperature', 'Pressure', 'Density', 'Volume','Internal Energy','Enthalpy', 'Entropy','Cv','Cp','Sound','Joule','Viscosity','Therm','Phase','epsilon','g','PMC Dipole','My Dipole','Nearest','diff']

#df_nist = pd.read_csv('../../TIP4P_Results/NIST_isotherms.csv')
#df = pd.read_csv('../../TIP4P_Results/5ns/Analysis/data.csv', header=None)
df = pd.read_csv(pcm_csvfile, header=None)
df.columns = ['Time', 'Temperature', 'Pressure', 'Density','Epsilon','Number of Honds','Percentage of H-bonds','Nearest Neighbour','Dipole Moment','State']
print(df_nist.columns)
print(df.columns)

cmap = plt.cm.rainbow
color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
#plt.figure(figsize=(3.5, 3))
plt.figure(figsize=(4, 3))

for T, color in zip(T_list, color_list):
    tmp = df_nist[df_nist['Temperature (K)'] == T]
    plt.plot(tmp['Pressure (bar)'], tmp['Density (kg/m3)'],
             color=color)

    tmp = df[df['Temperature']==T]
    plt.plot(tmp['Pressure'], tmp['Density'],
             marker='o', markersize=markersize, color=color,
             ls='None',
             label=f'$T={T} K$')
    print(f'T={T}')
    print(tmp['Pressure'])

    
plt.ylabel(r'density / kg m$^{-3}$',fontsize=fontsize)
plt.xlabel(r'pressure / bar',fontsize=fontsize)
plt.xscale('log')
plt.xlim([0.5, 1.0e4])
plt.ylim([0, 1200])
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)  
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
           prop={'size':legendsize})


################################################################################
plt.tight_layout()

plt.savefig('../density.png')
plt.savefig('../density.pdf')
if (show):
    plt.show()
else:    
    plt.clf()


################################################################################
################################################################################
################################################################################
