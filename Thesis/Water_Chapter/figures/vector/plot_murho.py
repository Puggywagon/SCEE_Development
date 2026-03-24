#!/usr/bin/python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob


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

#csvfile = '../../data/pcm_dipoles.csv'
#csvfile = '../../TIP4P_Results/5ns/Analysis/data.csv'
csvfile = '../../TIP4P_Results/pcm_dipoles.csv'


################################################################################
################################################################################
################################################################################
T_list = [298, 400, 500, 600, 700, 800, 900, 1000]

df = pd.read_csv(csvfile, header=None)
df.columns = ['Temperature', 'Pressure', 'Density', 'Volume','Internal Energy','Enthalpy', 'Entropy','Cv','Cp','Sound','Joule','Viscosity','Therm','Phase','epsilon','g','PMC Dipole','My Dipole','Nearest','diff']
#df = pd.read_csv('../../TIP4P_Results/5ns/Analysis/data.csv', header=None)
#df.columns = ['Time', 'Temperature', 'Pressure', 'Density','Epsilon','Number of Honds','Percentage of H-bonds','Nearest Neighbour','Dipole Moment','State']


cmap = plt.cm.rainbow
color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
plt.figure(figsize=(5.5, 3.5))


for T, color in zip(T_list, color_list):
    tmp = df[df['Temperature'] == T]
    plt.plot(tmp['Density'], tmp['My Dipole'], ls='None',
             color=color, marker='o', markersize=markersize,
             label=r'$T='+f'{T}$ K')
#    plt.plot(tmp['Density'], tmp['PMC Dipole'], ls='',
#             color=color, marker='o', markerfacecolor='none', markersize=markersize,
#             markeredgewidth=1.0)

plt.ylabel(r'$\bar{\mu}$ / D',fontsize=fontsize)
plt.ylim([1.6, 3.0])
plt.xlabel(r'density / kg m$^{-3}$',fontsize=fontsize)
plt.xlim([0, 1200])
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)  

plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})   

################################################################################
plt.tight_layout()

plt.savefig('../murho.png')
plt.savefig('../murho.pdf')
if (show):
    plt.show()
else:    
    plt.clf()


################################################################################
################################################################################
################################################################################
