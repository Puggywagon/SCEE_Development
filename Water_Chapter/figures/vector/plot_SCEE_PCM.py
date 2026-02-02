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

fontsize = 12
ticksize = 18

csvfile = '../../TIP4P_Results/pcm_dipoles.csv'
#csvfile = '../../data/pcm_dipoles.csv'


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
    plt.plot(tmp['PMC Dipole'], tmp['My Dipole'], ls='',
             color=color, marker='o',label=r'$T='+f'{T}$ K')

    
plt.plot([1.6, 3.0], [1.6, 3.0], color='black', ls='solid')
plt.xlabel(r'$\bar{\mu}_{\rm PCM}$ / D', fontsize=fontsize)
plt.ylabel(r'$\bar{\mu}_{\rm SCEE}$ / D', fontsize=fontsize)
plt.xlim([1.8, 2.2])
plt.ylim([1.6, 3.0])
#plt.xlim([0, 1200])
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)  


plt.tight_layout()

plt.savefig('../SCEE_PCM.png')
plt.savefig('../SCEE_PCM.pdf')
if (show):
    plt.show()
else:    
    plt.clf()


################################################################################
################################################################################
################################################################################
