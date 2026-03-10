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

# required data files
csvfile = '../../TIP4P_Results/data_2.csv'


################################################################################
################################################################################
################################################################################
df = pd.read_csv(csvfile, header=None)
df.columns = ['Time', 'Temperature', 'Pressure', 'Density','Epsilon','Number of HBonds','Percentage of H-bonds','Nearest Neighbour','Dipole Moment','State','donors','acceptors','nHB']


T_list = [298, 400, 500, 600, 700, 800, 900, 1000]


plt.figure(figsize=(3.5, 3))

cmap = plt.cm.rainbow
color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]


for T, color in zip(T_list, color_list):
    tmp = df[df['Temperature'] == T]
    plt.plot(tmp['nHB'], tmp['Dipole Moment'], ls='',
             color=color, marker='o', markersize=markersize, label=r'$T='+f'{T}$ K')



plt.xlabel(r'$\langle N_{\rm HB}\rangle$', fontsize=fontsize)    
plt.ylabel(r'$\bar{\mu}$ / D', fontsize=fontsize)
plt.xlim([0, 4])
plt.ylim([1.65, 3.75])
#plt.xlim([0, 1200])
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)  


plt.tight_layout()

plt.savefig('../dip_nHB.png')
plt.savefig('../dip_nHB.pdf')
if (show):
    plt.show()
else:    
    plt.clf()


plt2.figure(figsize=(3.5, 3))



for T, color in zip(T_list, color_list):
    tmp = df[df['Temperature'] == T]
    plt2.plot(tmp['donors'], tmp['Dipole Moment'], ls='',
             color=color, marker='o', markersize=markersize, label=r'$T='+f'{T}$ K')

plt2.xlabel(r'$\langle N_{\rm Donors}\rangle$', fontsize=fontsize)    
plt2.ylabel(r'$\bar{\mu}$ / D', fontsize=fontsize)
plt2.xlim([0, 4])
plt2.ylim([1.65, 3.75])
#plt.xlim([0, 1200])
plt2.xticks(fontsize=fontsize)
plt2.yticks(fontsize=fontsize)  


plt2.tight_layout()

plt2.savefig('../dip_donors.png')
plt2.savefig('../dip_donors.pdf')
if (show):
    plt2.show()
else:    
    plt2.clf()
    
plt3.figure(figsize=(3.5, 3))



for T, color in zip(T_list, color_list):
    tmp = df[df['Temperature'] == T]
    plt3.plot(tmp['acceptors'], tmp['Dipole Moment'], ls='',
             color=color, marker='o', markersize=markersize, label=r'$T='+f'{T}$ K')


plt3.xlabel(r'$\langle N_{\rm Acceptors}\rangle$', fontsize=fontsize)    
plt3.ylabel(r'$\bar{\mu}$ / D', fontsize=fontsize)
plt3.xlim([0, 4])
plt3.ylim([1.65, 3.75])
#plt.xlim([0, 1200])
plt3.xticks(fontsize=fontsize)
plt3.yticks(fontsize=fontsize)  


plt3.tight_layout()

plt3.savefig('../dip_acceptors.png')
plt3.savefig('../dip_acceptors.pdf')
if (show):
    plt3.show()
else:    
    plt3.clf()

################################################################################
################################################################################
################################################################################
