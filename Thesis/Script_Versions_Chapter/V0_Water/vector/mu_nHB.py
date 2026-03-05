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

fontsize = 10
ticksize = 10
markersize = 4
legendsize = 6

# required data files
csvfile = '../../TIP4P_Results/data_2.csv'

################################################################################
################################################################################
################################################################################
df = pd.read_csv(csvfile, header=None)
df.columns = ['Time', 'Temperature', 'Pressure', 'Density','Epsilon','Number of HBonds','Percentage of H-bonds','Nearest Neighbour','Dipole Moment','State','donors','acceptors','nHB','Cluster']
print(df)

nHB_min, nHB_max = 0, 4
donors_min, donors_max = 0, 4 
acceptors_min, acceptors_max = 0, 4


T_list = [298, 400, 500, 600, 700, 800, 900, 1000]
cmap = plt.cm.rainbow
color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
fig, ax = plt.subplots(figsize=(3.5, 3))

for T, color in zip(T_list, color_list):
    tmp = df[df['Temperature'] == T]
    ax.plot(tmp['nHB'], tmp['Dipole Moment'], ls='',
             color=color, marker='o', markersize=markersize, label=r'$T='+f'{T}$ K')

ax.set_xlabel(r'$\langle N_{\rm HB}\rangle$', fontsize=fontsize)    
ax.set_ylabel(r'$\bar{\mu}$ / D', fontsize=fontsize)
ax.set_xlim([0, 4])
ax.set_ylim([1.6, 3.0])
ax.set_xticks(range(nHB_min, nHB_max + 1))
ax.tick_params(axis="both", which="both", direction="in", labelsize=ticksize)  
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})

plt.tight_layout()
fig.subplots_adjust(hspace=0)

plt.savefig('../dip_nHB.png')
plt.savefig('../dip_nHB.pdf')
if (show):
    plt.show()
else:    
    plt.clf()


fig, ax = plt.subplots(figsize=(3.5, 3))

for T, color in zip(T_list, color_list):
    tmp = df[df['Temperature'] == T]
    ax.plot(tmp['donors'], tmp['Dipole Moment'], ls='',
             color=color, marker='o', markersize=markersize, label=r'$T='+f'{T}$ K')

ax.set_xlabel(r'$\langle N_{\rm Donors}\rangle$', fontsize=fontsize)    
ax.set_ylabel(r'$\bar{\mu}$ / D', fontsize=fontsize)
ax.set_xlim([0, 4])
ax.set_ylim([1.6, 3.0])
#ax.xlim([0, 1200])
ax.set_xticks(range(donors_min, donors_max + 1))
ax.tick_params(axis="both", which="both", direction="in", labelsize=ticksize)
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})  

plt.tight_layout()
fig.subplots_adjust(hspace=0)

plt.savefig('../dip_donors.png')
plt.savefig('../dip_donors.pdf')
if (show):
    plt.show()
else:    
    plt.clf()
    
    
fig, ax = plt.subplots(figsize=(3.5, 3))

for T, color in zip(T_list, color_list):
    tmp = df[df['Temperature'] == T]
    ax.plot(tmp['acceptors'], tmp['Dipole Moment'], ls='',
             color=color, marker='o', markersize=markersize, label=r'$T='+f'{T}$ K')

ax.set_xlabel(r'$\langle N_{\rm Acceptors}\rangle$', fontsize=fontsize)    
ax.set_ylabel(r'$\bar{\mu}$ / D', fontsize=fontsize)
ax.set_xlim([0, 4])
ax.set_ylim([1.6, 3.0])
#ax.xlim([0, 1200])
ax.set_xticks(range(acceptors_min, acceptors_max + 1))
ax.tick_params(axis="both", which="both", direction="in", labelsize=ticksize) 
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})   
plt.tight_layout()
fig.subplots_adjust(hspace=0)

plt.savefig('../dip_acceptors.png')
plt.savefig('../dip_acceptors.pdf')
if (show):
    plt.show()
else:    
    plt.clf()

################################################################################
################################################################################
################################################################################
