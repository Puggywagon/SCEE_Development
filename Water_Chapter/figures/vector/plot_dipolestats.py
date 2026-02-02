#!/usr/bin/python3

import argparse
import gromacs
import subprocess
import numpy as np
import os
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import glob


################################################################################
################################################################################
################################################################################
parser = argparse.ArgumentParser()

parser.add_argument('--show', type=str, default='True', help='plot figures')

args = parser.parse_args()

show = (args.show == 'True')

#fontsize = 22
fontsize = 12
ticksize = 12
labelsize = 18
markersize = 4

csvfile = '../../Replica_Dipoles/Dipoles2.csv'
#csvfile = '../../data/Dipoles2.csv'
lines_csvfile = './Dipoles2_lines.csv'


################################################################################
################################################################################
################################################################################
T_list = [298, 400, 500, 600, 700, 800, 900, 1000]

df = pd.read_csv(csvfile, header=None)
df.columns = ['Temperature', 'Pressure', 'Dipoles','SE','std means']
df_lines = pd.read_csv(lines_csvfile)

cmap = plt.cm.rainbow
color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]

#fig= plt.subplots(figsize=(15, 12))
fig, ax= plt.subplots(2, 1, figsize=(3.5, 5))
ax1 = ax[0]
ax2 = ax[1]
#ax1b = ax1.twinx()


################################################################################
for T, color in zip(T_list, color_list):
    tmp = df[df['Temperature'] == T]
    minval=tmp['SE']
    maxval=tmp['SE']
    corrections=tmp['Dipoles']/2.3054
    yerr=[minval,maxval]
    ax1.errorbar(tmp['Pressure'], tmp['Dipoles'], yerr=yerr, capsize=3, ecolor = "black", ls='none', color=color, marker='o', markersize=markersize)
#    ax1b.plot(tmp['Pressure'], corrections, ls='none', color=color, marker='s')
    plt.plot([], [], ls='', color=color, marker='o', markersize=markersize,
             label=r'$T='+f'{T}$ K')    
    if (T ==500):
        tmp = df_lines[(df_lines['T']==T) & (df_lines['P']>=10)]
        ax1.plot(tmp['P'], tmp['dipole line'], ls='dashed', color=color, lw=2)
        tmp = df_lines[(df_lines['T']==T) & (df_lines['P']<10)]
        ax1.plot(tmp['P'], tmp['dipole line'], ls='dashed', color=color, lw=2)
    else:
        tmp = df_lines[df_lines['T']==T]
        ax1.plot(tmp['P'], tmp['dipole line'], ls='dashed', color=color, lw=2)

#ax1.plot([0.5, 1.0e4], [2.3054, 2.3054], ls='dotted', color='black')
ax1.hlines(2.3054, 0.1, 1.0e4, ls='dotted', color='black')

ax1.text(0.9, 0.1, '(a)', transform=ax1.transAxes, fontsize=fontsize)
ax1.set_ylabel(r'$\bar{\mu}$ / D',fontsize=fontsize)
ax1.tick_params(axis='x', which='both', direction='in',
                labelsize=ticksize, labelbottom=False)
ax1.tick_params(axis='y',which='both',labelsize=ticksize) 
ax1.set_xscale('log')
ax1.set_xlim([0.5, 1.1e4])
#ax1b.set_ylabel(r'Dipole Moment / D',fontsize=fontsize)
#ax1b.set_ylabel(r'$\bar{\mu} / \mu_{\rm TIP4P}$',fontsize=fontsize)
#ax1b.set_xlabel(r'Pressure / bar',fontsize=fontsize)
#ax1b.tick_params(axis='y',which='both',labelsize=ticksize) 
#ax1b.set_ylim([min_correction, max_correction])
mn,mx=ax1.set_ylim([1.6, 3.0])
min_correction=mn/2.3054
max_correction=mx/2.3054


################################################################################
for T, color in zip(T_list, color_list):
    tmp = df[df['Temperature'] == T]
    ax2.plot(tmp['Pressure'],  tmp['std means'], ls='none', color=color,
             marker='o', markersize=markersize)
    ax2.plot([], [], ls='', color=color, marker='o', markersize=markersize,
             label=r'$T='+f'{T}$ K')
    if (T ==500):
        tmp = df_lines[(df_lines['T']==T) & (df_lines['P']>=10)]
        ax2.plot(tmp['P'], tmp['Stdev line'], ls='dashed', color=color)
        tmp = df_lines[(df_lines['T']==T) & (df_lines['P']<10)]
        ax2.plot(tmp['P'], tmp['Stdev line'], ls='dashed', color=color)
    else:
        tmp = df_lines[df_lines['T']==T]
        ax2.plot(tmp['P'], tmp['Stdev line'], ls='dashed', color=color)
#    
#    a,b=np.polyfit(tmp['Pressure'],tmp['std means'],1)
#    x=tmp['Pressure']    
#    ax2.plot(x,(a*x+b), ls='--',color=color)
#ax2.set_ylabel(r'Dipole Moment Distribution Width / D',fontsize=fontsize)
ax2.text(0.9, 0.1, '(b)', transform=ax2.transAxes, fontsize=fontsize)
ax2.set_ylabel(r'$\sigma_\mu$ / D',fontsize=fontsize)
ax2.set_xlabel(r'Pressure / bar',fontsize=fontsize)
ax2.tick_params(axis='x', which='both', direction='in',
                labelsize=ticksize, labelbottom=True)
ax2.tick_params(axis='y', which='both', labelsize=ticksize)
#plt.xticks(fontsize=ticksize)
#plt.yticks(fontsize=ticksize) 
ax2.set_xscale('log')
ax2.set_xlim([0.5, 1.1e4])
ax2.set_ylim([0, 0.4])



################################################################################
plt.tight_layout()
fig.subplots_adjust(hspace=0)

plt.savefig('../dipole_statistics.pdf')
plt.savefig('../dipole_statistics.png')

if (show):
    plt.show()
else:    
    plt.clf()



################################################################################
################################################################################
################################################################################
