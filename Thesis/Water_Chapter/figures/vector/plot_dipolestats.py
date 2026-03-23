#!/usr/bin/python3

import argparse
import numpy as np
import os
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import glob


################################################################################
################################################################################
parser = argparse.ArgumentParser()
parser.add_argument('--show', type=str, default='True', help='plot figures')
args = parser.parse_args()
show = (args.show == 'True')
 
fontsize = 12
ticksize = 12
labelsize = 18
markersize = 4
legendsize = 6
 
csvfile = '../../Replica_Dipoles/Dipoles2.csv'
lines_csvfile = './Dipoles2_lines.csv'
 
################################################################################
T_list = [298, 400, 500, 600, 700, 800, 900, 1000]
 
df = pd.read_csv(csvfile, header=None)
df.columns = ['Temperature', 'Pressure', 'Dipoles', 'SE', 'std means']
df_lines = pd.read_csv(lines_csvfile)
 
cmap = plt.cm.rainbow
color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
 
# Side-by-side layout — extra width on right for the legend
fig, ax = plt.subplots(1, 2, figsize=(7.5, 3.5))
ax1 = ax[0]
ax2 = ax[1]
 
################################################################################
# Panel (a) — mean dipole moment
for T, color in zip(T_list, color_list):
    tmp = df[df['Temperature'] == T]
    yerr = [tmp['SE'], tmp['SE']]
    ax1.errorbar(tmp['Pressure'], tmp['Dipoles'], yerr=yerr, capsize=3,
                 ecolor="black", ls='none', color=color, marker='o',
                 markersize=markersize)
    if T == 500:
        tmp_lines = df_lines[(df_lines['T'] == T) & (df_lines['P'] >= 10)]
        ax1.plot(tmp_lines['P'], tmp_lines['dipole line'], ls='dashed', color=color, lw=2)
        tmp_lines = df_lines[(df_lines['T'] == T) & (df_lines['P'] < 10)]
        ax1.plot(tmp_lines['P'], tmp_lines['dipole line'], ls='dashed', color=color, lw=2)
    else:
        tmp_lines = df_lines[df_lines['T'] == T]
        ax1.plot(tmp_lines['P'], tmp_lines['dipole line'], ls='dashed', color=color, lw=2)
 
ax1.hlines(2.3054, 0.1, 1.0e4, ls='dotted', color='black')
 
ax1.text(0.05, 0.95, '(a)', transform=ax1.transAxes, fontsize=fontsize,
         va='top', ha='left')
ax1.set_ylabel(r'$\bar{\mu}$ / D', fontsize=fontsize)
ax1.set_xlabel(r'Pressure / bar', fontsize=fontsize)
ax1.tick_params(axis='x', which='both', direction='in', labelsize=ticksize)
ax1.tick_params(axis='y', which='both', labelsize=ticksize)
ax1.set_xscale('log')
ax1.set_xlim([0.5, 1.1e4])
ax1.set_ylim([1.6, 3.0])
 
################################################################################
# Panel (b) — dipole distribution width
for T, color in zip(T_list, color_list):
    tmp = df[df['Temperature'] == T]
    ax2.plot(tmp['Pressure'], tmp['std means'], ls='none', color=color,
             marker='o', markersize=markersize,
             label=r'$T=' + f'{T}$ K')
    if T == 500:
        tmp_lines = df_lines[(df_lines['T'] == T) & (df_lines['P'] >= 10)]
        ax2.plot(tmp_lines['P'], tmp_lines['Stdev line'], ls='dashed', color=color)
        tmp_lines = df_lines[(df_lines['T'] == T) & (df_lines['P'] < 10)]
        ax2.plot(tmp_lines['P'], tmp_lines['Stdev line'], ls='dashed', color=color)
    else:
        tmp_lines = df_lines[df_lines['T'] == T]
        ax2.plot(tmp_lines['P'], tmp_lines['Stdev line'], ls='dashed', color=color)
 
ax2.text(0.05, 0.95, '(b)', transform=ax2.transAxes, fontsize=fontsize,
         va='top', ha='left')
ax2.set_ylabel(r'$\sigma_\mu$ / D', fontsize=fontsize)
ax2.set_xlabel(r'Pressure / bar', fontsize=fontsize)
ax2.tick_params(axis='x', which='both', direction='in', labelsize=ticksize)
ax2.tick_params(axis='y', which='both', labelsize=ticksize)
ax2.set_xscale('log')
ax2.set_xlim([0.5, 1.1e4])
ax2.set_ylim([0, 0.4])
 
# Legend outside to the right of panel (b)
ax2.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize},
           frameon=True, framealpha=0.8, borderaxespad=0)
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
