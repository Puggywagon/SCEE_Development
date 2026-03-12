#!/usr/bin/python3

import argparse
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
def read_data(T, p_list):
    csvfile_list = []
    for p in p_list:
        csvfile_list += glob.glob(f'../../Replica_Dipoles/Replica_*/{T}K/{p:.1f}Bar/Dipole.csv')
    df = pd.DataFrame()
    for csvfile in csvfile_list:
        tmp = pd.read_csv(csvfile)
        df = pd.concat([df, tmp])
    return df


################################################################################
################################################################################
################################################################################
parser = argparse.ArgumentParser()

parser.add_argument('--show', type=str, default='True', help='plot figures')

args = parser.parse_args()

show = (args.show == 'True')


fontsize = 22
ticksize = 16
labelsize = 18


################################################################################
################################################################################
################################################################################
T_list = [298, 400, 500, 600, 700, 800, 900, 1000]
p_list = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]


cmap = plt.cm.rainbow
color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]


################################################################################
T = 298
color = color_list[0]
print(f'temperature = {T} K')
p_list = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
df = read_data(T, p_list)
custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 0.2, num=20) 
counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
plt.plot(x_list, counts, color=color, label=f'T={T:.0f} K')


    
################################################################################
T = 400
color = color_list[1]
print(f'temperature = {T} K')
p_list = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
df = read_data(T, p_list)
custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 1, num=20) 
counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
plt.plot(x_list, counts, color=color, label=f'T={T:.0f} K')



################################################################################
T = 500
color = color_list[2]
p_list = [10, 20, 50, 100, 200, 500, 1000]
print(f'temperature = {T} K')
df = read_data(T, p_list)
custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 1, num=20)
counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
plt.plot(x_list, counts, color=color, label=f'T={T:.0f} K')

    
################################################################################
T = 600
color = color_list[4]
p_list = [100, 200, 500, 1000]
print(f'temperature = {T} K')
df = read_data(T, p_list)
custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 1, num=20)
counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
plt.plot(x_list, counts, color=color, label=f'T={T:.0f} K')
    
#ax3.text(0.9, 0.9, '(c)', transform=ax3.transAxes, fontsize=fontsize)
#ax3.tick_params(axis='x', which='both', labelsize=ticksize)
#ax3.tick_params(axis='y',which='both',labelsize=ticksize) 
#ax3.set_xlim([1.65,3.75])
plt.xlabel(r'$\mu$ / D', fontsize=fontsize)
plt.ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
plt.legend(loc='upper right')
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)  


################################################################################
plt.tight_layout()


plt.savefig('../dipole_distliq.pdf')
plt.savefig('../dipole_distliq.png')

if (show):
    plt.show()
else:    
    plt.clf()



################################################################################
################################################################################
################################################################################
