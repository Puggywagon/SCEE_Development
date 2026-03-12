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
def read_data(T, p):
    csvfile_list = glob.glob(f'../../Replica_Dipoles/Replica_*/{T}K/{p:.1f}Bar/Dipole.csv')
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


fontsize = 10
ticksize = 10
legendsize = 6
################################################################################
################################################################################
################################################################################
T_list = [298, 400, 500, 600, 700, 800, 900, 1000]
p_list = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]


cmap = plt.cm.rainbow
color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(p_list))]
################################################################################
fig, ax = plt.subplots(figsize=(3.5, 3))
T = 298
print(f'temperature = {T} K')
for p, color in zip(p_list, color_list):
    df = read_data(T, p)

    custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 0.2, num=20) 
    counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
    x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
    ax.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')

ax.tick_params(axis='x', which='both', direction='in',
                labelsize=ticksize, labelbottom=ticksize)
ax.tick_params(axis='y',which='both',labelsize=ticksize) 
ax.set_xlim([1.65,3.75])
ax.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
ax.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})
plt.tight_layout()


plt.savefig(f'../dipole_distribution_{T}.pdf')
plt.savefig(f'../dipole_distribution_{T}.png')

if (show):
    plt.show()
else:    
    plt.clf()
################################################################################
fig, ax = plt.subplots(figsize=(3.5, 3))
T = 400
print(f'temperature = {T} K')
for p, color in zip(p_list, color_list):
    df = read_data(T, p)

    custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 0.2, num=20) 
    counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
    x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
    ax.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')

ax.tick_params(axis='x', which='both', direction='in',
                labelsize=ticksize, labelbottom=ticksize)
ax.tick_params(axis='y',which='both',labelsize=ticksize) 
ax.set_xlim([1.65,3.75])
ax.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
ax.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})
plt.tight_layout()

plt.savefig(f'../dipole_distribution_{T}.pdf')
plt.savefig(f'../dipole_distribution_{T}.png')

if (show):
    plt.show()
else:    
    plt.clf()
################################################################################
fig, ax = plt.subplots(figsize=(3.5, 3))
T = 500
print(f'temperature = {T} K')
for p, color in zip(p_list, color_list):
    df = read_data(T, p)

    custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 0.2, num=20) 
    counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
    x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
    ax.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')

ax.tick_params(axis='x', which='both', direction='in',
                labelsize=ticksize, labelbottom=ticksize)
ax.tick_params(axis='y',which='both',labelsize=ticksize) 
ax.set_xlim([1.65,3.75])
ax.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
ax.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})
plt.tight_layout()

plt.savefig(f'../dipole_distribution_{T}.pdf')
plt.savefig(f'../dipole_distribution_{T}.png')

if (show):
    plt.show()
else:    
    plt.clf()
################################################################################
fig, ax = plt.subplots(figsize=(3.5, 3))
T = 600
print(f'temperature = {T} K')
for p, color in zip(p_list, color_list):
    df = read_data(T, p)
    custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 1, num=20) 
    counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
    x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
    ax.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')

ax.tick_params(axis='x', which='both', direction='in',
                labelsize=ticksize, labelbottom=ticksize)
ax.tick_params(axis='y',which='both',labelsize=ticksize) 
ax.set_xlim([1.65,3.75])
ax.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
ax.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})
plt.tight_layout()

plt.savefig(f'../dipole_distribution_{T}.pdf')
plt.savefig(f'../dipole_distribution_{T}.png')

if (show):
    plt.show()
else:    
    plt.clf()
################################################################################
fig, ax = plt.subplots(figsize=(3.5, 3))
T = 700
print(f'temperature = {T} K')
for p, color in zip(p_list, color_list):
    df = read_data(T, p)

    custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 0.2, num=20) 
    counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
    x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
    ax.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')

ax.tick_params(axis='x', which='both', direction='in',
                labelsize=ticksize, labelbottom=ticksize)
ax.tick_params(axis='y',which='both',labelsize=ticksize) 
ax.set_xlim([1.65,3.75])
ax.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
ax.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
ax.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})
plt.tight_layout()

plt.savefig(f'../dipole_distribution_{T}.pdf')
plt.savefig(f'../dipole_distribution_{T}.png')

if (show):
    plt.show()
else:    
    plt.clf()
################################################################################
fig, ax = plt.subplots(figsize=(3.5, 3))
T = 800
print(f'temperature = {T} K')
for p, color in zip(p_list, color_list):
    df = read_data(T, p)
    custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 1, num=20)
    counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
    x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
    ax.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')
    
ax.tick_params(axis='x', which='both', direction='in',
                labelsize=ticksize, labelbottom=ticksize)
ax.tick_params(axis='y',which='both',labelsize=ticksize) 
ax.set_xlim([1.65,3.75])
ax.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
ax.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})
plt.tight_layout()

plt.savefig(f'../dipole_distribution_{T}.pdf')
plt.savefig(f'../dipole_distribution_{T}.png')

if (show):
    plt.show()
else:    
    plt.clf()
################################################################################
fig, ax = plt.subplots(figsize=(3.5, 3))
T = 900
print(f'temperature = {T} K')
for p, color in zip(p_list, color_list):
    df = read_data(T, p)
    custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 1, num=20) 
    counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
    x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
    ax.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')

ax.tick_params(axis='x', which='both', direction='in',
                labelsize=ticksize, labelbottom=ticksize)
ax.tick_params(axis='y',which='both',labelsize=ticksize) 
ax.set_xlim([1.65,3.75])
ax.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
ax.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})
plt.tight_layout()

plt.savefig(f'../dipole_distribution_{T}.pdf')
plt.savefig(f'../dipole_distribution_{T}.png')

if (show):
    plt.show()
else:    
    plt.clf()
################################################################################
fig, ax = plt.subplots(figsize=(3.5, 3))
T = 1000
print(f'temperature = {T} K')
for p, color in zip(p_list, color_list):
    df = read_data(T, p)
    custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 1, num=20)
    counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
    x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
    ax.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')
    
ax.tick_params(axis='x', which='both', direction='in',
                labelsize=ticksize, labelbottom=ticksize)
ax.tick_params(axis='y',which='both',labelsize=ticksize) 
ax.set_xlim([1.65,3.75])
ax.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
ax.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})
plt.tight_layout()

plt.savefig(f'../dipole_distribution_{T}.pdf')
plt.savefig(f'../dipole_distribution_{T}.png')

if (show):
    plt.show()
else:    
    plt.clf()
################################################################################
