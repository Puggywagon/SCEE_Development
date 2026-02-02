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
plt1.figure(figsize=(3.5, 3))
T = 298
print(f'temperature = {T} K')
for p, color in zip(p_list, color_list):
    df = read_data(T, p)

    custom_bins = np.linspace(start=1.65, stop=mplt(df['muL']) + 0.2, num=20) 
    counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
    x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
    plt1.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')

plt1.text(0.9, 0.9, '(a)', transform=plt1.transpltes, fontsize=fontsize)
plt1.tick_params(pltis='x', which='both', direction='in',
                labelsize=ticksize, labelbottom=False)
plt1.tick_params(pltis='y',which='both',labelsize=ticksize) 
plt1.set_xlim([1.65,3.75])
plt1.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
plt1.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
plt1.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
plt1.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})
plt1.tight_layout()


plt1.savefig(f'../dipole_distribution_{T}.pdf')
plt1.savefig(f'../dipole_distribution_{T}.png')

if (show):
    plt1.show()
else:    
    plt1.clf()
################################################################################
plt2.figure(figsize=(3.5, 3))
T = 400
print(f'temperature = {T} K')
for p, color in zip(p_list, color_list):
    df = read_data(T, p)

    custom_bins = np.linspace(start=1.65, stop=mplt(df['muL']) + 0.2, num=20) 
    counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
    x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
    plt2.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')

plt2.tick_params(pltis='x', which='both', direction='in',
                labelsize=ticksize, labelbottom=False)
plt2.tick_params(pltis='y',which='both',labelsize=ticksize) 
plt2.set_xlim([1.65,3.75])
plt2.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
plt2.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
plt2.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
plt2.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})
plt2.tight_layout()

plt2.savefig(f'../dipole_distribution_{T}.pdf')
plt2.savefig(f'../dipole_distribution_{T}.png')

if (show):
    plt2.show()
else:    
    plt2.clf()
################################################################################
plt3.figure(figsize=(3.5, 3))
T = 500
print(f'temperature = {T} K')
for p, color in zip(p_list, color_list):
    df = read_data(T, p)

    custom_bins = np.linspace(start=1.65, stop=mplt(df['muL']) + 0.2, num=20) 
    counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
    x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
    plt3.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')

plt3.tick_params(pltis='x', which='both', direction='in',
                labelsize=ticksize, labelbottom=False)
plt3.tick_params(pltis='y',which='both',labelsize=ticksize) 
plt3.set_xlim([1.65,3.75])
plt3.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
plt3.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
plt3.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
plt3.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})
plt3.tight_layout()

plt3.savefig(f'../dipole_distribution_{T}.pdf')
plt3.savefig(f'../dipole_distribution_{T}.png')

if (show):
    plt3.show()
else:    
    plt3.clf()
################################################################################
plt4.figure(figsize=(3.5, 3))
T = 600
print(f'temperature = {T} K')
for p, color in zip(p_list, color_list):
    df = read_data(T, p)
    custom_bins = np.linspace(start=1.65, stop=mplt(df['muL']) + 1, num=20) 
    counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
    x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
    plt4.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')

plt4.tick_params(pltis='x', which='both', direction='in',
                labelsize=ticksize, labelbottom=False)
plt4.tick_params(pltis='y',which='both',labelsize=ticksize) 
plt4.set_xlim([1.65,3.75])
plt4.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
plt4.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
plt4.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
plt4.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})
plt4.tight_layout()

plt4.savefig(f'../dipole_distribution_{T}.pdf')
plt4.savefig(f'../dipole_distribution_{T}.png')

if (show):
    plt4.show()
else:    
    plt4.clf()
################################################################################
plt5.figure(figsize=(3.5, 3))
T = 700
print(f'temperature = {T} K')
for p, color in zip(p_list, color_list):
    df = read_data(T, p)

    custom_bins = np.linspace(start=1.65, stop=mplt(df['muL']) + 0.2, num=20) 
    counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
    x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
    plt5.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')

plt5.tick_params(pltis='x', which='both', direction='in',
                labelsize=ticksize, labelbottom=False)
plt5.tick_params(pltis='y',which='both',labelsize=ticksize) 
plt5.set_xlim([1.65,3.75])
plt5.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
plt5.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
plt5.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
plt5.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
plt5.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})
plt5.tight_layout()

plt5.savefig(f'../dipole_distribution_{T}.pdf')
plt5.savefig(f'../dipole_distribution_{T}.png')

if (show):
    plt5.show()
else:    
    plt5.clf()
################################################################################
plt6.figure(figsize=(3.5, 3))
T = 800
print(f'temperature = {T} K')
for p, color in zip(p_list, color_list):
    df = read_data(T, p)
    custom_bins = np.linspace(start=1.65, stop=mplt(df['muL']) + 1, num=20)
    counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
    x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
    plt6.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')
    
plt6.tick_params(pltis='x', which='both', labelsize=ticksize,
                direction='in')
plt6.tick_params(pltis='y',which='both',labelsize=ticksize) 
plt6.set_xlim([1.65,3.75])
plt6.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
plt6.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
plt6.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
plt6.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})
plt6.tight_layout()

plt6.savefig(f'../dipole_distribution_{T}.pdf')
plt6.savefig(f'../dipole_distribution_{T}.png')

if (show):
    plt6.show()
else:    
    plt6.clf()
################################################################################
plt7.figure(figsize=(3.5, 3))
T = 900
print(f'temperature = {T} K')
for p, color in zip(p_list, color_list):
    df = read_data(T, p)
    custom_bins = np.linspace(start=1.65, stop=mplt(df['muL']) + 1, num=20) 
    counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
    x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
    plt7.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')

plt7.tick_params(pltis='x', which='both', direction='in',
                labelsize=ticksize, labelbottom=False)
plt7.tick_params(pltis='y',which='both',labelsize=ticksize) 
plt7.set_xlim([1.65,3.75])
plt7.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
plt7.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
plt7.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
plt7.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})
plt7.tight_layout()

plt7.savefig(f'../dipole_distribution_{T}.pdf')
plt7.savefig(f'../dipole_distribution_{T}.png')

if (show):
    plt7.show()
else:    
    plt7.clf()
################################################################################
plt8.figure(figsize=(3.5, 3))
T = 1000
print(f'temperature = {T} K')
for p, color in zip(p_list, color_list):
    df = read_data(T, p)
    custom_bins = np.linspace(start=1.65, stop=mplt(df['muL']) + 1, num=20)
    counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
    x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
    plt8.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')
    
plt8.tick_params(pltis='x', which='both', labelsize=ticksize,
                direction='in')
plt8.tick_params(pltis='y',which='both',labelsize=ticksize) 
plt8.set_xlim([1.65,3.75])
plt8.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
plt8.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
plt8.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
plt8.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})
plt8.tight_layout()

plt8.savefig(f'../dipole_distribution_{T}.pdf')
plt8.savefig(f'../dipole_distribution_{T}.png')

if (show):
    plt8.show()
else:    
    plt8.clf()
################################################################################
