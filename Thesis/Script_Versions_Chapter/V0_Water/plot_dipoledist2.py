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

class plot_dipoledist2(object):
    def __init__(self):
        pass
################################################################################
################################################################################
################################################################################
    def read_data(self,T, p):
        csvfile_list = glob.glob(f'./Simulations/Replica_*/{T}K/{p:.1f}Bar/Dipole.csv')
        df = pd.DataFrame()
        for csvfile in csvfile_list:
            tmp = pd.read_csv(csvfile)
            df = pd.concat([df, tmp])
        return df


################################################################################
################################################################################
################################################################################
    def plot_dipoledist2(self):
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

        fig, ax = plt.subplots(5, 1, figsize=[4, 10])
        ax1 = ax[0]
        ax2 = ax[1]
        ax3 = ax[2]
        ax4 = ax[3]
        ax5 = ax[4]

################################################################################
        T = 400
        print(f'temperature = {T} K')
        for p, color in zip(p_list, color_list):
            df = read_data(T, p)
    
            custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 0.2, num=20) 
            counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
            x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
            ax1.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')
    
        ax1.text(0.9, 0.9, '(a)', transform=ax1.transAxes, fontsize=fontsize)
        ax1.tick_params(axis='x', which='both', direction='in',
                        labelsize=ticksize, labelbottom=False)
        ax1.tick_params(axis='y',which='both',labelsize=ticksize) 
        ax1.set_xlim([1.65,3.75])
        ax1.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
        ax1.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})
    
        
################################################################################
        T = 600
        print(f'temperature = {T} K')
        for p, color in zip(p_list, color_list):
            df = read_data(T, p)
            custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 1, num=20) 
            counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
            x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
            ax2.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')
    
        ax2.text(0.9, 0.9, '(b)', transform=ax2.transAxes, fontsize=fontsize)
        ax2.tick_params(axis='x', which='both', direction='in',
                        labelsize=ticksize, labelbottom=False)
        ax2.tick_params(axis='y',which='both',labelsize=ticksize) 
        ax2.set_xlim([1.65,3.75])
        ax2.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
    
    
################################################################################
        T = 800
        print(f'temperature = {T} K')
        for p, color in zip(p_list, color_list):
            df = read_data(T, p)
            custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 1, num=20)
            counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
            x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]    
            ax3.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')
        
        ax3.text(0.9, 0.9, '(c)', transform=ax3.transAxes, fontsize=fontsize)
        ax3.tick_params(axis='x', which='both', labelsize=ticksize,
                        direction='in')
        ax3.tick_params(axis='y',which='both',labelsize=ticksize) 
        ax3.set_xlim([1.65,3.75])
        ax3.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
        ax3.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
################################################################################
        T = 900
        print(f'temperature = {T} K')
        for p, color in zip(p_list, color_list):
            df = read_data(T, p)
            custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 1, num=20) 
            counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
            x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
            ax4.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')
    
        ax4.text(0.9, 0.9, '(d)', transform=ax4.transAxes, fontsize=fontsize)
        ax4.tick_params(axis='x', which='both', direction='in',
                        labelsize=ticksize, labelbottom=False)
        ax4.tick_params(axis='y',which='both',labelsize=ticksize) 
        ax4.set_xlim([1.65,3.75])
        ax4.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)


################################################################################
        T = 1000
        print(f'temperature = {T} K')
        for p, color in zip(p_list, color_list):
            df = read_data(T, p)
            custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 1, num=20)
            counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
            x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
            ax5.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')
        
        ax5.text(0.9, 0.9, '(e)', transform=ax5.transAxes, fontsize=fontsize)
        ax5.tick_params(axis='x', which='both', labelsize=ticksize,
                        direction='in')
        ax5.tick_params(axis='y',which='both',labelsize=ticksize) 
        ax5.set_xlim([1.65,3.75])
        ax5.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
        ax5.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
        for i in range(5):
            ax[i].set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
        
################################################################################
        plt.tight_layout()
        fig.subplots_adjust(hspace=0)


        plt.savefig('./Figures/dipole_distribution3.pdf')
        plt.savefig('./Figures/dipole_distribution3.png')

        if (show):
            plt.show()
        else:    
            plt.clf()


################################################################################
################################################################################
