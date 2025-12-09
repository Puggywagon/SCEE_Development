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
class plot_dipoledist(object):
    def __init__(self):
        pass
################################################################################
################################################################################
################################################################################
    def read_data(self):
        csvfile_list = glob.glob(f'./Results/Replica_*/*K/*Bar/Dipole.csv')
        df = pd.DataFrame()
        for csvfile in csvfile_list:
            tmp = pd.read_csv(csvfile)
            df = pd.concat([df, tmp])
        return df


################################################################################
################################################################################
################################################################################
    def plot_dipoledist(self):
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

        fig, ax = plt.subplots(1, 1, figsize=[20, 8])
        df = read_data() 
        custom_bins = np.linspace(start=1.65, stop=max(df['muL']) + 0.2, num=20) 
        counts,bins = np.histogram(df['muL'], bins=custom_bins, density=True)
        x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
        ax.plot(x_list, counts, color=color, label=f'p={p:.0f} bar')
    
        ax.text(0.8, 0.9, '(a)', transform=ax1.transAxes, fontsize=fontsize)
        ax.tick_params(axis='x', which='both', labelsize=ticksize,
                direction='in')
        ax.tick_params(axis='y',which='both',labelsize=ticksize) 
        ax.set_xlim([1.65,3.75])
        ax.set_ylabel(r'$p(\mu)$ / D$^{-1}$', fontsize=fontsize)
        ax.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
        ax.tick_params(axis='x', which='both', labelsize=ticksize,
            direction='in')
################################################################################
        plt.tight_layout()
        fig.subplots_adjust(hspace=0)


        plt.savefig('./Figures/dipole_distribution.pdf')
        plt.savefig('./Figures/dipole_distribution.png')

        if (show):
            plt.show()
        else:    
            plt.clf()
################################################################################
################################################################################
################################################################################
