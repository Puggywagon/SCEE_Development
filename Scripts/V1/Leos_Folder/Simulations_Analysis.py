#!/usr/bin/python3
import gromacs
import subprocess
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import socket
import time
import argparse
################################################################################
class Simulations_Analysis(object):
    def __init__(self):
        pass
################################################################################
    def get_dipole_model(self): # Think we move this one and the next one into a different section in the script. Maybe we combine the MD and SCEE loops so that it is all happening in one replica folder?
        gromacs.environment.flags['capture_output'] = True
        input_str = '0\n'
        tmp = gromacs.tools.Dipoles(f=f'nvt_vacuum2',
                                s=f'nvt_vacuum2',
                                input=input_str)
        chk, dipole_output, stderr = tmp.run()

        f = open(f'Dipole.txt', 'w')
        f.write(dipole_output)
        f.close()

        dipole_lines = dipole_output.splitlines()  # Split text into lines
        dipole_line = dipole_lines[8]  # Assuming epsilon is on the last line
        dipole_model = float(dipole_line.split()[2])  # Extract the numerical value
        return dipole_model
######################################################
    def get_dipole_model_liquid(self): # Add in 
        gromacs.environment.flags['capture_output'] = True
        input_str = '0\n'
        tmp = gromacs.tools.Dipoles(f=f'Pure_QMMM_md3.trr',
                                s=f'Pure_QMMM_md3.tpr',
                                input=input_str)
        chk, dipole_output, stderr = tmp.run()

        f = open(f'Dipole2.txt', 'w')
        f.write(dipole_output)
        f.close()

        dipole_lines = dipole_output.splitlines()  # Split text into lines
        dipole_line = dipole_lines[8]  # Assuming epsilon is on the last line
        Model_Dipole = float(dipole_line.split()[2])  # Extract the numerical value
        epsilon_line = dipole_lines[-1]  # Assuming epsilon is on the last line
        epsilon = float(epsilon_line.split()[2])  # Extract the numerical value
        
        return Model_Dipole, Epsilon
################################################################################
    def read_data(self):
        csvfile_list = glob.glob(f'Simulations/replica_*/*/*/Dipole_Calculations.csv')
        df = pd.DataFrame()
        for i, csvfile in enumerate(csvfile_list, start=1):
            tmp = pd.read_csv(csvfile)
            df = pd.concat([df, tmp])
        return df
######################################################
    def _iqr_masks(self,x, k=1.5):
        x = np.asarray(x, dtype=float)
        q1, q3 = np.percentile(x, [25, 75])
        iqr = q3 - q1
        lo = q1 - k * iqr
        hi = q3 + k * iqr
        out = (x < lo) | (x > hi)
        keep = ~out
        return keep, out, (q1, q3, iqr, lo, hi)
######################################################
    def _normal_pdf(self,x, mu, sigma):
        return (1.0 / (sigma * np.sqrt(2.0 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
######################################################
    def plot_dipoledist(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--show', type=str, default='True', help='plot figures')
        parser.add_argument('--bins', type=int, default=30, help='number of histogram bins')
        args = parser.parse_args()
        show = (args.show == 'True')
        df = self.read_data()
        x_all = df['mu_liquid'].astype(float).to_numpy()
        x_all = x_all[np.isfinite(x_all)]
        if x_all.size < 5:
            raise ValueError("Not enough dipole samples to form a distribution.")
    
        keep, out, (q1, q3, iqr, lo, hi) = self._iqr_masks(x_all, k=1.5)
        x_keep = x_all[keep]
        x_out = x_all[out]
        
        mu = x_keep.mean()
        sigma = x_keep.std(ddof=0)
        if sigma == 0.0:
            raise ValueError("Sigma is zero after outlier removal. Check your data.")
        fig, ax = plt.subplots(figsize=(8, 6))
        # Binning: match your previous style or use full range
        # This matches your old approach starting at 1.00 D but updates the upper limit safely
        xmin = 1.00
        xmax = max(x_all.max(), 1.00) + 0.2
        bin_edges = np.linspace(xmin, xmax, args.bins + 1)
        bin_width = bin_edges[1] - bin_edges[0]
        # Histograms as counts (frequency)
        ax.hist(x_keep, bins=bin_edges, alpha=0.9, label="Data (outliers removed)")
        if x_out.size:
            ax.hist(x_out, bins=bin_edges, alpha=0.9, label="Outliers (IQR rule)")
        # Gaussian overlay scaled to counts
        x_line = np.linspace(bin_edges[0], bin_edges[-1], 400)
        y_line = self._normal_pdf(x_line, mu, sigma) * x_keep.size * bin_width
        ax.plot(x_line, y_line, linewidth=2, label="Gaussian fit")
    
        ax.set_xlabel(r'$\mu_{liquid}$ / D', fontsize=16)
        ax.set_ylabel('Frequency', fontsize=16)
    
        ax.legend(frameon=False, fontsize=12)
        ax.tick_params(axis='both', which='major', labelsize=12)

        plt.tight_layout()
        plt.savefig(f'dipole_distribution.png', bbox_inches="tight", dpi=300)    
        plt.savefig(f'dipole_distribution.pdf', bbox_inches="tight", dpi=300)
    
        plt.close(fig)  
############################################################################################################
