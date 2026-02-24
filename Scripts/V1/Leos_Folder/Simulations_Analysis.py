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
    def write_filtered_data(self,k):
        csvfile_list = glob.glob(f'Simulations/replica_*/*/*/Dipole_Calculations.csv')
        for csvfile in csvfile_list:
            df = pd.read_csv(csvfile)
            x = df['mu_liquid'].astype(float).to_numpy()
            finite_mask = np.isfinite(x)
            x_finite = x[finite_mask]
            if x_finite.size < 5:
                continue
            keep, out, _ = self._iqr_masks(x_finite, k=k)
            full_keep = np.zeros(len(x), dtype=bool)
            full_keep[finite_mask] = keep
            df_filtered = df.copy()
            cols_to_blank = ['muL_SCEE', 'delta_mu', 'mu_liquid']
            df_filtered.loc[~full_keep, cols_to_blank] = np.nan
            out_dir = os.path.dirname(csvfile)
            out_path = os.path.join(out_dir, 'Dipole_Calculations_filtered.csv')
            df_filtered.to_csv(out_path, index=False)
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
    def plot_dipoledist(self,Molecule, bin_width=0.25, x_ceil=10.0):
        parser = argparse.ArgumentParser()
        parser.add_argument('--show', type=str, default='True', help='plot figures')
        args = parser.parse_args()
        show = (args.show == 'True')

        df = self.read_data()
        x_all = df['mu_liquid'].astype(float).to_numpy()
        x_all = x_all[np.isfinite(x_all)]
        if x_all.size < 5:
            raise ValueError("Not enough dipole samples to form a distribution.")

        keep, out, (q1, q3, iqr, lo, hi) = self._iqr_masks(x_all, k=1.5)
        x_keep = x_all[keep]
        x_out  = x_all[out]

        mu    = x_keep.mean()
        sigma = x_keep.std(ddof=0)
        if sigma == 0.0:
            raise ValueError("Sigma is zero after outlier removal. Check your data.")

        # Axis limits
        xmin = 0.0
        largest_populated = x_all.max()  # will be overridden below
        # find the largest bin edge that actually has data within x_ceil
        x_within_ceil = x_all[x_all <= x_ceil]
        if x_within_ceil.size > 0:
            largest_populated = np.ceil(x_within_ceil.max() / bin_width) * bin_width
        xmax = min(largest_populated, x_ceil)
        xmax = np.ceil(xmax / bin_width) * bin_width  # snap to bin edge

        bin_edges = np.arange(xmin, xmax + bin_width, bin_width)
        bw = bin_edges[1] - bin_edges[0]

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.hist(x_keep, bins=bin_edges, alpha=0.9, label="Data (outliers removed)")
        if x_out.size:
            # only plot outliers that fall within the axis range
            x_out_visible = x_out[(x_out >= xmin) & (x_out <= xmax)]
            if x_out_visible.size:
                ax.hist(x_out_visible, bins=bin_edges, alpha=0.9, label="Outliers (IQR rule)")

        x_line = np.linspace(xmin, xmax, 400)
        y_line = self._normal_pdf(x_line, mu, sigma) * x_keep.size * bw
        ax.plot(x_line, y_line, linewidth=2, label="Gaussian fit")

        ax.set_xlim(xmin, xmax)
        tick_spacing = 0.5 if (xmax - xmin) < 5 else 1.0
        ax.set_xticks(np.arange(xmin, xmax + tick_spacing, tick_spacing))
        ax.set_xlabel(r'$\mu_{liquid}$ / D', fontsize=16)
        ax.set_ylabel('Frequency', fontsize=16)
        ax.legend(frameon=False, fontsize=12)
        ax.tick_params(axis='both', which='major', labelsize=12)
        plt.tight_layout()
        plt.savefig(f'../figures/dipole_distribution_{Molecule}.png', bbox_inches="tight", dpi=300)
        plt.savefig(f'../figures/dipole_distribution_{Molecule}.pdf', bbox_inches="tight", dpi=300)
        plt.close(fig)

        # Summary statistics table
        stats = {
            'Molecule':         [Molecule],
            'N_total':          [x_all.size],
            'N_kept':           [x_keep.size],
            'N_outliers':       [x_out.size],
            'mean_D':           [round(mu, 4)],
            'std_D':            [round(sigma, 4)],
            'median_D':         [round(np.median(x_keep), 4)],
            'IQR_lo':           [round(lo, 4)],
            'IQR_hi':           [round(hi, 4)],
            'q1':               [round(q1, 4)],
            'q3':               [round(q3, 4)],
        }
        stats_df = pd.DataFrame(stats)
        stats_df.to_csv(f'../vectors/figures/dipole_stats_{Molecule}.csv', index=False)
        print(stats_df.to_string(index=False))
######################################################
    def read_filtered_data(v):
        csvfile_list = glob.glob(f'Simulations/replica_*/*/*/Dipole_Calculations_filtered.csv')
        df = pd.DataFrame()
        for csvfile in csvfile_list:
            tmp = pd.read_csv(csvfile)
            df = pd.concat([df, tmp])
        return df
######################################################
    def plot_dipoledist_filtered(self,Molecule, bin_width=0.25, x_ceil=10.0):
        parser = argparse.ArgumentParser()
        parser.add_argument('--show', type=str, default='True', help='plot figures')
        args = parser.parse_args()
        show = (args.show == 'True')

        df = self.read_filtered_data()
        x_all = df['mu_liquid'].astype(float).to_numpy()
        x_all = x_all[np.isfinite(x_all)]
        if x_all.size < 5:
            raise ValueError("Not enough dipole samples to form a distribution.")

        mu    = x_all.mean()
        sigma = x_all.std(ddof=0)
        if sigma == 0.0:
            raise ValueError("Sigma is zero. Check your data.")

        # Axis limits
        xmin = 0.0
        x_within_ceil = x_all[x_all <= x_ceil]
        if x_within_ceil.size > 0:
            largest_populated = np.ceil(x_within_ceil.max() / bin_width) * bin_width
        xmax = min(largest_populated, x_ceil)
        xmax = np.ceil(xmax / bin_width) * bin_width

        bin_edges = np.arange(xmin, xmax + bin_width, bin_width)
        bw = bin_edges[1] - bin_edges[0]
    
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.hist(x_all, bins=bin_edges, alpha=0.9, label="Data (outliers removed)")
    
        x_line = np.linspace(xmin, xmax, 400)
        y_line = self._normal_pdf(x_line, mu, sigma) * x_all.size * bw
        ax.plot(x_line, y_line, linewidth=2, label="Gaussian fit")

        ax.set_xlim(xmin, xmax)
        tick_spacing = 0.5 if (xmax - xmin) < 5 else 1.0
        ax.set_xticks(np.arange(xmin, xmax + tick_spacing, tick_spacing))
        ax.set_xlabel(r'$\mu_{liquid}$ / D', fontsize=16)
        ax.set_ylabel('Frequency', fontsize=16)
        ax.legend(frameon=False, fontsize=12)
        ax.tick_params(axis='both', which='major', labelsize=12)
        plt.tight_layout()
        plt.savefig(f'../figures/dipole_distribution_{Molecule}_filtered.png', bbox_inches="tight", dpi=300)
        plt.savefig(f'../figures/dipole_distribution_{Molecule}_filtered.pdf', bbox_inches="tight", dpi=300)
        plt.close(fig)

        # Summary statistics table
        stats = {
            'Molecule':   [Molecule],
            'N_total':    [x_all.size],
            'mean_D':     [round(mu, 4)],
            'std_D':      [round(sigma, 4)],
            'median_D':   [round(np.median(x_all), 4)],
        }
        stats_df = pd.DataFrame(stats)
        stats_df.to_csv(f'../figures/vectors/dipole_stats_{Molecule}_filtered.csv', index=False)
        print(stats_df.to_string(index=False))
############################################################################################################
