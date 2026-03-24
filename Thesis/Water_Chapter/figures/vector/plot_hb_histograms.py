#!/usr/bin/python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob


################################################################################
################################################################################
################################################################################
def read_data(T, p):
    csvfile_list = glob.glob(f'../../Replica_Dipoles/Replica_*/{T}K/{p:.1f}Bar/hbanalysis.csv')
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
 
fontsize     = 10
ticksize     = 10
legendsize   = 6
alpha        = 0.5
bins         = 20
panel_labels = ['(a)', '(b)', '(c)']
 
T_list = [298, 400, 500, 600, 700, 800, 900, 1000]
p_list = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
 
panels = [
    ('nHB',       r'$N_{\rm HB}$',       0, 5),
    ('donors',    r'$N_{\rm Donors}$',    0, 2),
    ('acceptors', r'$N_{\rm Acceptors}$', 0, 3),
]
 
cmaps_hb = {
    'nHB':       plt.cm.plasma,
    'donors':    plt.cm.viridis,
    'acceptors': plt.cm.cividis,
}
 
################################################################################
################################################################################
################################################################################
for T in T_list:
    print(f'temperature = {T} K')
 
    all_data = {}
    for p in p_list:
        df = read_data(T, p)
        if not df.empty:
            all_data[p] = df
 
    if not all_data:
        print(f'  No data found for T={T} K, skipping.')
        continue
 
    all_muL   = pd.concat([df['muL'] for df in all_data.values()])
    mu_min    = all_muL.min()
    mu_max    = all_muL.max()
    bin_edges = np.linspace(mu_min, mu_max, bins + 1)
 
    fig, axes = plt.subplots(1, 3, figsize=(6.5, 2.8), sharey=False)
 
    for ax, (col, col_label, val_min, val_max), panel_label in zip(axes, panels, panel_labels):
 
        # Only include physically valid values
        all_col = pd.concat(all_data.values())[col].dropna()
        unique_vals = sorted([v for v in all_col.unique() if val_min <= v <= val_max])
        n_vals      = len(unique_vals)
        hb_colors   = [cmaps_hb[col](x) for x in np.linspace(0.15, 0.85, n_vals)]
 
        for val, hb_color in zip(unique_vals, hb_colors):
            subset = pd.concat([
                df[df[col] == val]['muL'] for df in all_data.values()
            ])
            if subset.empty:
                continue
            # Plot as filled histogram bars rather than curves
            ax.hist(
                subset,
                bins      = bin_edges,
                density   = True,
                alpha     = alpha,
                color     = hb_color,
                edgecolor = 'none',
                label     = f'{col_label} $= {int(val)}$',
            )
 
        ax.set_xlabel(r'$\bar{\mu}_L$ / D', fontsize=fontsize)
        ax.set_ylabel(r'$p(\bar{\mu}_L)$ / D$^{-1}$', fontsize=fontsize)
        ax.set_xlim([1.65, 3.75])
        ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
        ax.tick_params(axis='both', which='both', direction='in', labelsize=ticksize)
 
        ax.text(0.02, 0.97, panel_label, transform=ax.transAxes,
                fontsize=fontsize, va='top', ha='left')
 
        ax.legend(loc='upper right', prop={'size': legendsize}, frameon=False)
 
    plt.tight_layout()
    plt.savefig(f'../hb_histogram_{T}.pdf')
    plt.savefig(f'../hb_histogram_{T}.png')
    print(f'  Saved: hb_histogram_{T}.png / .pdf')
 
    if show:
        plt.show()
    else:
        plt.clf()
 

################################################################################
################################################################################
################################################################################
