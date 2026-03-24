#!/usr/bin/python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob


################################################################################
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
markersize   = 2
alpha        = 0.6
panel_labels = ['(a)', '(b)', '(c)']
 
T_list = [298, 400, 500, 600, 700, 800, 900, 1000]
p_list = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
 
mu_min, mu_max = 1.65, 3.75
 
panels = [
    ('nHB',       r'$n_{\rm HB}$',       0, 5),
    ('donors',    r'$n_{\rm Donors}$',    0, 2),
    ('acceptors', r'$n_{\rm Acceptors}$', 0, 3),
]
 
cmap       = plt.cm.rainbow
color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(p_list))]
 
################################################################################
################################################################################
################################################################################
for T in T_list:
    print(f'temperature = {T} K')
 
    fig, axes = plt.subplots(1, 3, figsize=(6.5, 2.8))
 
    for ax, (col, col_label, ymin, ymax), panel_label in zip(axes, panels, panel_labels):
 
        for p, colour in zip(p_list, color_list):
            df = read_data(T, p)
 
            if df.empty or 'muL' not in df.columns or col not in df.columns:
                continue
 
            # Filter to physically valid values only
            df = df[(df[col] >= ymin) & (df[col] <= ymax)]
 
            ax.plot(
                df['muL'], df[col],
                ls='', marker='o', markersize=markersize,
                alpha=alpha, color=colour,
                label=f'p={p:.0f} bar',
            )
 
        ax.set_xlim([mu_min, mu_max])
        ax.set_ylim([ymin - 0.5, ymax + 0.5])
        ax.set_xlabel(r'$\mu$ / D', fontsize=fontsize)
        ax.set_ylabel(col_label, fontsize=fontsize)
        ax.set_xticks([1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
        ax.set_yticks(range(ymin, ymax + 1))
        ax.tick_params(axis='both', which='both', direction='in', labelsize=ticksize)
 
        ax.text(0.02, 0.97, panel_label, transform=ax.transAxes,
                fontsize=fontsize, va='top', ha='left')
 
    axes[-1].legend(bbox_to_anchor=(1.05, 1.0), loc='upper left',
                    prop={'size': legendsize}, frameon=False)
 
    plt.tight_layout()
    plt.savefig(f'../dipdist_hb_{T}.pdf', bbox_inches='tight')
    plt.savefig(f'../dipdist_hb_{T}.png', bbox_inches='tight')
    print(f'  Saved: dipdist_hb_{T}.png / .pdf')
 
    if show:
        plt.show()
    else:
        plt.close(fig)

################################################################################
################################################################################
################################################################################
