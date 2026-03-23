#!/usr/bin/python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob



################################################################################
################################################################################
################################################################################
def get_data(tag, T, p):
    data = np.loadtxt(f"../../TIP4P_Results/5ns/{T:.0f}K/{p:.1f}Bar/{tag}_RDF.xvg",
                      skiprows=26)
    r = data[:, 0]
    gr = data[:, 1]
    return r, gr
 
 
################################################################################
parser = argparse.ArgumentParser()
parser.add_argument('--show', type=str, default='True', help='plot figures')
parser.add_argument('--T', type=float, default=298.0, help='temperature / K')
args = parser.parse_args()
 
show = (args.show == 'True')
 
fontsize = 10
legendsize = 6
ticksize = 10
 
T = args.T
print(f'temperature = {T:.0f} K')
 
################################################################################
p_list = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
 
cmap = plt.cm.rainbow
color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(p_list))]
 
fig, ax = plt.subplots(1, 3, figsize=[10, 4])
 
################################################################################
# Panel (a) — g_OO(r)
i = 0
for p, color in zip(p_list, color_list):
    r, gr = get_data('OO', T, p)
    ax[i].plot(r, gr, color=color, label=f'p={p:.0f} bar')
ax[i].text(0.05, 0.95, '(a)', transform=ax[i].transAxes, fontsize=fontsize,
           va='top', ha='left')
ax[i].tick_params(axis='x', which='both', direction='in',
                  labelsize=ticksize, labelbottom=False)
ax[i].set_xlim([0, 1])
ax[i].set_ylabel('$g_{OO}(r)$', fontsize=fontsize)
 
# Panel (b) — g_OH(r)
i = 1
for p, color in zip(p_list, color_list):
    r, gr = get_data('OH', T, p)
    ax[i].plot(r, gr, color=color, label=f'p={p:.0f} bar')
ax[i].text(0.05, 0.95, '(b)', transform=ax[i].transAxes, fontsize=fontsize,
           va='top', ha='left')
ax[i].tick_params(axis='x', which='both', direction='in',
                  labelsize=ticksize, labelbottom=False)
ax[i].set_ylabel('$g_{OH}(r)$', fontsize=fontsize)
ax[i].set_xlim([0, 1])
 
# Panel (c) — g_HH(r)
i = 2
for p, color in zip(p_list, color_list):
    r, gr = get_data('HH', T, p)
    ax[i].plot(r, gr, color=color, label=f'p={p:.0f} bar')
ax[i].text(0.05, 0.95, '(c)', transform=ax[i].transAxes, fontsize=fontsize,
           va='top', ha='left')
ax[i].tick_params(axis='x', which='both', direction='in',
                  labelsize=ticksize, labelbottom=True)
ax[i].set_xlabel('r / nm', fontsize=fontsize)
ax[i].set_ylabel('$g_{HH}(r)$', fontsize=fontsize)
ax[i].set_xlim([0, 1])
 
# Legend outside to the right of panel (c)
ax[2].legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize},
             frameon=True, framealpha=0.8, borderaxespad=0)
 

################################################################################
plt.tight_layout()
fig.subplots_adjust(hspace=0)

plt.savefig(f'../rdf_{T:.0f}.pdf', dpi=300)
plt.savefig(f'../rdf_{T:.0f}.png', dpi=300)
if (show):
    plt.show()
else:    
    plt.clf()


################################################################################
################################################################################
################################################################################
