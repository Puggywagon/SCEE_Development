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
    file_list = glob.glob("../../TIP4P_Results/5ns/{T:.0f}K/{p:.1f}Bar/{tag}_RDF.xvg")
#    file_list += glob.glob("../../../TIP4P_Backup/Replicas/runs_leo/replica_*/5ns/{T:.0f}K/{p:.1f}Bar/{tag}_RDF.xvg")
    data = np.loadtxt(f"../../TIP4P_Results/5ns/{T:.0f}K/{p:.1f}Bar/{tag}_RDF.xvg",
                      skiprows=26) 
    r = data[:, 0]
    gr = data[:, 1]
    return r, gr


################################################################################
################################################################################
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
################################################################################
################################################################################
dirlist = f'./5ns/'
T_list = [298, 400, 500, 600, 700, 800, 900, 1000]
p_list = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]

cmap = plt.cm.rainbow
color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(p_list))]

fig, ax = plt.subplots(3, 1, figsize=[4, 10])


i = 0
for p, color in zip(p_list, color_list):
    r, gr = get_data('OO', T, p)
    ax[0].plot(r,gr, color=color, label=f'p={p:.0f} bar')
ax[i].text(0.9, 0.9, '(a)', transform=ax[i].transAxes, fontsize=fontsize)
ax[i].tick_params(axis='x', which='both', direction='in',
                  labelsize=ticksize, labelbottom=False)
ax[i].set_xlim([0, 1])
ax[i].set_ylabel('$g_{OO}(r)$',fontsize=fontsize)
ax[i].legend(loc='upper center', prop={'size': legendsize})


i = 1
for p, color in zip(p_list, color_list):
    r, gr = get_data('OH', T, p)
    ax[i].plot(r,gr, color=color, label=f'p={p:.0f} bar')
ax[i].text(0.9, 0.9, '(b)', transform=ax[i].transAxes, fontsize=fontsize)
ax[i].tick_params(axis='x', which='both', direction='in',
                  labelsize=ticksize, labelbottom=False)
ax[i].set_ylabel('$g_{OH}(r)$',fontsize=fontsize)
ax[i].set_xlim([0, 1])

i = 2
for p, color in zip(p_list, color_list):
    r, gr = get_data('HH', T, p)
    ax[i].plot(r,gr, color=color, label=f'p={p:.0f} bar')
ax[i].text(0.9, 0.9, '(c)', transform=ax[i].transAxes, fontsize=fontsize)
ax[i].tick_params(axis='x', which='both', direction='in',
                  labelsize=ticksize, labelbottom=True)
ax[i].set_xlabel('r / nm',fontsize=fontsize)
ax[i].set_ylabel('$g_{HH}(r)$',fontsize=fontsize)
ax[i].set_xlim([0, 1])


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
