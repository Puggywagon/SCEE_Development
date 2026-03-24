#!/usr/bin/python3

import argparse
import gromacs
import subprocess
import numpy as np
import os
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


################################################################################
################################################################################
################################################################################
def plot_experiment():
    df = pd.read_csv(nist_vp_csvfile)
#    df.columns = ['Temperature (K)', 'Pressure (bar)', 'Phase']
    Tc = 647.096      # Critical temperature (Tc) 647.096 K
    pc = 220.640      # Critical pressure (Pc)    220.640 bar
    rhoc = 17.873728  # Critical density (Dc)     17.873728 mol/l

    tmp = df[df['Phase']=='liquid']
    x_list = list(tmp['Temperature (K)']) + [Tc]
    y_list = list(tmp['Pressure (bar)']) + [pc]
    plt.plot(x_list, y_list, ls='solid', color='black')
    plt.plot([Tc], [pc], ls='-',marker='o', color='black',label='Experiment')

    plt.xlim([250, 1101])
    plt.ylim([0.5, 1.0e4])
    plt.xticks([300,500,700,1000])
    plt.yscale('log')
    plt.ylabel('pressure / bar')
    plt.xlabel('temperature / K')

    
################################################################################
def plot_vega():
    df1 = pd.read_csv(tip4p_vp_csvfile)
#    df1.columns = ['T / K', 'p / bar', 'rhol','rhog']
      
    Tc = 640
    pc = 146
    rhoc = 0.31   # g cm^{-3} 
    
    x_list =  [Tc] + list(df1['T / K'])
    y_list =  [pc] + list(df1['p / bar'])
    plt.plot(x_list, y_list, color='black', ls='dashed')
    plt.plot([Tc], [pc], ls='--',marker='o', color='black',markerfacecolor='none',
         markeredgecolor='k',label='TIP4P/2005')
    
    plt.xlim([250, 1101])
    plt.ylim([0.5, 1.0e4])
    plt.yscale('log')

    
################################################################################
################################################################################
################################################################################
parser = argparse.ArgumentParser()

parser.add_argument('--show', type=str, default='True', help='plot figures')

args = parser.parse_args()

show = (args.show == 'True')

fontsize = 10
legendsize = 6

nist_vp_csvfile = './NIST_vaporpressure.csv'
tip4p_vp_csvfile = './TIP4P_vaporpressure.csv'


################################################################################
def liquid_boundary(T,p):
    if T==298 and p==1.0:
        return "liquid"
    elif T==298 and p==2.0:
        return "liquid"
    elif T==298 and p==5.0:
        return "liquid"
    elif T==298 and p==10.0:
        return "liquid"
    elif T==298 and p==20.0:
        return "liquid"
    elif T==298 and p==50.0:
        return "liquid"
    elif T==298 and p==100.0:
        return "liquid"
    elif T==298 and p==200.0:
        return "liquid"
    elif T==298 and p==500.0:
        return "liquid"
    elif T==298 and p==1000.0:
        return "liquid"
    elif T==298 and p==2000.0:
        return "liquid"
    elif T==298 and p==5000.0:
        return "liquid"
    elif T==400 and p==1.0:
        return "liquid"
    elif T==400 and p==2.0:
        return "liquid"
    elif T==400 and p==5.0:
        return "liquid"
    elif T==400 and p==10.0:
        return "liquid"
    elif T==400 and p==20.0:
        return "liquid"
    elif T==400 and p==50.0:
        return "liquid"
    elif T==400 and p==100.0:
        return "liquid"
    elif T==400 and p==200.0:
        return "liquid"
    elif T==400 and p==500.0:
        return "liquid"
    elif T==400 and p==1000.0:
        return "liquid"
    elif T==400 and p==2000.0:
        return "liquid"
    elif T==400 and p==5000.0:
        return "liquid"
    elif T==500 and p==20.0:
        return "liquid"
    elif T==500 and p==50.0:
        return "liquid"
    elif T==500 and p==100.0:
        return "liquid"
    elif T==500 and p==200.0:
        return "liquid"
    elif T==500 and p==500.0:
        return "liquid"
    elif T==500 and p==1000.0:
        return "liquid"
    elif T==500 and p==2000.0:
        return "liquid"
    elif T==500 and p==5000.0:
        return "liquid"
    elif T==600 and p==100.0:
        return "liquid"
    elif T==600 and p==200.0:
        return "liquid"
    elif T==600 and p==500.0:
        return "liquid"
    elif T==600 and p==1000.0:
        return "liquid"
    elif T==600 and p==2000.0:
        return "liquid"
    elif T==600 and p==5000.0:
        return "liquid"
    else:
        return 0
################################################################################
def gas_boundary(T,p):
    if T==500 and p==1.0:
        return "gas"
    elif T==500 and p==2.0:
        return "gas"
    elif T==500 and p==5.0:
        return "gas"
    elif T==500 and p==10.0:
        return "gas"
    elif T==600 and p==1.0:
        return "gas"
    elif T==600 and p==2.0:
        return "gas"
    elif T==600 and p==5.0:
        return "gas"
    elif T==600 and p==10.0:
        return "gas"
    elif T==600 and p==20.0:
        return "gas"
    elif T==600 and p==50.0:
        return "gas"
    elif T==700 and p==1.0:
        return "gas"
    elif T==700 and p==2.0:
        return "gas"
    elif T==700 and p==5.0:
        return "gas"
    elif T==700 and p==10.0:
        return "gas"
    elif T==700 and p==20.0:
        return "gas"
    elif T==700 and p==50.0:
        return "gas"
    elif T==700 and p==100.0:
        return "gas"
    elif T==800 and p==1.0:
        return "gas"
    elif T==800 and p==2.0:
        return "gas"
    elif T==800 and p==5.0:
        return "gas"
    elif T==800 and p==10.0:
        return "gas"
    elif T==800 and p==20.0:
        return "gas"
    elif T==800 and p==50.0:
        return "gas"
    elif T==800 and p==100.0:
        return "gas"
    elif T==900 and p==1.0:
        return "gas"
    elif T==900 and p==2.0:
        return "gas"
    elif T==900 and p==5.0:
        return "gas"
    elif T==900 and p==10.0:
        return "gas"
    elif T==900 and p==20.0:
        return "gas"
    elif T==900 and p==50.0:
        return "gas"
    elif T==900 and p==100.0:
        return "gas"
    elif T==1000 and p==1.0:
        return "gas"
    elif T==1000 and p==2.0:
        return "gas"
    elif T==1000 and p==5.0:
        return "gas"
    elif T==1000 and p==10.0:
        return "gas"
    elif T==1000 and p==20.0:
        return "gas"
    elif T==1000 and p==50.0:
        return "gas"
    elif T==1000 and p==100.0:
        return "gas"
    else:
        return 0

    
####################################################################################
fontsize = 10
legendsize = 6

fig, ax = plt.subplots(figsize=(5.5, 3.5))

plot_experiment()
plot_vega()

label_liquid = 0
label_gas = 0
label_sc = 0
#time()
T_list = [298, 400, 500, 600, 700, 800, 900, 1000]
p_list = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]

for T in T_list:
    for p in p_list:
        if liquid_boundary(T, p) == "liquid":
            if label_liquid==0:
                plt.plot(T, p, ls='None', color='blue', marker='o',
                     label='Liquid')
            else:
                plt.plot(T, p, ls='None', color='blue', marker='o',
                     label=None)
            label_liquid += 1

        elif gas_boundary(T, p) == "gas":
            if label_gas==0:
                plt.plot(T, p, ls='None', color='green', marker='^',
                     label='Gas')
            else:
                plt.plot(T, p, ls='None', color='green', marker='^',
                     label=None)
            label_gas += 1

        else:
            if label_sc==0:
                plt.plot(T, p, ls='None', color='red', marker='*',
                     label='Supercritical')
            else:
                plt.plot(T, p, ls='None', color='red', marker='*',
                     label=None)
            label_sc += 1
################################################################################
plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', prop={'size': legendsize})  
plt.tight_layout()
fig.subplots_adjust(hspace=0)

plt.savefig(f'../states.pdf', dpi=300)
plt.savefig(f'../states.png', dpi=300)
if (show):
    plt.show()
else:    
    plt.clf()


################################################################################
################################################################################
################################################################################


