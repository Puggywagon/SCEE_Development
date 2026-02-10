#!/usr/bin/python3

import subprocess
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt


################################################################################
################################################################################
################################################################################
def plot_densitydipole(df2,T_list):
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
    plt.figure(figsize=(15, 12))
    for T, color in zip(T_list, color_list):
        tmp = df2[df2['Temperature'] == T]
        plt.plot(tmp['Density'], tmp['Dipole Moment'], ls='None',
                 color=color, marker='o',label=r'$T='+f'{T}$ K')

    plt.ylabel(r'Dipole Moment / D',fontsize=22)
    plt.ylim([1.6, 3.0])
    plt.xlabel(r'Density / kg m$^{-3}$',fontsize=22)
    plt.xlim([0, 1200])
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)  
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    Min=min(T_list)
    Max=max(T_list)
    plt.savefig(f'density_dipole_{Min}_{Max}.pdf', dpi=300)
    plt.savefig(f'density_dipole_{Min}_{Max}.png', dpi=300)
    plt.close()   
############################################################################
def analysis():
    T_list1=[298]
    for temp in range(400,1001,100):
        T_list1.append(temp)
    T_list2=[298]
    for temp in range(273,374,10):
        T_list2.append(temp)
    for T_list in [T_list1,T_list2]:
        enter_analysis()
        df1=read_nist()
        df2=read_data()
        plot_density(df1,df2,T_list)
        plot_epsilondipole(df2,T_list)
        plot_densitydipole(df2,T_list)
        exit_dir()
################################################################
def time():
    time_length = [5,10]
    for Time in time_length:
        enter_dirtime(Time)
        analysis()
        exit_dir()
################################################################



T_list = [298, 400, 500, 600, 700, 800, 900, 1000]

#df_nist = pd.read_csv('./NIST_isotherms.csv')
df = pd.read_csv('./5ns/Analysis/data.csv', header=None)
df.columns = ['Time', 'Temperature', 'Pressure', 'Density','Epsilon','Number of Honds','Percentage of H-bonds','Nearest Neighbour','Dipole Moment','State']
print(df_nist.columns)


cmap = plt.cm.rainbow
color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]

