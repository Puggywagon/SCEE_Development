import gromacs
import subprocess
import numpy as np
import os
import pandas as pd
import math
import matplotlib.pyplot as plt
import glob
######################################################
def read_data():
    df3 = pd.read_csv('Dipoles2.csv')#,header=False)
    df3.columns = ['Temperature', 'Pressure', 'Dipoles','SE','std means']
    return df3
##########################################################################################
def plt_dipole(df3,T_list):    
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
    fig,ax1= plt.subplots(figsize=(15, 12))
    ax2=ax1.twinx()
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    for T, color in zip(T_list, color_list):
        tmp = df3[df3['Temperature'] == T]
        minval=tmp['SE']
        maxval=tmp['SE']
        corrections=tmp['Dipoles']/2.3054
        yerr=[minval,maxval]
        ax1.errorbar(tmp['Pressure'], tmp['Dipoles'], yerr=yerr, capsize=3, ecolor = "black", ls='none', color=color, marker='o')
        #ax2.plot(tmp['Pressure'], corrections, ls='none', color=color, marker='s')
        plt.plot([], [], ls='', color=color, marker='o',label=r'$T='+f'{T}$ K')
        a,b=np.polyfit(tmp['Pressure'],tmp['Dipoles'],1)
        x=tmp['Pressure']
        ax1.plot(x,(a*x+b), ls='--',color=color)
    ax1.set_ylabel(r'Dipole Moment / D',fontsize=22)
    ax2.set_ylabel(r'Scaling Factor',fontsize=22)
    ax1.set_xlabel(r'Pressure / bar',fontsize=22)
    ax1.tick_params(axis='x', which='both', labelsize=18)
    ax1.tick_params(axis='y',which='both',labelsize=18) 
    ax2.tick_params(axis='y',which='both',labelsize=18) 
    plt.xscale('log')
    plt.xlim([0.5, 1.0e4])
    mn,mx=ax1.set_ylim([1.6, 3.0])
    min_correction=mn/2.3054
    max_correction=mx/2.3054
    ax2.set_ylim([min_correction, max_correction])
    #plt.legend(loc='center left', bbox_to_anchor=(1.04, 0.5))
    Min=min(T_list)
    Max=max(T_list)
    plt.savefig(f'Dipole_Corrections_{Min}_{Max}.pdf', dpi=300)
    plt.savefig(f'Dipole_Corrections_{Min}_{Max}.png', dpi=300)
##########################################################################################
def plt_stddipole(df3,T_list):    
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
    plt.figure(figsize=(15, 12))
    for T, color in zip(T_list, color_list):
        tmp = df3[df3['Temperature'] == T]
        plt.plot(tmp['Pressure'],  tmp['std means'], ls='none', color=color, marker='o')
        plt.plot([], [], ls='', color=color, marker='o',label=r'$T='+f'{T}$ K')
        a,b=np.polyfit(tmp['Pressure'],tmp['std means'],1)
        x=tmp['Pressure']
        plt.plot(x,(a*x+b), ls='--',color=color)
    plt.ylabel(r'Dipole Moment Distribution Width / D',fontsize=22)
    plt.xlabel(r'Pressure / bar',fontsize=22)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18) 
    plt.xscale('log')
    plt.xlim([0.5, 1.0e4])
    plt.ylim([0, 0.35])
    #plt.legend(loc='center left', bbox_to_anchor=(1.04, 0.5))
    Min=min(T_list)
    Max=max(T_list)
    plt.savefig(f'Dipole_Corrections_Stdev_{Min}_{Max}.pdf', dpi=300)
    plt.savefig(f'Dipole_Corrections_Stdev_{Min}_{Max}.png', dpi=300)
###############################################################################################
def plt_close():
    plt.close()
##########################################################################################
def plotting_dipole():
    T_list=[298]
    for temp in range(400,1001,100):
        T_list.append(temp)
    p_ref = [1,2,5]
    P_list = []
    for scale in np.logspace(0, 3, 4):
        P_list += [p*scale for p in p_ref]
    df3=read_data()
    plt_dipole(df3,T_list)
    plt_close()
    plt_stddipole(df3,T_list)
    plt_close()
##########################################################################################   
plotting_dipole()
