import subprocess
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
################################################################################
def enter_dirtime(Time):
    dire='{}ns'.format(Time)
    os.chdir(dire)
################################################################################
def enter_analysis():
    dire='Analysis'
    os.chdir(dire)
################################################################################
def exit_dir():
    dire='../'
    os.chdir(dire)
#################################################################################
def read_nist():
    df1 = pd.read_csv('../../nist_data.csv')
    df1.columns = ['Temperature', 'Pressure', 'Density']
    return df1
#################################################################################
def read_data():
    df2 = pd.read_csv('data.csv')
    df2.columns = ['Time', 'Temperature', 'Pressure', 'Density','Epsilon','Number of Honds','Percentage of H-bonds','Nearest Neighbour','Dipole Moment']
    return df2
#################################################################################
def plot_density(df1,df2,T_list):
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
    plt.figure(figsize=(14, 6))
    for T, color in zip(T_list, color_list):
        tmp = df1[df1['Temperature'] == T]
        plt.plot(tmp['Pressure'], tmp['Density'],
                 color=color)
    for T, color in zip(T_list, color_list):
        tmp = df2[df2['Temperature'] == T]
        plt.plot(tmp['Pressure'], tmp['Density'], ls='None',
                 color=color, marker='o',label=r'$T='+f'{T}$ K')

    plt.ylabel(r'Density / kg m$^{-3}$',fontsize=12)
    plt.xlabel(r'Pressure / bar',fontsize=12)
    plt.xscale('log')
    plt.xlim([0.5, 1.0e4])
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title('Density vs Pressure',fontsize=14,weight='bold')
    Min=min(T_list)
    Max=max(T_list)
    plt.savefig(f'density_{Min}_{Max}.pdf', dpi=300)
    plt.savefig(f'density_{Min}_{Max}.png', dpi=300)
    plt.close()
################################################################################
def enter_dirtime(Time):
    dire='{}ns'.format(Time)
    os.chdir(dire)
#################################################################################
def plot_epsilon(df2,T_list):
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
    plt.figure(figsize=(14, 6))
    for T, color in zip(T_list, color_list):
        tmp = df2[df2['Temperature'] == T]
        plt.plot(tmp['Pressure'], tmp['Epsilon'], ls='None',
                 color=color, marker='o',label=r'$T='+f'{T}$ K')

    plt.ylabel(r'Epsilon',fontsize=12)
    plt.xlabel(r'Pressure / bar',fontsize=12)
    plt.xscale('log')
    plt.xlim([0.5, 1.0e4])
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title('Deielectric Constant vs Pressure',fontsize=14,weight='bold')
    Min=min(T_list)
    Max=max(T_list)
    plt.savefig(f'Epsilon_{Min}_{Max}.pdf', dpi=300)
    plt.savefig(f'Epsilon_{Min}_{Max}.png', dpi=300)
    plt.close()
#################################################################################
def plot_nearestneighbour(df2,T_list):
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
    plt.figure(figsize=(14, 6))
    for T, color in zip(T_list, color_list):
        tmp = df2[df2['Temperature'] == T]
        plt.plot(tmp['Pressure'], tmp['Nearest Neighbour'], ls='None',
                 color=color, marker='o',label=r'$T='+f'{T}$ K')

    plt.ylabel(r'Nearest Neighbour',fontsize=12)
    plt.xlabel(r'Pressure / bar',fontsize=12)
    plt.xscale('log')
    plt.xlim([0.5, 1.0e4])
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title('Nearest Neighbour vs Pressure',fontsize=14,weight='bold')
    Min=min(T_list)
    Max=max(T_list)
    plt.savefig(f'Nearest_Neighbour_{Min}_{Max}.pdf', dpi=300)
    plt.savefig(f'Nearest_Neighbour_{Min}_{Max}.png', dpi=300)
    plt.close()
#################################################################################
def plot_dipole(df2,T_list):
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
    plt.figure(figsize=(14, 6))
    for T, color in zip(T_list, color_list):
        tmp = df2[df2['Temperature'] == T]
        plt.plot(tmp['Pressure'], tmp['Dipole Moment'], ls='None',
                 color=color, marker='o',label=r'$T='+f'{T}$ K')

    plt.ylabel(r'Dipole Moment / D',fontsize=12)
    plt.xlabel(r'Pressure / bar',fontsize=12)
    plt.xscale('log')
    plt.xlim([0.5, 1.0e4])
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title('Dipole Moment vs Pressure',fontsize=14,weight='bold')
    Min=min(T_list)
    Max=max(T_list)
    plt.savefig(f'Dipole_Moment_{Min}_{Max}.pdf', dpi=300)
    plt.savefig(f'Dipole_Moment_{Min}_{Max}.png', dpi=300)
    plt.close()
#################################################################################
def plot_densitydipole(df2,T_list):
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
    plt.figure(figsize=(14, 6))
    for T, color in zip(T_list, color_list):
        tmp = df2[df2['Temperature'] == T]
        plt.plot(tmp['Density'], tmp['Dipole Moment'], ls='None',
                 color=color, marker='o',label=r'$T='+f'{T}$ K')

    plt.ylabel(r'Dipole Moment / D',fontsize=12)
    plt.xlabel(r'Density / kg m$^{-3}$',fontsize=12)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title('Dipole Moment vs Density',fontsize=14,weight='bold')
    Min=min(T_list)
    Max=max(T_list)
    plt.savefig(f'density_dipole_{Min}_{Max}.pdf', dpi=300)
    plt.savefig(f'density_dipole_{Min}_{Max}.png', dpi=300)
    plt.close()

#################################################################################
def plot_nearestneighbourdipole(df2,T_list):
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
    plt.figure(figsize=(14, 6))
    for T, color in zip(T_list, color_list):
        tmp = df2[df2['Temperature'] == T]
        plt.plot(tmp['Nearest Neighbour'], tmp['Dipole Moment'], ls='None',
                 color=color, marker='o',label=r'$T='+f'{T}$ K')

    plt.ylabel(r'Dipole Moment / D',fontsize=12)
    plt.xlabel(r'Nearest Neighbour',fontsize=12)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title('Dipole Moment vs Nearest Neighbour',fontsize=14,weight='bold')
    Min=min(T_list)
    Max=max(T_list)
    plt.savefig(f'Nearest_Neighbour_Dipole_{Min}_{Max}.pdf', dpi=300)
    plt.savefig(f'Nearest_Neighbour_Dipole_{Min}_{Max}.png', dpi=300)
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
        plot_epsilon(df2,T_list)
        plot_nearestneighbour(df2,T_list)
        plot_dipole(df2,T_list)
        plot_densitydipole(df2,T_list)
        plot_nearestneighbourdipole(df2,T_list)
        
        exit_dir()
################################################################
def time():
    time_length = [5,10]
    for Time in time_length:
        enter_dirtime(Time)
        analysis()
        exit_dir()
################################################################            
time()



