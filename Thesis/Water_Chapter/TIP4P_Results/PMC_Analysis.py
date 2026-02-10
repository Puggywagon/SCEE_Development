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
def read_data():
    df2 = pd.read_csv('../../pcm_dipoles.csv')
    df2.columns = ['Temperature', 'Pressure', 'Density', 'Volume','Internal Energy','Enthalpy', 'Entropy','Cv','Cp','Sound','Joule','Viscosity','Therm','Phase','epsilon','g','PMC Dipole','My Dipole','Nearest','diff']
    return df2
################################################################################
def enter_dirtime(Time):
    dire='{}ns'.format(Time)
    os.chdir(dire)
#################################################################################
def plot_pmcdipole(df2,T_list):
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
    plt.figure(figsize=(15, 12))
    for T, color in zip(T_list, color_list):
        tmp = df2[df2['Temperature'] == T]
        plt.plot(tmp['PMC Dipole'], tmp['My Dipole'], ls='',
                 color=color, marker='o',label=r'$T='+f'{T}$ K')
        #plt.plot(tmp['Nearest'], equal, ls='-',color='black')
    extrapolated_points=[(min(df2['PMC Dipole']),1.85),(max(df2['PMC Dipole']),2.2)]
    plt.plot([point[0] for point in extrapolated_points],[point[1] for point in extrapolated_points], ls='-', color='black', linewidth=1) 
    plt.ylabel(r'SCEE Dipole Moment / D',fontsize=22)
    plt.xlabel(r'PMC Dipole Moment / D',fontsize=22)
    plt.xlim([min(df2['PMC Dipole']),max(df2['PMC Dipole'])])
    plt.ylim([1.6, 3.0])
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)  
    plt.xlim([1.84, 2.17])
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    Min=min(T_list)
    Max=max(T_list)
    plt.savefig(f'./PMC/PMC_Dipoles_Newline_{Min}_{Max}_not.pdf', dpi=300)
    plt.savefig(f'./PMC/PMC_Dipoles_Newline_{Min}_{Max}_not.png', dpi=300)
    plt.close()
#################################################################################
def plot_densitypmc(df2,T_list):
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
    plt.figure(figsize=(15, 12))
    for T, color in zip(T_list, color_list):
        tmp = df2[df2['Temperature'] == T]
        plt.plot(tmp['Density'], tmp['PMC Dipole'], ls='None',
                 color=color, marker='o',label=r'$T='+f'{T}$ K')

    plt.ylabel(r'PMC Dipole Moment / D',fontsize=22)
    #plt.ylim([1.6, 3.0])
    plt.xlabel(r'Density / kg m$^{-3}$',fontsize=22)
    plt.xlim([0, 1200])
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)  
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    Min=min(T_list)
    Max=max(T_list)
    plt.savefig(f'./PMC/density_pmc_{Min}_{Max}_not.pdf', dpi=300)
    plt.savefig(f'./PMC/density_pmc_{Min}_{Max}_not.png', dpi=300)
    plt.close()
 #################################################################################
def plot_nearestneighbourPMC(df2,T_list):
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(T_list))]
    plt.figure(figsize=(15, 12))
    for T, color in zip(T_list, color_list):
        tmp = df2[df2['Temperature'] == T]
        plt.scatter(tmp['Nearest'], tmp['PMC Dipole'], ls='None',
                 color=color, facecolor='none', edgecolors=color, marker='o',label=r'$T='+f'{T}$ K')
        plt.scatter(tmp['Nearest'], tmp['My Dipole'], ls='None',
                 color=color, marker='o')
    a,b=np.polyfit(df2['Nearest'],df2['My Dipole'],1)
    x=df2['Nearest']
    plt.plot(x,(a*x+b), ls='-',color='black',label='SCEE')

    plt.ylabel(r'Dipole Moment / D',fontsize=22)
    plt.ylim([1.6, 3.0])
    plt.xlabel(r'Nearest H-Bonded Neighbour',fontsize=22)
    plt.xlim([0, 4.0])
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)  
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    Min=min(T_list)
    Max=max(T_list)
    plt.savefig(f'./PMC/Nearest_Neighbour_PMC_{Min}_{Max}_not.pdf', dpi=300)
    plt.savefig(f'./PMC/Nearest_Neighbour_PMC_{Min}_{Max}_not.png', dpi=300)
    plt.close()       
############################################################################
def analysis():
    T_list=[298]
    for temp in range(400,1001,100):
        T_list.append(temp)
    enter_analysis()
    df2=read_data()
    plot_pmcdipole(df2,T_list)
    plot_densitypmc(df2,T_list)
    plot_nearestneighbourPMC(df2,T_list)
        
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



