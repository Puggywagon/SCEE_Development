import gromacs
import subprocess
import numpy as np
import os
import pandas as pd
import math
import matplotlib.pyplot as plt
import glob
################################################################################
def temp_dipdist(T):
    csv_list = []
    for p in [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0]:
        csv_list += glob.glob(f'./5ns/{T}K/{p}Bar/Dipole.csv')   
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(csv_list))]
    plt.figure(figsize=(14, 6))
    for csvfile, color in zip(csv_list, color_list):
        data = csvfile.split('/')
        T = float( data[-3].replace('K', '') )
        p = float( data[-2].replace('Bar', '') )
        tmp = pd.read_csv(csvfile)
        tmp.columns = ['config', 'dipole_l', 'dipole_m', 'dipole_h','muL']
        counts,bins = np.histogram(tmp['muL'], bins=10)
        x_list = [0.5 * (b1 + b2) for b1, b2 in zip(bins[:-1], bins[1:])]
        plt.plot(x_list, counts, color=color)
    
    plt.xlabel('dipole moment / D',fontsize=12)
    plt.ylabel('probability density / D$^{-1}$',fontsize=12)
    plt.title('{}K'.format(T),fontsize=14,weight='bold')
    plt.xlim(1.65,3.75)
    #plt.x(0.05)
    plt.margins(x=0,y=0)
    plt.autoscale(enable=True, axis='y')
    plt.savefig(f'./5ns/Analysis/Dipole_Distributions/Pressure_Distributions/Dip_dist_T_{T:.0f}K.pdf')
    plt.savefig(f'./5ns/Analysis/Dipole_Distributions/Pressure_Distributions/Dip_dist_T_{T:.0f}K.png')
    plt.close()    
################################################################              
for T in [298]:
    temp_dipdist(T)
       
       







