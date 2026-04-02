#!/usr/bin/python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

class plot_states(object):
    def __init__(self):
        pass
################################################################################
################################################################################
################################################################################
    def plot_experiment(self):
        df = pd.read_csv(nist_vp_csvfile)
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
        plt.yscale('log')
        plt.ylabel('pressure / bar')
        plt.xlabel('temperature / K')

    
################################################################################
    def plot_vega(self):
        df1 = pd.read_csv(tip4p_vp_csvfile)
#        df1.columns = ['T / K', 'p / bar', 'rhol','rhog']
      
        Tc = 640
        pc = 146
        rhoc = 0.31   # g cm^{-3} 
    
        x_list =  [Tc] + list(df1['T / K'])
        y_list =  [pc] + list(df1['p / bar'])
        plt.plot(x_list, y_list, color='black', ls='dashed')
        plt.plot([Tc], [pc], ls='--',marker='o', color='black',label='TIP4P/2005')
        
        plt.xlim([250, 1101])
        plt.ylim([0.5, 1.0e4])
        plt.yscale('log')

    
################################################################################
    def liquid_boundary(self,T,p):
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
    def gas_boundary(self,T,p):
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
################################################################################
################################################################################
################################################################################
    def plot_states(self):

        parser = argparse.ArgumentParser()

        parser.add_argument('--show', type=str, default='True', help='plot figures')

        args = parser.parse_args()

        show = (args.show == 'True')

        fontsize = 10
        legendsize = 6

        nist_vp_csvfile = './Results/NIST_vaporpressure.csv'
        tip4p_vp_csvfile = './Results/TIP4P_vaporpressure.csv'

    
####################################################################################
        fontsize = 10
        legendsize = 6
    
        plt.figure(figsize=(3.5, 3))

        plot_experiment()
        plot_vega()

################################################################################
        plt.tight_layout()

        plt.savefig(f'./Figures/states.pdf', dpi=300)
        plt.savefig(f'./Figures/states.png', dpi=300)
        if (show):
            plt.show()
        else:    
            plt.clf()


################################################################################
################################################################################
################################################################################


