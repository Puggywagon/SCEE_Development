#!/usr/bin/python3

import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pylab as plt

class plot_hydrogenbond(object):
    def __init__(self):
        pass
################################################################################
################################################################################
################################################################################
    def plot_vle(self):
#        header = sim[0]
#        sim.pop(0)
#        df_sim = pd.DataFrame(sim, columns=header)
    
        Tc = 647.096      # Critical temperature (Tc) 647.096 K
        pc = 220.640      # Critical pressure (Pc)    220.640 bar
        rhoc = 17.873728  # Critical density (Dc)     17.873728 mol/l
    
        tmp = df_nist[df_nist['Phase']=='liquid']
        x_list = list(tmp['Temperature (K)']) + [Tc]
        y_list = list(tmp['Pressure (bar)']) + [pc]
        ax.plot(x_list, y_list, ls='solid', color='black')
        ax.plot([Tc], [pc], marker='o', color='black')


    # Vega et al. 2006
        Tc = 640
        pc = 146
        rhoc = 0.31   # g cm^{-3}  
        x_list =  [Tc] + list(df_sim['T / K'])
        y_list =  [pc] + list(df_sim['p / bar'])    
#        ax.plot(x_list, y_list, marker='x', ms=5, color='black', ls='dashed')
        ax.plot(x_list, y_list, color='black', ls='dashed')
        ax.plot([Tc], [pc], marker='o', color='black')


################################################################################
################################################################################
################################################################################
    def plot_hydrogenbond(self):
        parser = argparse.ArgumentParser()

        parser.add_argument('--show', type=str, default='True', help='plot figures')

        args = parser.parse_args()
    
        show = (args.show == 'True')


        fontsize = 10
        ticksize = 10
        legendsize = 6
        
        csvfile = './Results/TIP4P_Results/pcm_dipoles.csv'
        nist_csvfile = './Results/NIST_vaporpressure.csv'
        tip4p_csvfile = './Results/TIP4P_vaporpressure.csv'
    

################################################################################
################################################################################
        df = pd.read_csv(csvfile, header=None)
        df.columns = ['Temperature', 'Pressure', 'Density', 'Volume','Internal Energy','Enthalpy', 'Entropy','Cv','Cp','Sound','Joule','Viscosity','Therm','Phase','epsilon','g','PMC Dipole','My Dipole','Nearest']


        df_nist = pd.read_csv(nist_csvfile)
        df_sim = pd.read_csv(tip4p_csvfile)
    

################################################################################
        cmap = plt.cm.rainbow
        fig, ax = plt.subplots(figsize=[3.5, 3])

        for index, row in df.iterrows():
            n_HB = row['Nearest']
            color = cmap(n_HB/4)
            x_data = [row['Temperature']]
            y_data = [row['Pressure']]
            ax.plot(x_data, y_data, ls='None', marker='o', color=color)


        norm = mpl.colors.Normalize(vmin=0, vmax=4)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        colorbar = plt.colorbar(sm, ax=ax)
        colorbar.set_label(r'$\langle N_{\rm neig}\rangle$')


        plot_vle()


        # interpolation
        df['log pressure'] = np.log10( df['Pressure'] )
        T_list = df['Temperature'].drop_duplicates().to_list()
        T_list.sort()
        x_list = []
        y_list = []
        for T in T_list:
            tmp = df[df['Temperature']==T]
            if(tmp['Nearest'].min() < 1.0):
                log_p = np.interp(x=1.0, xp=tmp['Nearest'], fp=tmp['log pressure'])
                p = 10**(log_p)
                x_list.append(T)
                y_list.append(p)
        T_min = min(x_list)
        print(T_min)
        tmp = df_sim[df_sim['T / K'] <= T_min]
        tmp = tmp.sort_values('T / K')
        tmp_list = tmp['T / K'].to_list()
        print(tmp_list)
        x_list = tmp_list + x_list
        tmp_list = tmp['p / bar'].to_list()
        y_list = tmp_list + y_list
        plt.fill_between(x_list, y_list, 0, color='purple', alpha=0.5)
        print(x_list)
        print(y_list)

        plt.yscale('log')
        plt.ylim([0.5, 1.0e4])
        plt.xlabel('temperature / K')
        plt.ylabel('pressure / bar')
        plt.tight_layout()

        plt.savefig('./Figures/hydrogenbond.pdf')
        plt.savefig('./Figures/hydrogenbond.png')
        if (show):
            plt.show()
        else:    
            plt.clf()


################################################################################
################################################################################
################################################################################
