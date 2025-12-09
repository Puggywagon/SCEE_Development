#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
################################################################################
class plot_rdf(object):
    def __init__(self):
        pass
################################################################################
################################################################################
################################################################################
    def get_data(self,tag):
        file_list = glob.glob("./Results/Replica_*/*K/*Bar/{tag}_RDF.xvg")
        
        data = np.loadtxt(f"../Results/Replica_*/*K/*Bar/{tag}_RDF.xvg",
                          skiprows=26) 
        r = data[:, 0]
        gr = data[:, 1]
        return r, gr


################################################################################
################################################################################
################################################################################
    def plot_rdf(self):
        parser = argparse.ArgumentParser()
    
        parser.add_argument('--show', type=str, default='True', help='plot figures')
    
        args = parser.parse_args()
    
        show = (args.show == 'True')

        fontsize = 10
        legendsize = 6
        ticksize = 10

################################################################################
################################################################################
################################################################################
        fig, ax = plt.subplots(3, 1, figsize=[3.5, 10])

        i = 0

        r, gr = get_data('OO')
        ax[0].plot(r,gr, color='red')
        ax[i].text(0.9, 0.9, '(g)', transform=ax[i].transAxes, fontsize=fontsize)
        ax[i].tick_params(axis='x', which='both', direction='in',
                  labelsize=ticksize, labelbottom=False)
        ax[i].set_xlim([0, 1])
        ax[i].set_ylabel('$g_{OO}(r)$',fontsize=fontsize)
        ax[i].legend(loc='upper center', prop={'size': legendsize})


        i = 1
        r, gr = get_data('OH')
        ax[i].plot(r,gr, color='blue')
        ax[i].text(0.9, 0.9, '(h)', transform=ax[i].transAxes, fontsize=fontsize)
        ax[i].tick_params(axis='x', which='both', direction='in',
                      labelsize=ticksize, labelbottom=False)
        ax[i].set_ylabel('$g_{OH}(r)$',fontsize=fontsize)
        ax[i].set_xlim([0, 1])

        i = 2
        r, gr = get_data('HH')
        ax[i].plot(r,gr, color='green')
        ax[i].text(0.9, 0.9, '(i)', transform=ax[i].transAxes, fontsize=fontsize)
        ax[i].tick_params(axis='x', which='both', direction='in',
                          labelsize=ticksize, labelbottom=True)
        ax[i].set_xlabel('r / nm',fontsize=fontsize)
        ax[i].set_ylabel('$g_{HH}(r)$',fontsize=fontsize)
        ax[i].set_xlim([0, 1])


################################################################################
        plt.tight_layout()
        fig.subplots_adjust(hspace=0)

        plt.savefig(f'./Figures/rdf.pdf', dpi=300)
        plt.savefig(f'./Figures/rdf.png', dpi=300)
        if (show):
            plt.show()
        else:    
            plt.clf()


################################################################################
################################################################################
################################################################################
