import gromacs
import subprocess
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
################################################################################
def enter_dirP(p):
    dire='{}Bar'.format(p)
    os.chdir(dire)
################################################################################
def enter_dirT(T):
    dire='{}K'.format(T)
    os.chdir(dire)
################################################################################
def enter_dirtime(Time):
    dire='{}ns'.format(Time)
    os.chdir(dire)
################################################################################
def exit_dir():
    dire=f'../'
    os.chdir(dire)
################################################################################
def make_ndx():
    input_1='del 2 \n del 1 \n del 0 \n a O* \n name 0 Oxygen \n a H* \n name 1 Hydrogen \n q'
    tmp = gromacs.tools.Make_ndx(f='Water_QMMM_md3.tpr',
                                o='index.ndx',
                                input=input_1)
    chk, makendx_output, stderr = tmp.run()
################################################################################
def rdf_OO():
    input_str='0\n0\n'
    tmp = gromacs.tools.Rdf(f='Water_QMMM_md3.xtc', 
                                s='Water_QMMM_md3.tpr',
                                n='index.ndx',
                                excl=[],
                                bin=0.004,
                                o='OO_RDF.xvg',
                                input=input_str)
    chk, rdf_output, stderr = tmp.run()
################################################################################
def rdf_OH():
    input_str='0\n 1\n'
    tmp = gromacs.tools.Rdf(f='Water_QMMM_md3.xtc', 
                                s='Water_QMMM_md3.tpr',
                                n='index.ndx',
                                bin=0.004,
                                excl=[],
                                o='OH_RDF.xvg',
                                input=input_str)
    chk, rdf_output, stderr = tmp.run()
################################################################################
def rdf_HH():
    input_str='1\n 1\n'
    tmp = gromacs.tools.Rdf(f='Water_QMMM_md3.xtc', 
                                s='Water_QMMM_md3.tpr',
                                n='index.ndx',
                                bin=0.004,
                                excl=[],
                                o='HH_RDF.xvg',
                                input=input_str)
    chk, rdf_output, stderr = tmp.run()
################################################################################
def plot_rdf_OO(T,p,color):
    data = np.loadtxt("OO_RDF.xvg", skiprows=26) 
    r = data[:, 0]
    gr = data[:, 1]
    plt.plot(r,gr, color=color, label=f'p={p:.0f} bar')
    plt.xlim(0,1.0)
    plt.xlabel('r / nm',fontsize=22)
    plt.ylabel('$g_{OO}(r)$',fontsize=22)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)  
    #plt.title('{}K'.format(T),fontsize=14,weight='bold')
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
################################################################################
def plot_rdf_OH(T,p,color):
    data = np.loadtxt("OH_RDF.xvg", skiprows=26) 
    r = data[:, 0]
    gr = data[:, 1]
    plt.plot(r,gr, color=color, label=f'p={p:.0f} bar')
    
    plt.xlim(0,1.0)
    plt.xlabel('r / nm',fontsize=22)
    plt.ylabel('$g_{OH}(r)$',fontsize=22)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)  
    #plt.title('{}K'.format(T),fontsize=14,weight='bold')
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
################################################################################
def plot_rdf_HH(T,p,color):
    data = np.loadtxt("HH_RDF.xvg", skiprows=25) 
    r = data[:, 0]
    gr = data[:, 1]
    plt.plot(r,gr, color=color, label=f'p={p:.0f} bar')
    
    plt.xlim(0,1.0)
    plt.xlabel('r / nm',fontsize=22)
    plt.ylabel('$g_{HH}(r)$',fontsize=22)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)  
    #plt.title('{}K'.format(T),fontsize=14,weight='bold')
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
################################################################################
def plt_close_OO(T):
    plt.savefig(f'../Analysis/RDFs/OO/{T}_OO_RDF.pdf', dpi=300)
    plt.savefig(f'../Analysis/RDFs/OO/{T}_OO_RDF.png', dpi=300)
    plt.close()
################################################################################
def plt_close_OH(T):
    plt.savefig(f'../Analysis/RDFs/OH/{T}_OH_RDF.pdf', dpi=300)
    plt.savefig(f'../Analysis/RDFs/OH/{T}_OH_RDF.png', dpi=300)
    plt.close()
################################################################################
def plt_close_HH(T):
    plt.savefig(f'../Analysis/RDFs/HH/{T}_HH_RDF.pdf', dpi=300)
    plt.savefig(f'../Analysis/RDFs/HH/{T}_HH_RDF.png', dpi=300)
    plt.close()
##############################################################################################
def pressure_rdf_OO(T):
    p_ref = [1,2,5]
    P_list = []
    for scale in np.logspace(0, 3, 4):
       P_list += [p*scale for p in p_ref]
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(P_list))]
    plt.figure(figsize=(15, 12))
    for p, color in zip(P_list, color_list):
        enter_dirP(p)
        #make_ndx()
        #rdf_OO()
        plot_rdf_OO(T,p,color)
        exit_dir()
##############################################################################################
def pressure_rdf_OH(T):
    p_ref = [1,2,5]
    P_list = []
    for scale in np.logspace(0, 3, 4):
       P_list += [p*scale for p in p_ref]
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(P_list))]
    plt.figure(figsize=(15, 12))
    for p, color in zip(P_list, color_list):
        enter_dirP(p)
        #rdf_OH()
        plot_rdf_OH(T,p,color)
        exit_dir()
##############################################################################################
def pressure_rdf_HH(T):
    p_ref = [1,2,5]
    P_list = []
    for scale in np.logspace(0, 3, 4):
       P_list += [p*scale for p in p_ref]
    cmap = plt.cm.rainbow
    color_list = [cmap(x) for x in np.linspace(0.0, 1.0, len(P_list))]
    plt.figure(figsize=(15, 12))
    for p, color in zip(P_list, color_list):
        enter_dirP(p)
        #rdf_HH()
        plot_rdf_HH(T,p,color)
        exit_dir()
##############################################################################################
def temperature():
    T_list=[298] #Temp in K
    for temp in range(400,1001,100):
        T_list.append(temp)
    for T in T_list:
        enter_dirT(T)        
        pressure_rdf_OO(T)
        plt_close_OO(T)
        pressure_rdf_OH(T)
        plt_close_OH(T)
        pressure_rdf_HH(T)
        plt_close_HH(T)
        exit_dir()
################################################################
def time():
    time_length = [5]
    for Time in time_length:
        enter_dirtime(Time)
        temperature()
        exit_dir()
################################################################            
time()
                                
