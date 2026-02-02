#!/usr/bin/python3
import gromacs
import subprocess
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
################################################################################
def create_repdir(i):
    cmd1='mkdir replica_{}'.format(i)
    subprocess.run(cmd1, shell=True)
################################################################################
def enter_repdir(i):
    dire='replica_{}'.format(i)
    os.chdir(dire)
################################################################################
def create_dirP(p):
    cmd1='mkdir {}Bar'.format(p)
    subprocess.run(cmd1, shell=True)
################################################################################
def enter_dirP(p):
    dire='{}Bar'.format(p)
    os.chdir(dire)
################################################################################
def create_dirT(T):
    cmd1='mkdir {}K'.format(T)
    subprocess.run(cmd1, shell=True)
################################################################################
def enter_dirT(T):
    dire='{}K'.format(T)
    os.chdir(dire)
################################################################################
def create_dirTime(Time):
    cmd1='mkdir {}ns'.format(Time)
    subprocess.run(cmd1, shell=True)
################################################################################
def enter_dirtime(Time):
    dire='{}ns'.format(Time)
    os.chdir(dire)
################################################################################
def exit_dir():
    dire=f'../'
    os.chdir(dire)
################################################################################
def run_md(mdpfile):
    print('performing molecular dynamics')

    mdpfile = f'junk.mdp'
    grofile = f'../../../../initial.gro'
    topol = f'../../../../system.top'

    runname = 'Water_QMMM_md3'
    edrfile = f'{runname}.edr'
    groout = f'{runname}.gro'
    xtcfile = f'{runname}.xtc'
    trrfile = f'{runname}.trr'
    tprfile = f'{runname}.tpr'
    # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
    gromacs.grompp(f=mdpfile, c=grofile, p=topol, o=tprfile, maxwarn=-2)

    # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
    gromacs.mdrun(v=True, s=tprfile, c=groout, o=trrfile, x=xtcfile, e=edrfile)
################################################################################
def create_mdpfile(mdpfile, T, p, steps):
    T_kelvin = T  # temperature / K
    p_bar = p  # pressure / bar
    replace_dict = {'TEMPERATURE': f'{T_kelvin}',
                    'PRESSURE': f'{p_bar}',
                    'STEPS': f'{steps}'}
    search_text = 'TEMPERATURE', 'PRESSURE', 'STEPS'
    replace_text = f'{T_kelvin}', f'{p_bar}', f'{steps}'
    with open('../../../../template.mdp', 'r') as file:
        data = file.read()
        for search_text, replace_text in replace_dict.items():
            data = data.replace(search_text, replace_text)
    with open(mdpfile, 'w') as file:
        file.write(data)


################################################################################
def process_trajectory():
    print('processing trajectory')
    input_str = '0\n'
    tmp = gromacs.tools.Trjconv(f='Water_QMMM_md3.xtc',
                                s='Water_QMMM_md3.tpr',
                                o=f'conf_.gro',
                                input=input_str,
                                b=1000, dt=20, sep=True)
    chk, stdout, stderr = tmp.run()


################################################################################
def process_gro():
    print('processing *.gro files')
    logfile = open('junk.log', 'w')
    cmd5 = ["gfortran -o Shell_Oniom Shell_Oniom.f90 && ./Shell_Oniom"]
    subprocess.call(cmd5, shell=True, stdout=logfile, stderr=logfile)
################################################################################
def plot_Pressure(p,T,Time):
    plt.plot(T,p,marker='o', markerfacecolor='none')
    
    plt.xlabel(r'Temperature / K')
    plt.xlim([200, 1101])    
    
    plt.ylabel(r'Pressure / bar')
    plt.ylim([0.5, 1.0e4])
    plt.yscale('log')
    
    plt.savefig(f'../../../../Sim_States_{Time}.pdf')
    plt.savefig(f'../../../../Sim_States_{Time}.png')
################################################################################
def plot_experiment():
    df = pd.read_csv('../experimental_results.csv')
    df.columns = ['Temperature (K)', 'Pressure (bar)', 'Phase']
    Tc = 647.096      # Critical temperature (Tc) 647.096 K
    pc = 220.640      # Critical pressure (Pc)    220.640 bar
    rhoc = 17.873728  # Critical density (Dc)     17.873728 mol/l

    tmp = df[df['Phase']=='liquid']
    x_list = list(tmp['Temperature (K)']) + [Tc]
    y_list = list(tmp['Pressure (bar)']) + [pc]
    plt.plot(x_list, y_list, ls='solid', color='black')
    plt.plot([Tc], [pc], marker='o', color='black')

    plt.xlim([200, 800])
    plt.ylim([0.01, 1.0e4])
    plt.yscale('log')
    plt.ylabel('pressure / bar')
    plt.xlabel('temperature / K')
################################################################################
def plot_vega():
    df = pd.read_csv('../vega_results.csv')
    df.columns = ['T / K', 'p / bar', 'rhol','rhog']
      
    Tc = 640
    pc = 146
    rhoc = 0.31   # g cm^{-3} 
    
    x_list =  [Tc] + list(df['T / K'])
    y_list =  [pc] + list(df['p / bar'])
    plt.plot(x_list, y_list, color='black', ls='dashed')
    plt.plot([Tc], [pc], marker='o', color='black')
################################################################################
def plot_Update(p,T,Time):    
    plt.plot(T,p,marker='o', color='red')
    
    plt.xlabel(r'Temperature / K')
    plt.xlim([200, 1101])    
    
    plt.ylabel(r'Pressure / bar')
    plt.yscale('log')
    plt.ylim([0.5, 1.0e4])
    
    plt.savefig(f'../../../Sim_States_{Time}.pdf')
    plt.savefig(f'../../../Sim_States_{Time}.png')
################################################################################
def plot_close():
    plt.show()
    plt.close()
################################################################################
def run_project(steps,Time,t,p):
    create_mdpfile('junk.mdp',steps,Time,t,p)
    run_md('junk.mdp')
    process_trajectory()
    process_gro()
##############################################################################################
def pressure(steps,Time,T):
    p_ref = [1,2,5]
    P_list = []
    for scale in np.logspace(0, 3, 4):
       P_list += [p*scale for p in p_ref]
    for p in P_list:
        create_dirP(p)
        enter_dirP(p)
        plot_Pressure(p,T,Time)
        #run_project(steps,Time,T,p)
        plot_Update(p,T,Time)
        exit_dir()
##############################################################################################
def temperature(steps,Time):
    T_list=[298] #Temp in K
    for temp in range(400,1001,100):
        T_list.append(temp)
    for temp in range(273,374,10):
        T_list.append(temp)
    for T in T_list:
        create_dirT(T)
        enter_dirT(T)
        pressure(steps,Time,T)
        exit_dir()
################################################################
def time(steps, time_step):
    times = (time_step * steps) / 1000
    index_levels = ['Time', 'T', 'p']
    time_length =[5]

    for Time in time_length:
        if Time == times:
             create_dirTime(Time)
             enter_dirtime(Time)
             plot_experiment()
             plot_vega()
             temperature(steps, Time)
             plot_close()
             exit_dir()
        else:
            return 0
####################################################################################
def steps():
    no_steps = [2500000]
    time_step = 0.002
    for steps in no_steps:
        time(steps, time_step)
####################################################################################
Replica_list=[1,2,3]
for i in Replica_list:
    create_repdir(i)
    enter_repdir(i)
    steps()
    exit_dir()
