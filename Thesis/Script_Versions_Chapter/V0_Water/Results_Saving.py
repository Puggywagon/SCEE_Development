import gromacs
import subprocess
import numpy as np
import os
import pandas as pd
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
def create_analysis():
    cmd1 = 'mkdir Analysis'
    subprocess.run(cmd1, shell=True)
################################################################################
def exit_dir():
    dire=f'../'
    os.chdir(dire)
################################################################################
# Dipole Moment and Dielectric Constant
def get_epsilon():
    gromacs.environment.flags['capture_output'] = True
    input_str = '0\n'
    tmp = gromacs.tools.Dipoles(f='Water_QMMM_md3',
                                s='Water_QMMM_md3',
                                input=input_str)
    chk, dipole_output, stderr = tmp.run()

    f = open(f'Dipole.txt', 'w')
    f.write(dipole_output)
    f.close()

    dipole_lines = dipole_output.splitlines()  # Split text into lines
    epsilon_line = dipole_lines[-1]  # Assuming epsilon is on the last line
    epsilon = float(epsilon_line.split()[2])  # Extract the numerical value

    return epsilon
################################################################################
# Density
def get_density():
    gromacs.environment.flags['capture_output'] = True
    input_str = 'Density 0\n'
    tmp = gromacs.tools.Energy(f='Water_QMMM_md3',
                               s='Water_QMMM_md3',
                               input=input_str)
    chk, density_output, stderr = tmp.run()

    f = open(f'Density.txt', 'w')
    f.write(density_output)
    f.close()

    density_lines = density_output.splitlines()
    density_line = density_lines[9]  # Assuming density is on the 5th line
    density = float(density_line.split()[1])

    return density
########################################################################################################################
# Hbonds
def get_Hbonds():
    gromacs.environment.flags['capture_output'] = True
    input_str = '0 0\n'
    tmp = gromacs.tools.Hbond(f='Water_QMMM_md3',
                              s='Water_QMMM_md3',
                              input=input_str)
    chk, hbond_output, stderr = tmp.run()
    f = open(f'HBond.txt', 'w')
    f.write(hbond_output)
    f.close()

    lines = hbond_output.splitlines()
    average_hbonds_line = lines[-1]  # Assuming average is on the second-to-last line
    hbond = float(average_hbonds_line.split()[-5])  # Extract the numerical value
    total_hbond = 1800
    percent_hbond = (hbond / total_hbond) * 100
    Nearest_neighbour=(hbond / total_hbond)*4

    return hbond, percent_hbond, Nearest_neighbour
###############################################################################################
def get_dipole():
    if os.path.isfile('Dipole.csv'):
        df = pd.read_csv('Dipole.csv')
        df.columns = ['config', 'dipole_l', 'dipole_m', 'dipole_h','muL']
        dipole_moments=df['muL']
        dipole=np.mean(dipole_moments)
    else:
        dipole=1.855
    return dipole
###############################################################################################
def build_df(Time,T,p,density,epsilon, hbond,percent_hbond, Nearest_neighbour,dipole):
    data = {
        'Time': Time,
        'Temperature': T,
        'Pressure': p,
        'Density': density,
        'Epsilon': epsilon,
        'Number of H-Bonds': hbond,
        'Percent of H-Bonds': percent_hbond,
        'Nearest Neighbour': Nearest_neighbour,
        'Dipole Moment': dipole
    },
    df=pd.DataFrame(data)
    return df
###############################################################################################
def df_csv(df):
    with open(f"../../Analysis/data.csv", "a", newline="") as file:  # Open in append mode
        df.to_csv(file, header=False, index=False)
        return df
###############################################################################################
def run_results(Time,T,p):
    density = get_density()
    epsilon = get_epsilon()
    hbond,percent_hbond, Nearest_neighbour = get_Hbonds()
    dipole = get_dipole()
    df=build_df(Time,T,p,density,epsilon,hbond,percent_hbond,Nearest_neighbour,dipole)
    df_csv(df)
##############################################################################################
def pressure(Time,T):
    p_ref = [1,2,5]
    P_list = []
    for scale in np.logspace(0, 3, 4):
       P_list += [p*scale for p in p_ref]
    for p in P_list:
        enter_dirP(p)
        run_results(Time,T,p)
        exit_dir()
##############################################################################################
def temperature(Time):
    T_list=[298] #Temp in K
    for temp in range(400,1001,100):
        T_list.append(temp)
    for temp in range(273,374,10):
        T_list.append(temp)
    for T in T_list:
        enter_dirT(T)
        pressure(Time,T)
        exit_dir()
################################################################
def time():
    time_length = [5,10]
    for Time in time_length:
        enter_dirtime(Time)
        create_analysis()
        temperature(Time)
        exit_dir()
################################################################            
time()

