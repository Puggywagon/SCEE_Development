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
    
    for line in dipole_lines[::-1]:  # Loop through lines in reverse order
        if "G_k" in line:
            words = line.split()
            index = words.index("G_k") + 2  # Get index of element after "G_k"
            big_G = float(words[index])
            break  # Exit loop after finding the value
    
    for line in dipole_lines[::-1]:  # Loop through lines in reverse order
        if "g_k" in line:
            words = line.split()
            index = words.index("g_k") + 2  # Get index of element after "G_k"
            Little_g = float(words[index])
            break  # Exit loop after finding the value

    return epsilon,big_G,Little_g
###############################################################################################
def build_df(T,p,epsilon, big_G, Little_g):
    data = {
        'Temperature': T,
        'Pressure': p,
        'Epsilon': epsilon,
        'big G': big_G,
        'Little g': Little_g
    },
    df=pd.DataFrame(data)
    return df
###############################################################################################
def df_csv(df):
    with open(f"../../Analysis/G_Factors.csv", "a", newline="") as file:  # Open in append mode
        df.to_csv(file, header=False, index=False)
        return df
###############################################################################################
def run_results(T,p):
    epsilon, big_G, Little_g = get_epsilon()
    df=build_df(T,p,epsilon, big_G, Little_g)
    df_csv(df)
##############################################################################################
def pressure(T):
    p_ref = [1,2,5]
    P_list = []
    for scale in np.logspace(0, 3, 4):
       P_list += [p*scale for p in p_ref]
    for p in P_list:
        enter_dirP(p)
        run_results(T,p)
        exit_dir()
##############################################################################################
def temperature():
    T_list=[298] #Temp in K
    for temp in range(400,1001,100):
        T_list.append(temp)
    for temp in range(273,374,10):
        T_list.append(temp)
    for T in T_list:
        enter_dirT(T)
        pressure(T)
        print('complete')
        exit_dir()
################################################################
def time():
    time_length = [5,10]
    for Time in time_length:
        enter_dirtime(Time)
        create_analysis()
        temperature()
        exit_dir()
################################################################            
time()
