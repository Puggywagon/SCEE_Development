#!/usr/bin/python3
import subprocess
from subprocess import Popen, check_call, PIPE
import numpy as np
import shutil
import glob
import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import os
import socket
import yaml

import Gro_Builder
import Pre_Eq_Simulations
import Oniom_Generation
import MD_Simulations
import SCEE
import Vacuum
import PCM1
import PCM2


################################################################################
def load_yaml(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)
################################################################################
def create_dir_Figures():
    cmd1='mkdir Figures'
    subprocess.run(cmd1, shell=True)
################################################################################
def create_dir_Simulations():
    cmd1='mkdir Simulations'
    subprocess.run(cmd1, shell=True)
    dire='Simulations'
    os.chdir(dire)
################################################################################
def create_dir_reps(i):
    j=i+1
    cmd1=f'mkdir replica_{i}'
    subprocess.run(cmd1, shell=True)
    dire=f'replica_{i}'
    os.chdir(dire)
################################################################################
def create_dir_temps(T):
    cmd1=f'mkdir {T}K'
    subprocess.run(cmd1, shell=True)
    dire=f'{T}K'
    os.chdir(dire)
################################################################################
def create_dir_press(p):
    cmd1=f'mkdir {p}Bar'
    subprocess.run(cmd1, shell=True)
    dire=f'{p}Bar'
    os.chdir(dire)
################################################################################
def enter_dir_Simulations():
    cmd1='mkdir Simulations'
    subprocess.run(cmd1, shell=True)
    dire='Simulations'
    os.chdir(dire)
################################################################################
def enter_dir_reps(i):
    j=i+1
    dire=f'replica_{i}'
    os.chdir(dire)
################################################################################
def enter_dir_temps(T):
    dire=f'{T}K'
    os.chdir(dire)
################################################################################
def enter_dir_press(p):
    dire=f'{p}Bar'
    os.chdir(dire)
################################################################################
def exit_dir():
    dire=f'../'
    os.chdir(dire)
################################################################################
Gro_Builder=Gro_Builder.Gro_Builder()
Pre_Eq_Simulations=Pre_Eq_Simulations.Pre_Eq_Simulations()
Oniom_Generation=Oniom_Generation.Oniom_Generation()
MD_Simulations=MD_Simulations.MD_Simulations()
Vacuum = Vacuum.Vacuum()
Molecule = SCEE.SCEE()
PCM1 = PCM1.PCM1()
PCM2 = PCM2.PCM2()


######################################################################################################
print ()

print('With this information are you ready to get started? Yes or No?')
ready=input()

if ready == 'Yes':
    print('Please continue with the next steps')
elif ready == 'No':
    print('Please restart this script once you are ready to start.')
    exit(0)


print('With this information are you ready to get started? Yes or No?')
System_Build=input()

print('With this information are you ready to get started? Yes or No?')    
Gro_Build=input()

print('Is this system a mixture? Yes or No?')    
MixtureLoop=input()


advanced_settings = load_yaml("Advanced_Settings.yml")
cfg1 = advanced_settings
R = cfg1["configuration"]["cluster_radius_nm"]
Configurations = cfg1["sampling"]["n_configurations"]
scaling = cfg1["electrostatics"]["charge_scaling"]
qr1 = scaling["qr1"]
qr2 = scaling["qr2"]
qr3 = scaling["qr3"]



if System_Build == 'No':   # We do this if they have given us a preformated system in a given box
    yes_gro = load_yaml("Yes_Gro.yml")
    cfg2 = yes_gro 
        
else:
    if Gro_Build == 'No': # We do this if they have given us one or two gro files but haven't done the solvent in a box
        no_gro = load_yaml("No_Gro.yml")
        cfg2 = no_gro 
        
    else:
        build_gro = load_yaml("Build_Gro.yml")    
        cfg2 = build_gro
        
        Gro_Builder.AA_Structure(mol,system_title,resname) #solute
        Gro_Builder.UA_Structure(solmol,system_title,solresname) #solvent 

    L = cfg["configuration"]["box_length_nm"]
    
    avo=6.022*10**23
    N=(density*(((L*10**-7)**3))/(mol_mass))*avo
    initial_molecules=int(np.ceil((N+0.05*N))-1)
    print(initial_molecules)
    

    #Pre-sim for pure solvent
    
    Pre_Eq_Simulations.Pre_Eq_Solute(central,L)
    #dipole_model=Pre_Eq_Simulations.get_dipole_model() # I dunno if we need this to do the dielectric constant values
    Pre_Eq_Simulations.insert_Molecules(central,solvent,system_title,initial_molecules) 
    Pre_Eq_Simulations.Write_to_top(system_title,initial_molecules,solresnametop) 


Pre_Eq_Simulations.Pre_Eq_System(system_title)        #Just because they gave us a box I don't trust that they have minimised it well.
            
cal_diconst=round(exp_diconst-(ref_ind**2)+1,3)
print(cal_diconst)

######################################################################################################
replicas=range(0,f'{Replicas+1}',1)
T_list=Temperature
P_list = Pressure
Cut_Off_Radius=R
################################################################################



# Dunno how to make the pcm1, pcm2 and vacuum steps for this





# Got an idea of how to do this but will be a long step and not needed for this publication
if Mixture_Loop=='Yes':
    Pre_Eq_Simulations.Pre_Eq_Solute(central,L)
    #dipole_model=Pre_Eq_Simulations.get_dipole_model() # I dunno if we need this to do the dielectric constant values
    Pre_Eq_Simulations.insert_Molecules(central,solvent,system_title,initial_molecules) 
    Pre_Eq_Simulations.Write_to_top(system_title,initial_molecules,solresnametop) 
    Pre_Eq_Simulations.Pre_Eq_System(system_title)   





        
