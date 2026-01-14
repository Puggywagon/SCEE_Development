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

import Gro_Builder
import Pre_Eq_Simulations
import Oniom_Generation
import MD_Simulations
import SCEE
import Vacuum
import PCM1
import PCM2

pipe = Popen("/usr/local/gromacs/bin/GMXRC.bash; env", stdout=PIPE, \
shell=True)
output = pipe.communicate()[0]
env = dict((line.decode('utf8').split("=", 1) for line in output.splitlines()))
print(env)
os.environ.update(env)

logfile = open('junk.log', 'w')
cmd = ['which', 'gmx']
check_call(cmd, stdout=logfile, stderr=logfile, env=env)
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
Pre_Eq_Simulations=Pre_Eq_Simulations.Pre_Eq_Simulations()
MD_Simulations=MD_Simulations.MD_Simulations()


# We then minimise and equalibriate the generated solute molecule through a sequence of simulation steps designed using Cecilia's MD Tutorial (I can send you these). We start with a single molecule and insert this molecule into a box containing the calculated number of molecules from the previous step. Then we minimise this box and obtain the dipole moment of the model for the vacuum phase.
Pre_Eq_Simulations.Pre_Eq_Solute(central,L)
dipole_model=Pre_Eq_Simulations.get_dipole_model() #I have discovered that this step stops later steps from being able to use -v
print(dipole_model)
Pre_Eq_Simulations.insert_Molecules(central,solvent,system_title,initial_molecules)
Pre_Eq_Simulations.Write_to_top(system_title,initial_molecules,solresnametop) # If having to rerun the presim I comment this section out as if it keeps running the topology file becomes in correct.
Pre_Eq_Simulations.Pre_Eq_System(system_title)


HOMEDIR = os.getcwd()
################################################################################
# These steps performs the production run for the MD simulations. We use template.mdp to loop through temperature (T) and pressure (p) for i number of replicas. We then extract the ~200 configs from the trajectory and process them in preparation for the gaussian calculation. Some molecules require a lower timestep or a different barrostat to run stably and so I have two versions of the template mdp.
create_dir_Simulations()
for i in replicas:
    create_dir_reps(i)
    for T in T_list:
        create_dir_temps(T)
        for p in P_list:
            create_dir_press(p)
            MD_Simulations.create_mdpfile(HOMEDIR,'junk.mdp',T,p)
            MD_Simulations.run_md('junk.mdp',HOMEDIR,system_title)
            MD_Simulations.process_trajectory(system_title)
            MD_Simulations.process_gro()
            exit_dir()
        exit_dir()
    exit_dir()
exit_dir()
################################################################################  

