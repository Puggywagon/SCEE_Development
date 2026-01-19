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
    if not isinstance(data, dict):
        raise ValueError(f"{path} did not parse to a YAML mapping (dict).")
    return data
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
def exit_dir():
    dire=f'../'
    os.chdir(dire)
################################################################################
def combine_group(g):
    mu_i = g["mu_liquid"].to_numpy()
    se_i = g["mu_se"].to_numpy()

    # Drop rows where se is missing
    mask = np.isfinite(mu_i) & np.isfinite(se_i)
    mu_i = mu_i[mask]
    se_i = se_i[mask]

    N = len(mu_i)
    if N == 0:
        return pd.Series({"mu_liquid_mean": np.nan, "mu_liquid_se": np.nan})

    mu_mean = mu_i.mean()
    mu_se = np.sqrt(np.sum(se_i**2)) / N

    return pd.Series({"mu_liquid_mean": mu_mean, "mu_liquid_se": mu_se})
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
print('With this information are you ready to get started? Yes or No?')
ready=input()
if ready == 'Yes':
    print('Please continue with the next steps')
elif ready == 'No':
    print('Please restart this script once you are ready to start.')
    exit(0)

config= load_yaml("Settings.yml")
user_settings = config["User_Settings"]
State_Conditions = user_settings["State_Conditions"]
T_list=State_Conditions["Temperature"]
P_list = State_Conditions["Pressure"]
Replicas = State_Conditions["Replicas"]
replicas=range(0,f'{Replicas+1}',1)
advanced_settings = config["Advanced_Settings"]

# May need to do something different for when Yes_Gro and Yes for a mixture my steps as that will follow will act as a better workflow better for the other mixture options... May need a different structure a little bit.
if user_settings["Mode"] == "Yes_Gro":
    Yes_Gro=user_settings["Yes_Gro"]
    Gro_File=Yes_Gro["Gro_File"]
    Topology_File=Yes_Gro["Topology_File"]
    solvent=Yes_Gro["solvent"]
    Box_Build='No'
    
elif user_settings["Mode"] == "No_Gro":
    No_Gro=user_settings["No_Gro"]
    Gro_File=No_Gro["Gro_File"]
    Topology_File=No_Gro["Topology_File"]
    solvent=No_Gro["solvent"]
    density=No_Gro["density"]
    mol_mass=No_Gro["mol_mass"]
    solresnametop=No_Gro["solresnametop"]
    Box_Build='Yes'

elif user_settings["Mode"] == "Build_Gro":
    Build_Gro=user_settings["Yes_Gro"]
    Gro_File=Build_Gro["Gro_File"]
    Topology_File=Build_Gro["Topology_File"]
    solvent=Build_Gro["solvent"]
    density=Build_Gro["density"]
    mol_mass=Build_Gro["mol_mass"]
    solresnametop=Build_Gro["solresnametop"]
    solresname=Build_Gro["solresname"]
    solmol=Build_Gro["solmol"]
    Box_Build='Yes'
    Gro_File=Gro_Builder.AA_Structure(solmol,solvent,solresname) #solute #Rename this to handle UA and AA models rather than how I have been doing things..., get this to return the gro file name
    
#####################################################################################
# Where is gromacs
pipe = Popen("/usr/local/gromacs/bin/GMXRC.bash; env", stdout=PIPE, \
shell=True)
output = pipe.communicate()[0]
env = dict((line.decode('utf8').split("=", 1) for line in output.splitlines()))
print(env)
os.environ.update(env)

logfile = open('junk.log', 'w')
cmd = ['which', 'gmx']
check_call(cmd, stdout=logfile, stderr=logfile, env=env)

#####################################################################################
if Box_Build == 'Yes'
    L = advanced_settings["configuration"]["box_length_nm"]
    avo=6.022*10**23
    N=(density*(((L*10**-7)**3))/(mol_mass))*avo
    initial_molecules=int(np.ceil((N+0.05*N))-1)
    print(initial_molecules)
    Pre_Eq_Simulations.Pre_Eq_Solute(Gro_File,Topology_File,L)
    #dipole_model=Pre_Eq_Simulations.get_dipole_model() # I dunno if we need this to do the dielectric constant values
    Pre_Eq_Simulations.insert_Molecules(Gro_File,initial_molecules) 
    Pre_Eq_Simulations.Write_to_top(Topology_File,initial_molecules,solresnametop) 
System_Gro=Pre_Eq_Simulations.Pre_Eq_System(Topology_File)        #Just because they gave us a box I don't trust that they have minimised it well.
######################################################################################################
#Here is the gaussian for pure liquids and mix
if user_settings["Mixture_Loop"]=='No':
    pure_solvent='Yes'
elif user_settings["Mixture_Loop"]=='Yes':
    if user_settings["Mode"] == "No_Gro" or user_settings["Mode"] == "Builds_Gro": 
        pure_solvent='Yes'
    elif user_settings["Mode"] == "Yes_Gro:
        pure_solvent='No'

R = advanced_settings["configuration"]["cluster_radius_nm"]
Configurations = advanced_settings["sampling"]["n_configurations"]
scaling = advanced_settings["electrostatics"]["charge_scaling"]
qr1 = scaling["qr1"]
qr2 = scaling["qr2"]
qr3 = scaling["qr3"]

Di_Const = State_Conditions["Di_Const"]
Ref_Ind = State_Conditions["Ref_Ind"]
cal_diconst=round(Di_Const-(Ref_Ind**2)+1,3)
print(cal_diconst)

#######################################################################
#Some Gaussian needed info here
Molecule.natom=Total_Atoms
Vacuum.natom=Total_Atoms
PCM1.natom=Total_Atoms
PCM2.natom=Total_Atoms

max_jobs=State_Conditions["max_jobs"]
Gaussian_Location =State_Conditions["g09root"]
Scratch_Location =State_Conditions["GAUSS_SCRDIR"]
nproc=8
mem='5GB'


PCM1.nproc = nproc
PCM1.mem = mem
PCM1.max_jobs = max_jobs

PCM1.g09root= Gaussian_Location
PCM1.GAUSS_SCRDIR = Scratch_Location
PCM1.rundir = 'PCM1/'

PCM2.nproc = nproc
PCM2.mem = mem
PCM2.max_jobs = max_jobs

PCM2.g09root= Gaussian_Location
PCM2.GAUSS_SCRDIR = Scratch_Location
PCM2.rundir = 'PCM2/'

Vacuum.nproc = nproc
Vacuum.mem = mem
Vacuum.max_jobs = max_jobs

Vacuum.g09root= Gaussian_Location
Vacuum.GAUSS_SCRDIR = Scratch_Location
Vacuum.rundir = 'Vacuum/'

Molecule.nproc = nproc
Molecule.mem = mem
Molecule.max_jobs = max_jobs

Molecule.g09root= Gaussian_Location
Molecule.GAUSS_SCRDIR = Scratch_Location
