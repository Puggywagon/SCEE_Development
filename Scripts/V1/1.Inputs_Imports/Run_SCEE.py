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

#This was for telling it what version of fromacs is present.
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
Gro_Builder=Gro_Builder.Gro_Builder()
Pre_Eq_Simulations=Pre_Eq_Simulations.Pre_Eq_Simulations()
Oniom_Generation=Oniom_Generation.Oniom_Generation()
MD_Simulations=MD_Simulations.MD_Simulations()
Vacuum = Vacuum.Vacuum()
Molecule = SCEE.SCEE()
PCM1 = PCM1.PCM1()
PCM2 = PCM2.PCM2()

######################################################################################################

print('Before starting the calculations there is some information and formatting that is required from you. This is important for ensuring the scripts being applied function as intended.')
print("We are able to create the system of interest through our python scripts. If this is a requirement, you will need to provide the gro files containing the AA {central}_AA.gro and UA {solvent}_UA.gro and associated topology files for these molecules. The central atom will be equilibrated, minimised and have the number of requested solvent atoms added to a 5*5*5 box. Please ensure the central and solvent atoms are given different names in this case. If this isn't required the gro file for the system should be titled either {central}.gro for a system of all the same molecules or {central}_in_{solvent}.gro if for a mixture.")
print("For systems that require the application of both UA and AA molecules that appears differently in the topology file, we request that these are separated into different .itp files labelled in the following format:f'{central}_AA.itp' and f'{solvent}_UA.itp', listed in the topology file using the #include notation. If your system is a mixture it should be labeled using f'{central}_in_{solvent}.top'. Otherwise, if using using a system that appears the same in the topology file when split into UA and AA we request that these are included as f'{central}.top'. You will be asked for the name of your central and solvent molecule names for this purposes.")
print('There are some additional requirements for the formating of your topology file. If a dummy atom is present, ensure the [atomtypes] section, either within the topology file or associated forcefield contains the following columns ;name typ MW q dummy V@nm W@kj/mol where A represents atom and D represents the dummy. If not, the following column structure is accepted ;name typ MW q V@nm W@kj/mol. This is require for the generation of the oniom file.') # The approach was designed to be able to use a different central molecule from the surrounding solvent for if we wanted to do mixtures.

print('With this information are you ready to get started? yes or no?')
#ready=input()
ready='yes'
# My attempt to give a user an out before going through the different steps if the script.
if ready == 'yes':
    print('Please continue with the next steps')
elif ready == 'no':
    print('Please restart this script once you are ready to start.')
    exit(0)

print('Does your system require separated .itp files? yes/no?')
#split=input()
split='yes'

print('What is your central and solvent molecues?')
print('Central molecules?')
#central=input()
central='Aniline'
# A smile string is a letter representation of a molecular structure say CCO would be ethanol while, c1ccccc1 is benzene. This is used in the building of the input files. Do need to rework this a little bit so we could do if you need it fill in some of these steps, if not don't. Best place to obtain these is pubmed, AI doesn't generate these particularly well yet
print('Its SMILE String?')
#mol = Chem.MolFromSmiles(input())
mol = 'c1ccccc1N'

print('Its residue name?')
#resname=input()
resname='ANI'
solresnametop='Ani'
print('Surrounding Solvent?')
#solvent=input()
solvent='Aniline'

print('Its SMILE String?')
#solmol = Chem.MolFromSmiles(input())
solmol = 'c1ccccc1N'

print('Its residue name?') 
#solresname=input()
solresname='ANI'
# This step tells the script if it is a mixture or not i.e methanol in water, the calculation would need to run a bit differently for this, run SCEE for the water, insert the methanol and then run a slightly different approach for methanol in the water
if central == solvent:
    system_title=f'{central}'
else:
    system_title=f'{central}_in_{solvent}'

print('What is your solvent keyword for Gaussian?') #Gaussian requires a solvent keyword to run the calculation but doesn't always have one for the molecule you are running. I have gotten away with a close approximate i.e. benzylalcohol for phenol but this isn't always the case. If you can find a structure for how these work I can help make a script to generate these...
#sol_keyword=input()
sol_keyword='Aniline'

print('Does the central atom you are using require the application of a dummy atom?')
#dummy=input()
dummy='no' # Some molecules like water require the addition of a dummy atom to balance the charges for the MD, but these are not needed in the QM. This little bit tells the oniom generation bit of the script if the dummy is there or not and changes how the oniom input is generated based on this.

print('What is the scaling ratio for your charges?') # These are used for the scaling of charges, I have not changed the values when doing my simulations and shouldn't impact your calculation results too much
print('qr1:')
#qr1=input()
qr1=1.0

print('qr2:')
#qr2=input()
qr2=1.35

print('qr3:')
#qr3=input()
qr3=1.7

print('Do you need the gro system generated?')
#gen_system=input()
gen_system='yes' # This activates the script that generates gro files for you, you can change this to no if you are using software to generate these for you

print('What is the experimental density in g cm$^-^1$?')
#density=input()
density=1.022 # This is used in the calculation of the number of molecules required in your box.

print('What is the molar mass of your?')
#mol_mass=input()
mol_mass=93.13 # The molar mass is something you can either google or workout yourself (also known as gfm)

print('What is the experimental dielectric constant of your molecule?')
#exp_diconst=input()
exp_diconst=6.890 # This is something used in calculating the simulation dielectric constant which is something my scripts don't do - ask Miguel/Leo

print('What is the experimental refractive index of your molecule?')
#ref_ind=input()
ref_ind=1.586 # This is something used in calculating the simulation dielectric constant which is something my scripts don't do - ask Miguel/Leo

print('What is the experimental dipole moment of your molecule in the gas phase?')
#gas_dipole=input()
gas_dipole=1.49 #This should be the experimental value but these can be difficult to find in the literature depending on your molecule

print('How many replicas would you like to run?')
#reps=input()
replicas=range(3,6,1) #How many times you want to repeat your simulation, this gives 5 replicas
print('What temperatures (in K) would you like to run your simulations to run at? (please use 1, 2, ... format)')
#temps=input()
temps='298.0' #If you wanted to explore a range of temperatures, change temps to '298.0, new_value' is in Kelvin
T_list=[t.strip()) for t in temps.split(',')]

print('What pressures (in K) would you like to run your simulations to run at? (please use 1, 2, ... format)')
#press=input()
press='1.0' #Similarly if you want to explore a range of temperature change press to '1.0, new_value' is in Bar
P_list = [p.strip()) for p in press.split(',')]

################################################################################
#This step will open the topology file and counts how many heavy atoms are present in the topology file
f = open(f'{system_title}.top')
Topology = f.read()
f.close()
if split == 'yes':
    itp_list = Oniom_Generation.get_included_files(Topology)
    for filename in itp_list:
        if filename == f'{central}_AA.itp':
            f1 = open(f'{filename}')
            Solute = f1.read()
            heavy_atoms=Oniom_Generation.Calc_Heavys(Solute,dummy,split)
elif split == 'No':
    heavy_atoms=Oniom_Generation.Calc_Heavys(Solute,dummy,split)
# We use the smile strings and collected information in the generation of the gro files of the solute and solvent molecules
Gro_Builder.AA_Structure(mol,system_title,resname)
Gro_Builder.UA_Structure(solmol,system_title,solresname)
# This is a calculation step of the dielectric constant, please ask Miguel about this bit as I am not working on dielectric constants
cal_diconst=round(exp_diconst-(ref_ind**2)+1,3)
print(cal_diconst)
# Here we calculate the number solvent molecules required for a simulation box of L*L*L dimensions according the experimental density.
#Previously the calculation was
#R=1.1*heavy_atoms*0.2    #This was used for calculating the cutoff radius but we found it was better using fixed values
#L=2R+0.1  #Calculating the required length of the box based on the cutoff radius

L=6.2
avo=6.022*10**23
N=(density*(((L*10**-7)**3))/(mol_mass))*avo
initial_molecules=int(np.ceil((N+0.05*N))-1)
print(initial_molecules)

Configurations=200
Cut_Off_Radius=2.8 #=R

#Note that I have an equivalent for this bit for Vacuum, PCM1 and PCM2. Figured this is something that should be kept for the inputs section though.
Molecule.g09root='/home/zoe/Software/Gaussian/g09_pgi/'
Molecule.GAUSS_SCRDIR = '/home/zoe/Research/Gaussian/scratch'
Molecule.rundir = 'opt-test-SCEE/'

hostname = socket.gethostname()
print(f"Running on host: {hostname}")
if 'Tower3' in hostname:
    # Tower3: 48 cores, 125 GiB
    Molecule.nproc = 8
    Molecule.mem = '5GB'
    Molecule.max_jobs = 4

    # (Optionally match Vacuum/PCM1/PCM2 to something similar later)
else:
    # Home / office desktop: 24 cores, 15 GiB
    Molecule.nproc = 8
    Molecule.mem = '3GB'
    Molecule.max_jobs = 3

