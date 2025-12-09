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
print("We are able to create the system of interest through our python scripts. If this is a requirement, you will need to provide the gro files containing the AA {central}_AA.gro and UA {solvent}_UAS.gro and associated topology files for these molecules. The central atom will be equilibrated, minimised and have the number of requested solvent atoms added to a 5*5*5 box. Please ensure the central and solvent atoms are given different names in this case. If this isn't required the gro file for the system should be titled either {central}.gro for a system of all the same molecules or {central}_in_{solvent}.gro if for a mixture.")
print("For systems that require the application of both UA and AA molecules that appears differently in the topology file, we request that these are separated into different .itp files labelled in the following format:f'{central}_AA.itp' and f'{solvent}_UA.itp', listed in the topology file using the #include notation. If your system is a mixture it should be labeled using f'{central}_in_{solvent}.top'. Otherwise, if using using a system that appears the same in the topology file when split into UA and AA we request that these are included as f'{central}.top'. You will be asked for the name of your central and solvent molecule names for this purposes.")
print('There are some additional requirements for the formating of your topology file. If a dummy atom is present, ensure the [atomtypes] section, either within the topology file or associated forcefield contains the following columns ;name typ MW q dummy V@nm W@kj/mol where A represents atom and D represents the dummy. If not, the following column structure is accepted ;name typ MW q V@nm W@kj/mol. This is require for the generation of the oniom file.')

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
T_list=[temps]

print('What pressures (in K) would you like to run your simulations to run at? (please use 1, 2, ... format)')
#press=input()
press='1.0' #Similarly if you want to explore a range of temperature change press to '1.0, new_value' is in Bar
P_list = [press]

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
L=6.2
avo=6.022*10**23
N=(density*(((L*10**-7)**3))/(mol_mass))*avo
initial_molecules=int(np.ceil((N+0.05*N))-1)
print(initial_molecules)
# We then minimise and equalibriate the generated solute molecule through a sequence of simulation steps designed using Cecilia's MD Tutorial (I can send you these). We start with a single molecule and insert this molecule into a box containing the calculated number of molecules from the previous step. Then we minimise this box and obtain the dipole moment of the model for the vacuum phase.
#Pre_Eq_Simulations.Pre_Eq_Solute(central,L)
dipole_model=Pre_Eq_Simulations.get_dipole_model()
print(dipole_model)
#Pre_Eq_Simulations.insert_Molecules(central,solvent,system_title,initial_molecules)
#Pre_Eq_Simulations.Write_to_top(system_title,initial_molecules,solresnametop)
#Pre_Eq_Simulations.Pre_Eq_System(system_title)


HOMEDIR = os.getcwd()
################################################################################
# These steps performs the production run for the MD simulations. We use template.mdp to loop through temperature (T) and pressure (p) for i number of replicas. We then extract the ~200 configs from the trajectory and process them in preparation for the gaussian calculation. Some molecules require a lower timestep or a different barrostat to run stably and so I have two versions of the template mdp.
#create_dir_Simulations()
#for i in replicas:
#    create_dir_reps(i)
#    for T in T_list:
#        create_dir_temps(T)
#        for p in P_list:
#            create_dir_press(p)
#            MD_Simulations.create_mdpfile(HOMEDIR,'junk.mdp',T,p)
#            MD_Simulations.run_md('junk.mdp',HOMEDIR,system_title)
#            MD_Simulations.process_trajectory(system_title)
#            MD_Simulations.process_gro()
#            exit_dir()
#        exit_dir()
#    exit_dir()
#exit_dir()
################################################################################  
# Here we generate the input file for our inhouse Oniom Script using information that can be extracted through the topology files. Here we select the number of configurations and cluster radius (r). We also extract the total number of atoms within the box
Configurations=200
Cut_Off_Radius=2.8
Oniom='oniom.inp'
solvent_molecules=initial_molecules
Oniom_Generation.Gen_File(Oniom,Configurations, system_title, Cut_Off_Radius)
if split == 'yes':
    itp_list = Oniom_Generation.get_included_files(Topology)
    for filename in itp_list:
        if filename == f'{central}_AA.itp':
            f1 = open(f'{filename}')
            Solute = f1.read() 
            qmax,Total_Atoms=Oniom_Generation.Calc_Qmax(Solute,Oniom,dummy,split)   
            Oniom_Generation.QM_Inputs(Solute,Oniom,dummy,split,qr1,qr2,qr3,qmax)     
        elif filename == f'{solvent}_UA.itp':
            f1 = open(f'{filename}')
            Solvent = f1.read()
            f.close()
            Oniom_Generation.MM_Inputs(Solvent,Oniom,split,qr1,qr2,qr3,qmax)  
            Oniom_Generation.Counting_Molecules(Solvent,Oniom,initial_molecules)
elif split == 'No':
    qmax,Total_Atoms=Oniom_Generation.Calc_Qmax(Solute,Oniom,dummy,split)
    Oniom_Generation.QM_Inputs(Topology,Oniom,dummy,split,qmax)
    Oniom_Generation.MM_Inputs(Topology,Oniom,split,qmax)
    Oniom_Generation.Counting_Molecules(Oniom,initial_molecules)
################################################################################  
# Here we calculate the model dipole moment of the production run box and use this in calculating the scaled charges
Model_Dipole=Pre_Eq_Simulations.get_dipole_model_liquid(system_title)
ratio=Model_Dipole / qmax
num1=qr1*ratio*qmax
num2=qr2*ratio*qmax
num3=qr3*ratio*qmax
Molecule.natom=Total_Atoms
Vacuum.natom=Total_Atoms
PCM1.natom=Total_Atoms
PCM2.natom=Total_Atoms
################################################################################  
# This step calculates the dipole moment of the vacuum using the output of the presim steps using a single molecule
Vacuum.g09root='/home/zoe/Software/Gaussian/g09_pgi'
Vacuum.GAUSS_SCRDIR = '/home/zoe/Research/Gaussian/scratch'

Vacuum.rundir = 'opt-test_vacuum/'
hostname = socket.gethostname()
print(f"Running on host: {hostname}")

if 'Tower3' in hostname:
    # Tower3: 48 cores, 125 GiB
    Vacuum.nproc = 8
    Vacuum.mem = '5GB'
    Vacuum.max_jobs = 4
else:
    # Home / office desktop: 24 cores, 15 GiB
    Vacuum.nproc = 8
    Vacuum.mem = '3GB'
    Vacuum.max_jobs = 3
#process *.gro files
Vacuum.init_v0(system_title,sol_keyword)
Vacuum.run_gaussian(step=0)
    
df4 = Vacuum.get_multipole_statistics()
print(df4.describe())
print(df4.head())
muG_list=[]
for index, row in df4.iterrows():
    x = np.array([num1, num2, num3])
    y = np.array([row['dipole_l'], row['dipole_m'], row['dipole_h']])
    coeff = np.polyfit(x, y, 2)
    b = coeff[1] - 1
    if (coeff[0] >= 1.0e-6):
         muG = (-b - np.sqrt(b*b-4*coeff[0]*coeff[2])) / (2*coeff[0])
    else:
        coeff = np.polyfit(x, y, 1)
        b = coeff[0] - 1
        muG = - coeff[1]/b
    muG_list.append(muG)
df4['muG'] = muG_list
df4.to_csv('Dipole_Vacuum.csv', index=False)
mu_Vacuum=df4['muG']
################################################################################  
# PCM1 and PCM2 calculate the dipole moment using a dielectric continuum I will explain why further down
PCM1.g09root='/home/zoe/Software/Gaussian/g09_pgi/'
PCM1.GAUSS_SCRDIR = '/home/zoe/Research/Gaussian/scratch'
PCM1.rundir = 'opt-test-PCM1/'
hostname = socket.gethostname()
print(f"Running on host: {hostname}")
if 'Tower3' in hostname:
    # Tower3: 48 cores, 125 GiB
    PCM1.nproc = 8
    PCM1.mem = '5GB'
    PCM1.max_jobs = 4
else:
    # Home / office desktop: 24 cores, 15 GiB
    PCM1.nproc = 8
    PCM1.mem = '3GB'
    PCM1.max_jobs = 3
#PCM1
PCM1.init_v0(system_title,sol_keyword,exp_diconst)
PCM1.run_gaussian(step=0)
          
df2 = PCM1.get_multipole_statistics()
print(df2.describe())
print(df2.head())
muL_list1 = []
for index, row in df2.iterrows():
    x = np.array([num1, num2, num3])
    y = np.array([row['dipole_l'], row['dipole_m'], row['dipole_h']])
    coeff = np.polyfit(x, y, 2)
    b = coeff[1] - 1
    if (coeff[0] >= 1.0e-6):
        muL = (-b - np.sqrt(b*b-4*coeff[0]*coeff[2])) / (2*coeff[0])
    else:
        coeff = np.polyfit(x, y, 1)
        b = coeff[0] - 1
        muL = - coeff[1]/b
    muL_list1.append(muL)
    
df2['muL_PCM1'] = muL_list1
df2.to_csv('Dipole_PCM1.csv', index=False)
PCM1= df2['muL_PCM1']
################################################################################ 
PCM2.g09root='/home/zoe/Software/Gaussian/g09_pgi/'
PCM2.GAUSS_SCRDIR = '/home/zoe/Research/Gaussian/scratch'
PCM2.rundir = 'opt-test-PCM2/'
hostname = socket.gethostname()
print(f"Running on host: {hostname}")

if 'Tower3' in hostname:
    # Tower3: 48 cores, 125 GiB
    PCM2.nproc = 8
    PCM2.mem = '5GB'
    PCM2.max_jobs = 4

else:
    # Home / office desktop: 24 cores, 15 GiB
    PCM2.nproc = 8
    PCM2.mem = '3GB'
    PCM2.max_jobs = 3
#PCM2 
PCM2.init_v0(system_title,sol_keyword,cal_diconst)
PCM2.run_gaussian(step=0)

df3 = PCM2.get_multipole_statistics()
print(df3.describe())
print(df3.head())
muL_list2 = []
for index, row in df3.iterrows():
    x = np.array([num1, num2, num3])
    y = np.array([row['dipole_l'], row['dipole_m'], row['dipole_h']])
    coeff = np.polyfit(x, y, 2)
    b = coeff[1] - 1
    if (coeff[0] >= 1.0e-6):
         muL = (-b - np.sqrt(b*b-4*coeff[0]*coeff[2])) / (2*coeff[0])
    else:
        coeff = np.polyfit(x, y, 1)
        b = coeff[0] - 1
        muL = - coeff[1]/b
    muL_list2.append(muL)
   
df3['muL_PCM2'] = muL_list2
df3.to_csv('Dipole_PCM2.csv', index=False)
PCM2= df3['muL_PCM2']
################################################################################ 
# These steps calculates the dipole moment using the SCEE approach for the ~200 configurations in each replica folder. Sometimes you will run into a problem in one replica and so the structure of the dir_list steps allows you to comment out the replicas you have completed and restart from here

dir_list = glob.glob('./Simulations/replica_1/*K/*.0Bar')
dir_list += glob.glob('./Simulations/replica_2/*K/*.0Bar')
dir_list += glob.glob('./Simulations/replica_3/*K/*.0Bar')
dir_list += glob.glob('./Simulations/replica_4/*K/*.0Bar')
dir_list += glob.glob('./Simulations/replica_5/*K/*.0Bar')
print(dir_list)

HOMEDIR = os.getcwd()
print(f'original directory: {HOMEDIR}')

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


for rundir in dir_list:
    # Here we copy the files into directory we are running these steps
    os.chdir('./' + rundir)
    cwd = os.getcwd()
    print(cwd)
    #process *.gro files
    shutil.copy2(HOMEDIR +'/' + 'oniom.inp', '.')
    shutil.copy2(HOMEDIR +'/' + 'shellO', '.')
    shutil.copy2(HOMEDIR +'/' + 'Shell_Oniom.f90', '.')
    
    #SCEE
    #Here we use our oniom input file and oniom inhouse script to convert the configurations extracted from the format from gromacs to a format accepted by gaussian. It then runs the dipole moment calculation on each config using each of the scaled charges
    Molecule.process_gro(exe=HOMEDIR +'/' + 'shellO',inp=HOMEDIR +'/' +'oniom.inp') 
    Molecule.init_v0()
    Molecule.run_gaussian(step=0)
    Molecule.init_v1()
    Molecule.run_gaussian(step=1)
    df1 = Molecule.get_multipole_statistics()
    print(df1.describe())
    print(df1.head())
    # Here we extract the dipole moments calculated and perform a quadratic fit to obtain the SCEE dipole moment which is then saved to the Dipole_SCEE.csv. Occasionally the value is outwith this fit and results in a blank box in the csv. If your molecule does not have a dipole like methane etc, change this to a linear fit otherwise you will start getting some strange numbers when the dipole moment should be ~0 D.
    muL_list = []
    for index, row in df1.iterrows():
        x = np.array([num1, num2, num3])
        y = np.array([row['dipole_l'], row['dipole_m'], row['dipole_h']])
        coeff = np.polyfit(x, y, 2)
        b = coeff[1] - 1
        if (coeff[0] >= 1.0e-6):
             muL = (-b - np.sqrt(b*b-4*coeff[0]*coeff[2])) / (2*coeff[0])
        else:
            coeff = np.polyfit(x, y, 1)
            b = coeff[0] - 1
            muL = - coeff[1]/b
        muL_list.append(muL)
    df1['muL_SCEE'] = muL_list
    df1.to_csv('Dipole_SCEE.csv', index=False)
    
    # Depending on the basis set used in the SCEE approach, it is more accurate to calculate the induced dipole moment (delta_mu) and use this to calculate the liquid dipole moment (mu_liquid) using the experimental gas phase values, rather than reporting the SCEE dipole moment directly (mu_SCEE). And so I have this step here to calculate this. If you have questions about this bit I recommend speaking with Miguel about it.
    delta_mu_list=[]
    mu_liquid_list=[]
    for i in df1['muL_SCEE']:
        mu_SCEE=i
        delta_mu=(mu_SCEE*(PCM1[0]/PCM2[0]))-mu_Vacuum[0]
        mu_liquid=delta_mu+gas_dipole
        delta_mu_list.append(delta_mu)
        mu_liquid_list.append(mu_liquid)
    
    data={}
    df=pd.DataFrame(data)
    
    df['config']=df1['config']
    df['dipole_l']=df1['dipole_l']
    df['dipole_m']=df1['dipole_m']
    df['dipole_h']=df1['dipole_h']
    df['muL_SCEE']=df1['muL_SCEE']
    df['delta_mu'] = delta_mu_list
    df['mu_liquid'] = mu_liquid_list
    
    df.to_csv('Dipole_Calculations.csv', index=False)  
    os.chdir(HOMEDIR)  
################################################################################ 
