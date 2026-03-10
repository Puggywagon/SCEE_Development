#!/usr/bin/python3
import subprocess
import numpy as np
import os
import shutil
import glob
import matplotlib.pyplot as plt
import pandas as pd

import Pre_Eq_Simulations
import Oniom_Generation_New
import MD_Simulations
import Di_const
import Density
import Dipole_distributions
import H_Bond
import RDF_Script
import Density_Average
import Di_const_Average
import Dipole_distributions_Average
import H_Bond_Average
import plot_rdf
import plot_dipoledist
import SCEE
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
replicas=[1,2,3,4,5,6,7,8,9,10]
T_list=[298] 
P_list = [1.0]

create_dir_Figures()
################################################################################

Pre_Eq_Simulations=Pre_Eq_Simulations.Pre_Eq_Simulations()
Oniom_Generation=Oniom_Generation_New.Oniom_Generation_New()
MD_Simulations=MD_Simulations.MD_Simulations()

Di_const=Di_const.Di_const()
Density=Density.Density()
Dipole_distributions=Dipole_distributions.Dipole_distributions()
H_Bond=H_Bond.H_Bond()
RDF_Script=RDF_Script.RDF_Script()

Density_Average=Density_Average.Density_Average()
Di_const_Average=Di_const_Average.Di_const_Average()
Dipole_distributions_Average=Dipole_distributions_Average.Dipole_distributions_Average()
H_Bond_Average=H_Bond_Average.H_Bond_Average()

plot_rdf=plot_rdf.plot_rdf()
plot_dipoledist=plot_dipoledist.plot_dipoledist()

######################################################################################################

print('Before starting the calculations there is some information and formatting that is required from you. This is important for ensuring the scripts being applied function as intended.')
print("We are able to create the system of interest through our python scripts. If this is a requirement, you will need to provide the gro files containing the AA {central}_AA.gro and UA {solvent}_UAS.gro and associated topology files for these molecules. The central atom will be equilibrated, minimised and have the number of requested solvent atoms added to a 5*5*5 box. Please ensure the central and solvent atoms are given different names in this case. If this isn't required the gro file for the system should be titled either {central}.gro for a system of all the same molecules or {central}_in_{solvent}.gro if for a mixture.")
print("For systems that require the application of both UA and AA molecules that appears differently in the topology file, we request that these are separated into different .itp files labelled in the following format:f'{central}_AA.itp' and f'{solvent}_UA.itp', listed in the topology file using the #include notation. If your system is a mixture it should be labeled using f'{central}_in_{solvent}.top'. Otherwise, if using using a system that appears the same in the topology file when split into UA and AA we request that these are included as f'{central}.top'. You will be asked for the name of your central and solvent molecule names for this purposes.")
print('There are some additional requirements for the formating of your topology file. If a dummy atom is present, ensure the [atomtypes] section, either within the topology file or associated forcefield contains the following columns ;name typ MW q dummy V@nm W@kj/mol where A represents atom and D represents the dummy. If not, the following column structure is accepted ;name typ MW q V@nm W@kj/mol. This is require for the generation of the oniom file.')

print('With this information are you ready to get started? yes or no?')
#ready=input()
ready='yes'

if ready == 'yes':
    print('Please continue with the next steps')
elif ready == 'no':
    print('Please restart this script once you are ready to start.')
    exit(0)

print('How many solvent molecules are surrounding your central atom? 500 - 1000 molecules.')
#Solvent_Molecules=input()
Solvent_Molecules=602

print('Does your system require separated .itp files? yes/no?')
#split=input()
split='yes'

print('What is your central and solvent molecues?')
print('Central atom?')
#central=input()
central='methanol'
print('Surrounding Solvent?')
#solvent=input()
solvent='methanol'

if central == solvent:
    system_title=f'{central}'
else:
    system_title=f'{central}_in_{solvent}'

Configurations=200
Cut_Off_Radius=2.0
Oniom='oniom.inp'
Oniom_Generation.Gen_File(Oniom,Configurations, system_title, Cut_Off_Radius)

f = open(f'{system_title}.top')
Topology = f.read()
f.close()

print('Does the central atom you are using require the application of a dummy atom?')
#dummy=input()
dummy='no'

print('What is the scaling ratio for your charges?')
print('qr1:')
#qr1=input()
print('qr2:')
#qr2=input()
print('qr3:')
#qr3=input()

qr1=1.0
qr2=1.4
qr3=1.8

######################################################################################################
print('Do you need the gro system generated?')
#gen_system=input()
gen_system='yes'

#if gen_system == 'yes':
#    Pre_Eq_Simulations.Pre_Eq_Solute(central)
#    Pre_Eq_Simulations.insert_Molecules(central,solvent,system_title,Solvent_Molecules)
#Pre_Eq_Simulations.Pre_Eq_System(system_yitle,Topology_file)

#dipole_model=Pre_Eq_Simulations.get_dipole_model()
HOMEDIR = os.getcwd()
################################################################################
#create_dir_Simulations()
#for i in replicas:
#    create_dir_reps(i)
#    for T in T_list:
#        create_dir_temps(T)
#        for p in P_list:
#            create_dir_press(p)
            #MD_Simulations.create_mdpfile(HOMEDIR,'junk.mdp',T,p)
            #MD_Simulations.run_md('junk.mdp',HOMEDIR,system_title)
            #MD_Simulations.process_trajectory(system_title)

#            MD_Simulations.process_gro()
#            exit_dir()
#        exit_dir()
#    exit_dir()
#exit_dir()


Model_Dipole=Pre_Eq_Simulations.get_dipole_model_liquid(system_title)

################################################################################ 
if split == 'yes':
    itp_list = Oniom_Generation.get_included_files(Topology)
    for filename in itp_list:
        if filename == f'{central}_AA.itp':
            f1 = open(f'{filename}')
            Solute = f1.read() 
            qmax=Oniom_Generation.Calc_Qmax(Solute,Oniom,dummy,split)   
            Oniom_Generation.QM_Inputs(Solute,Oniom,dummy,split,qr1,qr2,qr3,qmax)   
        elif filename == f'{solvent}_UA.itp':
            f1 = open(f'{filename}')
            Solvent = f1.read()
            f.close()
            Oniom_Generation.MM_Inputs(Solvent,Oniom,split,qr1,qr2,qr3,qmax)  
            Oniom_Generation.Counting_Molecules(Solvent,Oniom,Solvent_Molecules)

elif split == 'No':
    qmax=Oniom_Generation.Calc_Qmax(Solute,Oniom,dummy,split)
    Oniom_Generation.QM_Inputs(Topology,Oniom,dummy,split,qmax)
    Oniom_Generation.MM_Inputs(Topology,Oniom,split,qmax)
    Oniom_Generation.Counting_Molecules(Oniom,Solvent_Molecules)
################################################################################         

ratio=Model_Dipole / qmax
num1=qr1*ratio
num2=qr2*ratio
num3=qr3*ratio

print(f'{ratio},{Model_Dipole},{num1},{num2},{num3}')
################################################################################  
Molecule = SCEE.SCEE()

Molecule.g09root='/home/zoe/Software/Gaussian/g09_pgi/'
Molecule.GAUSS_SCRDIR = '/home/zoe/Research/Gaussian/scratch'

Molecule.rundir = './opt-test/'
dir_list = glob.glob('./Simulations/replica_1/298K/1.0Bar')
print(dir_list)

HOMEDIR = os.getcwd()
print(f'original directory: {HOMEDIR}')



for rundir in dir_list:
    os.chdir('./' + rundir)
    cwd = os.getcwd()
    print(cwd)

    #process *.gro files
    shutil.copy2(HOMEDIR +'/' + 'oniom.inp', '.')
    shutil.copy2(HOMEDIR +'/' + 'shellO', '.')
    shutil.copy2(HOMEDIR +'/' + 'Shell_Oniom.f90', '.')
    Molecule.process_gro(exe=HOMEDIR + '/' + 'shellO',inp=HOMEDIR + '/' + 'oniom.inp') 
    #Molecule.init_v0()
    #Molecule.run_gaussian(step=0)

    #Molecule.init_v1()
    #Molecule.run_gaussian(step=1)
        
    
    df = Molecule.get_multipole_statistics()
    print(df.describe())
    print(df.head())
    muL_list = []
    for index, row in df.iterrows():
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
    
    df['muL'] = muL_list
    df.to_csv('Dipole1.csv', index=False)

    Q1 = df['muL'].quantile(0.25)
    Q3 =df['muL'].quantile(0.75)
    IQR = Q3 - Q1
    lower = Q1 - 1.5*IQR
    upper = Q3 + 1.5*IQR
    print(Q1, Q3, IQR)
    print(lower, upper)
    os.chdir(HOMEDIR)
    cwd = os.getcwd()
    plt.hist(df['dipole_l'], 20)
    plt.hist(df['dipole_m'], 20)
    plt.hist(df['dipole_h'], 20)
    plt.show()
exit(0)
enter_dir_Simulations()
for i in replicas:
    enter_dir_reps(i)
    for T in T_list:
        enter_dir_temps(T)
        for p in P_list:   
            enter_dir_press(p)  
            Di_const.get_epsilon(system_title,i,T,p)
            Density.get_density(system_title,i,T,p)
            Dipole_distributions.get_dipole(i,T,p)
            H_Bond.get_Hbonds(system_title,i,T,p)
            RDF_Script.make_ndx(system_title)
            RDF_Script.rdf_OO(system_title)
            RDF_Script.rdf_OH(system_title)
            RDF_Script.rdf_HH(system_title)
            exit_dir()
        exit_dir()
    exit_dir()
exit_dir()

################################################################################
#Density_Average.Density_Average()
#Di_const_Average.Di_const_Average()
#Dipole_distributions_Average().Dipole_distributions_Average()
#H_Bond_Average.H_bond_Average()
################################################################################
#plot_murho.plot_murho()
#plot_states.plot_states()
#plot_rdf.plot_rdf()
#plot_dipoledist.plot_dipoledist()
################################################################################
