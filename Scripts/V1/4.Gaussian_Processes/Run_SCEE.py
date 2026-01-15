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

import SCEE
import Vacuum
import PCM1
import PCM2


Vacuum = Vacuum.Vacuum()
Molecule = SCEE.SCEE()
PCM1 = PCM1.PCM1()
PCM2 = PCM2.PCM2()


Model_Dipole=Pre_Eq_Simulations.get_dipole_model_liquid(system_title)
ratio=Model_Dipole / qmax
num1=qr1*ratio*qmax
num2=qr2*ratio*qmax
num3=qr3*ratio*qmax
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
