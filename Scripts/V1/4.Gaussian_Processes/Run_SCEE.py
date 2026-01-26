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

import Gaussian_Calculations


Gauss=Gaussian_Calculations.Gaussian_Calculations()


Model_Dipole=Pre_Eq_Simulations.get_dipole_model_liquid(system_title)
################################################################################  
# This step calculates the dipole moment of the vacuum using the output of the presim steps using a single molecule
Gauss.g09root='/home/zoe/Software/Gaussian/g09_pgi'
Gauss.GAUSS_SCRDIR = '/home/zoe/Research/Gaussian/scratch'

hostname = socket.gethostname()
print(f"Running on host: {hostname}")

mu_Vacuum=Gauss.init(Gaus='Vacuum',sol_keyword)
#PCM1
PCM1= Gauss.init(Gaus='PCM1',sol_keyword,exp_diconst)
################################################################################ 
PCM2= Gauss.init(Gaus='PCM2',sol_keyword,cal_diconst)
################################################################################ 
# These steps calculates the dipole moment using the SCEE approach for the ~200 configurations in each replica folder. Sometimes you will run into a problem in one replica and so the structure of the dir_list steps allows you to comment out the replicas you have completed and restart from here

dir_list = glob.glob('./Simulations/replica_1/*K/*.0Bar')
dir_list += glob.glob('./Simulations/replica_2/*K/*.0Bar')
dir_list += glob.glob('./Simulations/replica_3/*K/*.0Bar')
dir_list += glob.glob('./Simulations/replica_4/*K/*.0Bar')
dir_list += glob.glob('./Simulations/replica_5/*K/*.0Bar')
print(dir_list)


for rundir in dir_list:
    # Here we copy the files into directory we are running these steps
    os.chdir('./' + rundir)
    cwd = os.getcwd()
    print(cwd)
    
    ratio=Model_Dipole / qmax
    num1=qr1*ratio*qmax
    num2=qr2*ratio*qmax
    num3=qr3*ratio*qmax
    Gauss.process_gro(exe=HOMEDIR +'/' + 'shellO',inp=HOMEDIR +'/' +'oniom.inp') 
    Gauss.init(Gaus='SCEE_V0')
    SCEE=Gauss.init(Gaus='SCEE_V1')
    
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
