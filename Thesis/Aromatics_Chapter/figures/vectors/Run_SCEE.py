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
import sys
import socket
from shutil import which
import yaml
from pathlib import Path

import Simulations_Analysis
import Gaussian_Calculations

#######################################################################
def expand_topology_with_itps(topology_file: str) -> str:
    topology_path = Path(topology_file).resolve()
    seen = set()

    def walk(file_path: Path) -> str:
        file_path = file_path.resolve()

        if file_path in seen:
            return ""
        seen.add(file_path)

        text = file_path.read_text()
        includes = Analysis.get_included_files(text)

        blocks = []
        if file_path.name != "forcefield.itp":
            blocks.append(text)

        for inc in includes:
            inc_path = (file_path.parent / inc).resolve()
            if file_path.name == "forcefield.itp":
                if inc_path.name == "ffnonbonded.itp":
                    blocks.append(inc_path.read_text())
                continue

            blocks.append(walk(inc_path))

        return "\n\n".join([b for b in blocks if b.strip()])
    return walk(topology_path)
#######################################################################


Analysis=Simulations_Analysis.Simulations_Analysis() #Make this into an analysis script rather than what I have it now.
Gauss=Gaussian_Calculations.Gaussian_Calculations()

#######################################################################
max_jobs= 3
Gaussian_Location ="/home/zoe/Software/Gaussian/g09_pgi"
Scratch_Location = "/home/zoe/Research/Gaussian/scratch"
Gromacs_Location ="/home/zoe/Software/gromacs-2024.2/bin/GMXRC.bash"
nproc=8
mem='5GB'

Gauss.nproc = nproc
Gauss.mem = mem
Gauss.max_jobs = max_jobs

Gauss.g09root= Gaussian_Location
Gauss.GAUSS_SCRDIR = Scratch_Location

pipe = Popen(f"{Gromacs_Location}; env", stdout=PIPE, \
shell=True)
output = pipe.communicate()[0]
env = dict((line.decode('utf8').split("=", 1) for line in output.splitlines()))
os.environ.update(env)

logfile = open('junk.log', 'w')
cmd = ['which', 'gmx']
check_call(cmd, stdout=logfile, stderr=logfile, env=env) 

#######################################################################
qr1 = 1.0
qr2 = 1.35
qr3 = 1.7


collection={'Molecule':[],'Temperatures':[],'Pressures':[],'Replicas':[],'Number of Configurations':[],'Exp Di Constant':[],'Refractive Index':[],'vacuum dipole moment model':[],'liquid dipole moment model':[],'Experimental gas dipole':[],'cal_diconst':[],'mu_Vacuum':[],'PCM1':[],'PCM2':[],'muL_SCEE':[],'delta_mu':[],'mu_liquid':[],'Model Epsilon':[],'Dielectric Constant Correction':[],'STDEV':[],'SE':[]}
collection=pd.DataFrame(collection)

Molecules=["Methylamine","Ethylamine","Propylamine","Butylamine","1-Butanol","1-Propanol","2-Butanol","2-Propanol","Ethanol","Acetone","Butanone","Phenol","Acetophenone","Aniline"]

T_list=[298.0,298.0,298.0,298.0,298.0,298.0,298.0,298.0,298.0,298.0,298.0,330.15,298.0,298.0]
P_list = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
Di_Consts = [4.9912,6.3,4.9912,4.62,17.332,20.524,15.944,19.264,24.852,20.493,18.246,9.78,17.44,6.890]
Ref_Inds = [1.37,1.3663, 1.388,1.401,1.399,1.384,1.600,1.377,1.361,1.3589,1.379,1.542,1.534,1.586]
gas_dipoles=[1.31,1.22,1.17,1.00,1.660,1.680,1.7,1.580,1.690,2.88,2.78,1.38,3.06,1.49]

MolDIR = os.getcwd()

print(len(T_list),len(P_list),len(Di_Consts),len(Ref_Inds),len(gas_dipoles))

for Molecule,T,P,Di_Const,Ref_Ind,gas_dipole in zip(Molecules,T_list,P_list,Di_Consts,Ref_Inds,gas_dipoles):
    print(Molecule)
    os.chdir(f'../../{Molecule}')
    collection['Molecule']=Molecule
    collection['Temperatures']=T
    collection['Pressures']=P
    collection['Exp Di Constant']=Di_Const
    collection['Refractive Index']=Ref_Ind
    collection['Experimental gas dipole']=gas_dipole

    HOMEDIR = os.getcwd()
    Topology_File=f'{Molecule}.top'
    Topologys = expand_topology_with_itps(Topology_File)
    with open("Concat_Top.top", "w") as f:
        f.write(Topologys)

    Oniom='oniom.inp'
    f = open("Concat_Top.top")
    Topology = f.read()
    Total_Atoms,qmax=Analysis.Calc_Qmax(Topology) 
    f.close()

    cal_diconst=round(Di_Const-(Ref_Ind**2)+1,3)
    dipole_model=Analysis.get_dipole_model() 

    Gauss.natom=Total_Atoms
    mu_Vacuum=Gauss.init1(Gaus='Vacuum')
    PCM1= Gauss.init2(exp_diconst,Gaus='PCM1')
    PCM2= Gauss.init3(cal_diconst,Gaus='PCM2')


    collection['cal_diconst']=cal_diconst
    collection['vacuum dipole moment model']=dipole_model
    collection['mu_Vacuum']=mu_Vacuum
    collection['PCM1']=PCM1
    collection['PCM2']=PCM2
    dir_list = glob.glob('./Simulations/replica_*/*K/*.0Bar')    



    Model_Dipole, Epsilon=Analysis.get_dipole_model_liquid(system_title)
    collection['liquid dipole moment model']=Model_Dipole

    collection['Model Epsilon']=Epsilon
    plot={}
    plot=pd.DataFrame(plot)
    count=0
    for rundir in dir_list:
        os.chdir('./' + rundir)
        cwd = os.getcwd()
        print(cwd)
        Gauss.process_gro(exe=HOMEDIR +'/' + 'shellO',inp=HOMEDIR +'/' +'oniom.inp') 
        SCEE=Gauss.init4(Molecule,Gaus='SCEE_V1')

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
    
        plot['muL_SCEE']=df1['muL_SCEE']
        plot['delta_mu'] = delta_mu_list
        plot['mu_liquid'] = mu_liquid_list
        count+=1
        df.to_csv('Dipole_Calculations.csv', index=False)  
        os.chdir(HOMEDIR)  
    os.chdir(MolDIR)  
    

    collection['Replicas']=count
    collection['Number of Configurations']=len(plot['mu_liquid'])     
  
    collection['muL_SCEE']=plot['muL_SCEE'].mean()    
    collection['delta_mu']=plot['delta_mu'].mean()    
    collection['mu_liquid']=plot['mu_liquid'].mean()  

    collection['Dielectric Constant Correction'] = Ref_Ind**2 + (collection['mu_liquid'] / Model_Dipole) * (epsilon + 1)
    collection['STDEV']=np.std(collection['mu_liquid'])
    collection['SE']=collection['STDEV']/count  
    Analysis.plot_dipoledist(Molecule)
    
collection.to_csv("../Results.csv", index=False)    
