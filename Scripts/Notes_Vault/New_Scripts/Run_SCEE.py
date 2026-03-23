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


import Gro_Builder
import Gro_Simulations
import Simulations_Analysis
import Oniom_Generation
import Gaussian_Calculations

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

Gro_Builder=Gro_Builder.Gro_Builder()
MD=Gro_Simulations.Gro_Simulations()
Analysis=Simulations_Analysis.Simulations_Analysis() 
Oniom_Generation=Oniom_Generation.Oniom_Generation()
Gauss=Gaussian_Calculations.Gaussian_Calculations()

######################################################################################################
print('To start the process, you need to complete the Yaml script based on your system requirements. If you have already completed this please type Yes, If you need to go back and finish doing this please respond No.')
#ready=input()
ready='Yes'
if ready == 'Yes':
    print('This script will continue with the process.')
elif ready == 'No':
    print('Please restart this script once you are ready to start.')
    exit(0)

config= load_yaml("Settings.yml")
user_settings = config["User_Settings"]
State_Conditions = user_settings["State_Conditions"]
T_list=State_Conditions["Temperature"]
P_list = State_Conditions["Pressure"]
Replicas = State_Conditions["Replicas"]+1
replicas=range(0,Replicas,1)
advanced_settings = config["Advanced_Settings"]
print('here:', user_settings['Mode'])
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
    print('Building *.gro file')
    Build_Gro=user_settings["Build_Gro"]
    Topology_File=Build_Gro["Solvent_Topology_File"]
    solvent=Build_Gro["solvent"]
    density=Build_Gro["density"]
    mol_mass=Build_Gro["mol_mass"]
    solresnametop=Build_Gro["solresnametop"]
    solresname=Build_Gro["solresname"]
    solmol=Build_Gro["solmol"]
    Box_Build='Yes'
    Gro_File=Gro_Builder.AA_Structure(solmol,solvent,solresname) #solute #Rename this to handle UA and AA models rather than how I have been doing things..., get this to return the gro file name #Need to figure out how to maintain consistency between gro structure and top file.
    print('RDKit generates an image of the molecule you have generated with your smile string. Please take a minute to check the structure matches your desired molecule and topology file. If you are ready to continue type \'Yes\', if you need to make adjustments to your smile string please type \'No\'.')
    #ready=input()

    ready='Yes'
    if ready == 'Yes':
        print('The script will continue on with the next steps of the system.')
    elif ready == 'No':
        print('Please restart this script after you have made your adjustments.')
        exit(0)

        
print('setting GROMACS run environment')
# LL: I don't have the GMXRC.bash script in my GROMACS installation
Gromacs_Location =State_Conditions["Grom_Location"]
print(f'Gromacs_Location={Gromacs_Location}')

pipe = Popen(f"{Gromacs_Location}; env", stdout=PIPE, \
shell=True)
output = pipe.communicate()[0]
env = dict((line.decode('utf8').split("=", 1) for line in output.splitlines()))
os.environ.update(env)

logfile = open('junk.log', 'w')
cmd = ['which', 'gmx']
check_call(cmd, stdout=logfile, stderr=logfile, env=env)    

#####################################################################################
if Box_Build == 'Yes':
    L = advanced_settings["configuration"]["box_length_nm"]
    avo=6.022e23
    N=(density*(((L*10**-7)**3))/(mol_mass))*avo
    initial_molecules=int(np.ceil((N+0.05*N))-1)
    print(initial_molecules)
    MD.run_md1(Gro_File,Topology_File,L,initial_molecules,solresnametop,Mixture='No', MD='Vacuum')
    dipole_model=Analysis.get_dipole_model() #I have discovered that this step stops later steps from being able to use -v
    print(dipole_model)
MD.run_md2(Topology_File,Mixture='No', MD='Box')        #Just because they gave us a box I don't trust that they have minimised it well.
######################################################################################################

#Here is the gaussian for pure liquids and mix
if user_settings["Mixture_Loop"]=='No':
    pure_solvent='Yes'
elif user_settings["Mixture_Loop"]=='Yes':
    if user_settings["Mode"] == "No_Gro" or user_settings["Mode"] == "Builds_Gro": 
        pure_solvent='Yes'
    elif user_settings["Mode"] == "Yes_Gro":
        pure_solvent='No'

R = advanced_settings["configuration"]["cluster_radius_nm"]
print(R)
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
max_jobs=State_Conditions["max_jobs"]
Gaussian_Location =State_Conditions["g09root"]
Scratch_Location =State_Conditions["GAUSS_SCRDIR"]
nproc=8
mem='5GB'

Gauss.nproc = nproc
Gauss.mem = mem
Gauss.max_jobs = max_jobs

Gauss.g09root= Gaussian_Location
Gauss.GAUSS_SCRDIR = Scratch_Location

#######################################################################
#Making oniom here now
print(f'topology file: {Topology_File}')
Topologys = Oniom_Generation.expand_topology_with_itps(Topology_File)
with open("Concat_Top.top", "w") as f:
    f.write(Topologys)


solvent_molecules=initial_molecules
Cut_Off_Radius=R
Oniom='oniom.inp'

f = open("Concat_Top.top")
Topology = f.read()
Oniom_Generation.Gen_File(Oniom,Configurations,solvent,Cut_Off_Radius)
Total_Atoms,qmax=Oniom_Generation.QM_Inputs(Topology,Oniom,qr1,qr2,qr3)     #We will need to add something here so it searches for the [atoms] for the solvent but I don't really know how to go about doing this script wise
Oniom_Generation.MM_Inputs(Topology,Oniom,qr1,qr2,qr3)  
Oniom_Generation.Counting_Molecules(Oniom,initial_molecules)
f.close()


print('tested up to here')  # Only testing to here as later steps will require simulations to actually be performed before we go much further.

hostname = socket.gethostname()
print(f"Running on host: {hostname}")
Gauss.natom=Total_Atoms
mu_Vacuum=Gauss.init0(sol_keyword,Gaus='Vacuum')
PCM1= Gauss.init1(sol_keyword,exp_diconst,Gaus='PCM1')
PCM2= Gauss.init2(sol_keyword,cal_diconst,Gaus='PCM2')

HOMEDIR = os.getcwd()
#Gaussian Process for pure liquids
if pure_solvent == 'Yes':
    create_dir_Simulations()
    rows = []
    for i in replicas:
        create_dir_reps(i)
        for T in T_list:
            create_dir_temps(T)
            for p in P_list:
                create_dir_press(p)
                mdfile='Junk.mdp'
                MD.run_md3(mdpfile,HOMEDIR,system_title,T,p,Mixture='No', MD='Production')
                exit_dir()
            exit_dir()
        exit_dir()
    exit_dir()
    
    dir_list = glob.glob('./Simulations/replica_*/*K/*.0Bar')
    count=0
    plot={}
    plot=pd.DataFrame(plot)
    for rundir in dir_list:
        os.chdir('./' + rundir)
        cwd = os.getcwd()
        print(cwd)
        Gauss.process_gro(exe=HOMEDIR +'/' + 'shellO',inp=HOMEDIR +'/' +'oniom.inp') 
        Gauss.init3(Gaus='SCEE_V0')
        SCEE=Gauss.init4(Gaus='SCEE_V1')
    
        delta_mu_list=[]
        mu_liquid_list=[]
        for i in df1['muL_SCEE']:
            mu_SCEE=i
            delta_mu=(mu_SCEE*(PCM1/PCM2))-mu_Vacuum
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
    
        df['replica'] = rundir  # or count if you prefer an integer label
        df['Simulation Density']=Analysis.Sim_Density()
        df['Liquid Dipole Moment Model'],df['Model Epsilon']=Analysis.get_dipole_model_liquid()
        plot = pd.concat([plot, df], ignore_index=True)
        
        count+=1
          
    print(f'Plotting distribution for {solvent}')
    IQR_K = 2.0
    Analysis.write_filtered_data(k=IQR_K)
    Analysis.plot_dipoledist(solvent,k=IQR_K)
    Analysis.plot_dipoledist_filtered(solvent,k=IQR_K)


    plot_filtered = read_filtered_data()
    x_filtered = plot_filtered['mu_liquid'].astype(float).to_numpy()
    x_filtered = x_filtered[np.isfinite(x_filtered)]
    n_filtered = len(x_filtered)

    collection={}
    collection = pd.DataFrame()
    collection['Solvent']                          = [f'{solvent}']
    collection['Temperatures']                     = T_list
    collection['Pressures']                        = P_list
    collection['Replicas']                         = [count]
    collection['Number of Configurations']         = [len(plot['mu_liquid'])]
    collection['Number of Configurations Filtered']= [n_filtered]

    collection['Scaling value 1']=[qr1]
    collection['Scaling value 2']=[qr2]
    collection['Scaling value 3']=[qr3]
    
    collection['Scaled Charge 1']=[num1] 
    collection['Scaled Charge  2']=[num2]
    collection['Scaled Charge  3']=[num3]
    
    collection['Experimental Density']=[density]
    collection['Simulation Density']=[plot['Simulation Density'].mean()]

    collection['Experimental Gas Dipole Moment']   = [gas_dipole]
    collection['Vacuum Dipole Moment Model']       = [dipole_model]
    collection['Liquid Dipole Moment Model']       = [plot['Liquid Dipole Moment Model'].mean()]
    collection['Experimental Refractive Index']    = [ref_ind]
    collection['Experimental Epsilon']             = [exp_diconst]
    collection['Model Epsilon']                    = [plot['Model Epsilon'].mean()]
    collection['cal_diconst']                      = [cal_diconst]
    collection['mu_Vacuum']                        = [mu_Vacuum[0]]
    collection['PCM1']                             = [PCM1[0]]
    collection['PCM2']                             = [PCM2[0]]
    
    # Unfiltered averages
    collection['muL_SCEE']                         = [plot['muL_SCEE'].mean()]
    collection['delta_mu']                         = [plot['delta_mu'].mean()]
    collection['mu_liquid']                        = [plot['mu_liquid'].mean()]

    x_all = plot['mu_liquid'].astype(float).to_numpy()
    finite_mask = np.isfinite(x_all)
    x_finite = x_all[finite_mask]
    keep, out, (q1, q3, iqr, lo, hi) = Analysis._iqr_masks(x_finite, k=IQR_K)
    full_keep = np.zeros(len(x_all), dtype=bool)
    full_keep[finite_mask] = keep
    plot_filtered = plot[full_keep].copy()
    
    replica_means = plot.groupby('replica')['mu_liquid'].mean()
    mu_mean_unfiltered = replica_means.mean()
    collection['STDEV'] = replica_means.std(ddof=1)  # ddof=1 for sample std
    collection['SE']  = collection['STDEV'] / np.sqrt(count)
    print(collection['mu_liquid'] > collection['Liquid Dipole Moment Model'])

    if collection['mu_liquid'][0] > Model_Dipole:
        collection['Dielectric Constant Correction'] = (collection['Experimental Refractive Index']**2 +(collection['mu_liquid'] /collection['Liquid Dipole Moment Model']) *     (collection['Model Epsilon'] + 1))
    else:
        collection['Dielectric Constant Correction'] = np.nan
    # Filtered averages
    collection['muL_SCEE_filtered']                = [plot_filtered['muL_SCEE'].astype(float).mean()]
    collection['delta_mu_filtered']                = [plot_filtered['delta_mu'].astype(float).mean()]
    collection['mu_liquid_filtered']               = [float(np.nanmean(x_filtered))]

    replica_means_filt = plot_filtered.groupby('replica')['mu_liquid'].apply(
        lambda x: np.nanmean(x.astype(float))
    )
    mu_mean_filtered = replica_means_filt.mean()
    collection['STDEV_filtered']  = replica_means_filt.std(ddof=1)
    collection['SE_filtered']  = collection['STDEV_filtered'] / np.sqrt(count)

    if collection['mu_liquid_filtered'][0] > Model_Dipole:
        collection['Dielectric Constant Correction Filtered'] = (collection['Experimental Refractive Index']**2 +(collection['mu_liquid_filtered'] /collection['Liquid Dipole Moment Model']) * (collection['Model Epsilon'] + 1))
    else:
        collection['Dielectric Constant Correction Filtered'] = np.nan

    collection.to_csv(f'../figures/vectors/Results_{solvent}.csv', index=False) 
    print(f'Data saved in ../figures/vectors/Results_{solvent}.csv')
    
