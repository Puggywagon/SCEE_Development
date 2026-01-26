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
MD=Gro_Simulations.Gro_Simulations()
Analysis=Simulations_Analysis.Simulations_Analysis() #Make this into an analysis script rather than what I have it now.
Oniom_Generation=Oniom_Generation.Oniom_Generation()
Gauss=Gaussian_Calculations.Gaussian_Calculations()

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
    MD.run_md(self, MD='Vacuum',Gro_File,Topology_File,L,initial_molecules,solresnametop)
    dipole_model=Analysis.get_dipole_model() #I have discovered that this step stops later steps from being able to use -v
    print(dipole_model)
MD.run_md(self, MD='Box',Topology_File)        #Just because they gave us a box I don't trust that they have minimised it well.
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
Gauss.natom=Total_Atoms
Vacuum.natom=Total_Atoms
PCM1.natom=Total_Atoms
PCM2.natom=Total_Atoms

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
################################################################################ 
hostname = socket.gethostname()
print(f"Running on host: {hostname}")

mu_Vacuum=Gauss.init(Gaus='Vacuum',sol_keyword)
PCM1= Gauss.init(Gaus='PCM1',sol_keyword,exp_diconst)
PCM2= Gauss.init(Gaus='PCM2',sol_keyword,cal_diconst)
#######################################################################
#Making oniom here now
f = open(Topology_File)
Topology = f.read()
f.close()

solvent_molecules=initial_molecules
Cut_Off_Radius=R
Oniom='oniom.inp'

Oniom_Generation.Gen_File(Oniom,Configurations, Cut_Off_Radius)
itp_list = Oniom_Generation.get_included_files(Topology)
Total_Atoms,qmax=Oniom_Generation.QM_Inputs(Solvent,Topology,qr1,qr2,qr3)     #We will need to add something here so it searches for the [atoms] for the solvent but I don't really know how to go about doing this script wise
Oniom_Generation.MM_Inputs(Solvent,Topology,Oniom,qr1,qr2,qr3)  
Oniom_Generation.Counting_Molecules(Solvent,Oniom,initial_molecules)

hostname = socket.gethostname()
print(f"Running on host: {hostname}")
Gauss.natom=Total_Atoms
mu_Vacuum=Gauss.init(Gaus='Vacuum',sol_keyword)
PCM1= Gauss.init(Gaus='PCM1',sol_keyword,exp_diconst)
PCM2= Gauss.init(Gaus='PCM2',sol_keyword,cal_diconst)

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
                MD.run_md(self, MD='Production',mdpfile,HOMEDIR,Topology_File,T,p)
                exit_dir()
            exit_dir()
        exit_dir()
    exit_dir()
    
    dir_list = glob.glob('./Simulations/replica_*/*K/*.0Bar')

    for rundir in dir_list:
        # Here we copy the files into directory we are running these steps
        os.chdir('./' + rundir)
        cwd = os.getcwd()
        print(cwd)
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
        
        values = df["mu_liquid"].dropna()
        n_cfg = values.count()

        mu_mean = values.mean()
        mu_stdev = values.std(ddof=1) if n_cfg > 1 else np.nan

        diconst_corr=(Ref_Ind**2)+((mu_mean /Model_Dipole)*(epsilon+1)

        row = {
            "Replica": i,
            "Temperature": T,
            "Pressure": p,
            "muL_SCEE": df["muL_SCEE"].mean(),
            "delta_mu": df["delta_mu"].mean(),
            "mu_liquid": mu_mean,
            "mu_stdev": mu_stdev,
            }
                
        diconst_corr = Ref_Ind**2 + (mu_mean / Model_Dipole) * (epsilon + 1)  #Can we check this???
        row["diconst_corr"] = diconst_corr
        mu_se = (mu_stdev / np.sqrt(n_cfg)) if n_cfg > 1 else np.nan
        row["mu_se"] = mu_se
        df.to_csv('Dipole_Calculations.csv', index=False)  
        os.chdir(HOMEDIR)  
        
    df5 = pd.DataFrame(rows)
    df5.to_csv(HOMEDIR + "Results_1.csv", index=False)
    df6=pd.read_csv('Results_1.csv')
    df6.columns=['Replica','Temperature','Pressure','muL_SCEE','delta_mu','mu_liquid','mu_stdev','diconst_corr']
    df6 = pd.read_csv(HOMEDIR + "Results_1.csv")
    
    df7 = (
        df6.groupby(["Temperature", "Pressure"], as_index=False)
           .apply(combine_group)
           .reset_index(drop=True))

    df7.to_csv(HOMEDIR + "Results_2.csv", index=False)
    



#Once you have reached this point you can stop as this is my mixture ideas.

#Gaussian process for mixtures
if user_settings["Mixture_Loop"]=='Yes':
    if user_settings["Mode"] == "Yes_Gro":
        Yes_Gro=user_settings["Yes_Gro"]
        Mixture_Loop=Yes_Gro["Mixture_Settings"]
        Solute=Mixture_Loop["solute"]
        
        Box_Build='No'
        
    elif user_settings["Mode"] == "No_Gro":
        No_Gro=user_settings["No_Gro"]
        Mixture_Loop=No_Gro["Mixture_Settings"]
        Solute=Mixture_Loop["solute"]
   
   
                exit_dir()
            exit_dir()
        exit_dir()
    exit_dir()     Solute_Gro_File=Mixture_Loop["Solute_Gro_File"]
        Solute_Topology_File=Mixture_Loop["Solute_Topology_File"]
        resnametop=Mixture_Loop["resnametop"]
        
        Box_Build='Yes'
        
    elif user_settings["Mode"] == "Builds_Gro":
        Builds_Gro=user_settings["Builds_Gro"]
        Mixture_Loop=Builds_Gro["Mixture_Settings"]
        Solute=Mixture_Loop["solute"]
        Solute_Gro_File=Mixture_Loop["Solute_Gro_File"]
        Solute_Topology_File=Mixture_Loop["Solute_Topology_File"]
        resnametop=Mixture_Loop["resnametop"]
        resname=Mixture_Loop["resname"]
        mol=Mixture_Loop["mol"]
        
        Gro_File=Gro_Builder.AA_Structure(mol,solute,resname) #solute #Rename this to handle UA and AA models rather than how I have been doing things..., get this to return
        Box_Build='Yes'
        

    if Box_Build == 'Yes'
    
        L = advanced_settings["configuration"]["box_length_nm"]
        solute_molecules=1
    
        Pre_Eq_Simulations.Pre_Eq_Solute(Sol_Gro_File,Solute_Topology_File,L)
        #dipole_model=Pre_Eq_Simulations.get_dipole_model() # I dunno if we need this to do the dielectric constant values
        Pre_Eq_Simulations.insert_Molecules(System_Gro,Sol_Gro_File,solute_molecules) 
        Pre_Eq_Simulations.Write_to_top(Topology_File,solute_molecules,resnametop) 
    
    Mixture_System_Gro=Pre_Eq_Simulations.Pre_Eq_System(Topology_File)        #Just because they gave us a box I don't trust that they have minimised it well.

    #Same loop for MD and use mixture gaussian 
