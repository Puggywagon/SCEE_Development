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
################################################################################ 
#These calculations aren't done in the best way, will need to adjust
Vacuum.init(sol_keyword)
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
#PCM1
PCM1.init(sol_keyword,exp_diconst)
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
###############################################################################
#PCM2 
PCM2.init(sol_keyword,cal_diconst)
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
                MD_Simulations.create_mdpfile(HOMEDIR,'junk.mdp',T,p)
                system_title=MD_Simulations.run_md('junk.mdp',HOMEDIR,System_Gro,Topology_File, pure_solvent)
                MD_Simulations.process_trajectory(system_title)
                MD_Simulations.process_gro()
                               
                Model_Dipole,epsilon=Pre_Eq_Simulations.get_dipole_model_liquid(system_title) 
                ratio=Model_Dipole / qmax
                num1=qr1*ratio*qmax
                num2=qr2*ratio*qmax
                num3=qr3*ratio*qmax
                
                
                Molecule.rundir = f'SCEE-{system_title}/'
                
                
               shutil.copy2(HOMEDIR +'/' + 'oniom.inp', '.')
               shutil.copy2(HOMEDIR +'/' + 'shellO', '.')
               shutil.copy2(HOMEDIR +'/' + 'Shell_Oniom.f90', '.')
    
               Molecule.process_gro(exe=HOMEDIR +'/' + 'shellO',inp=HOMEDIR +'/' +'oniom.inp') 
               Molecule.init_v0()
               Molecule.run_gaussian(step=0)
               Molecule.init_v1()
               Molecule.run_gaussian(step=1)
               df1 = Molecule.get_multipole_statistics()
               print(df1.describe())
               print(df1.head())
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
                
                values = df["mu_liquid"].dropna()
                n_cfg = values.count()

                mu_mean = values.mean()
                mu_stdev = values.std(ddof=1) if n_cfg > 1 else np.nan

                diconst_corr=Ref_Ind**2+((mu_mean /Model_Dipole)*(epsilon+1)

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
   
                exit_dir()
            exit_dir()
        exit_dir()
    exit_dir()
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
        Solute_Gro_File=Mixture_Loop["Solute_Gro_File"]
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
