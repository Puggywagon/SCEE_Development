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

def setup_gromacs_env():
    home = os.path.expanduser("~")
    gromacs_root = os.path.join(home, "Software", "gromacs-2024.2")
    gromacs_bin = os.path.join(gromacs_root, "bin")

    if gromacs_bin not in os.environ["PATH"]:
        os.environ["PATH"] += os.pathsep + gromacs_bin

    if which("gmx") is None:
        raise EnvironmentError(f"'gmx' not found in PATH. Expected in: {gromacs_bin}")

    py_version = f"python{sys.version_info.major}.{sys.version_info.minor}"
    gmxapi_dir = os.path.join(gromacs_root, "lib", py_version, "site-packages")
    if os.path.isdir(gmxapi_dir) and gmxapi_dir not in sys.path:
        sys.path.insert(0, gmxapi_dir)

setup_gromacs_env()

try:
    import gromacs
except ImportError:
    try:
        import gmx as gromacs
    except ImportError:
        try:
            import gmxapi as gromacs
        except ImportError:
            raise ImportError("Could not import 'gromacs', 'gmx', or 'gmxapi'. Check your PYTHONPATH and installation.")


import Bash_steps
import Gro_Builder
import Pre_Eq_Simulations
import Oniom_Generation
import MD_Simulations
import SCEE
import Vacuum
import PCM1
import PCM2
import Analysis
import plot_rdf
import plot_dipoledist
import argparse
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
def read_data():
    csvfile_list = glob.glob(f'Simulations/replica_*/*/*/Dipole_Calculations.csv')
    df = pd.DataFrame()
    for i, csvfile in enumerate(csvfile_list, start=1):
        tmp = pd.read_csv(csvfile)
        df = pd.concat([df, tmp])
    return df
######################################################
def _iqr_masks(x, k=1.5):
    x = np.asarray(x, dtype=float)
    q1, q3 = np.percentile(x, [25, 75])
    iqr = q3 - q1
    lo = q1 - k * iqr
    hi = q3 + k * iqr
    out = (x < lo) | (x > hi)
    keep = ~out
    return keep, out, (q1, q3, iqr, lo, hi)
######################################################
def _normal_pdf(x, mu, sigma):
    return (1.0 / (sigma * np.sqrt(2.0 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
######################################################
def plot_dipoledist():
    parser = argparse.ArgumentParser()
    parser.add_argument('--show', type=str, default='True', help='plot figures')
    parser.add_argument('--bins', type=int, default=30, help='number of histogram bins')
    args = parser.parse_args()
    show = (args.show == 'True')
    df = read_data()
    x_all = df['mu_liquid'].astype(float).to_numpy()
    x_all = x_all[np.isfinite(x_all)]
    if x_all.size < 5:
        raise ValueError("Not enough dipole samples to form a distribution.")

    keep, out, (q1, q3, iqr, lo, hi) = _iqr_masks(x_all, k=1.5)
    x_keep = x_all[keep]
    x_out = x_all[out]
    
    mu = x_keep.mean()
    sigma = x_keep.std(ddof=0)
    if sigma == 0.0:
        raise ValueError("Sigma is zero after outlier removal. Check your data.")
    fig, ax = plt.subplots(figsize=(8, 6))
    # Binning: match your previous style or use full range
    # This matches your old approach starting at 1.00 D but updates the upper limit safely
    xmin = 1.00
    xmax = max(x_all.max(), 1.00) + 0.2
    bin_edges = np.linspace(xmin, xmax, args.bins + 1)
    bin_width = bin_edges[1] - bin_edges[0]
    # Histograms as counts (frequency)
    ax.hist(x_keep, bins=bin_edges, alpha=0.9, label="Data (outliers removed)")
    if x_out.size:
        ax.hist(x_out, bins=bin_edges, alpha=0.9, label="Outliers (IQR rule)")
    # Gaussian overlay scaled to counts
    x_line = np.linspace(bin_edges[0], bin_edges[-1], 400)
    y_line = _normal_pdf(x_line, mu, sigma) * x_keep.size * bin_width
    ax.plot(x_line, y_line, linewidth=2, label="Gaussian fit")

    ax.set_xlabel(r'$\mu_{liquid}$ / D', fontsize=16)
    ax.set_ylabel('Frequency', fontsize=16)

    # Keep these only if they are physically meaningful for this chapter
    ax.legend(frameon=False, fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=12)

    plt.tight_layout()
    plt.savefig(f'dipole_distribution.png', bbox_inches="tight", dpi=300)    
    plt.savefig(f'dipole_distribution.pdf', bbox_inches="tight", dpi=300)
    
    plt.close(fig)  
######################################################
Gro_Builder=Gro_Builder.Gro_Builder()
Bash_steps=Bash_steps.Bash_steps()
Pre_Eq_Simulations=Pre_Eq_Simulations.Pre_Eq_Simulations()
Oniom_Generation=Oniom_Generation.Oniom_Generation()
MD_Simulations=MD_Simulations.MD_Simulations()
Vacuum = Vacuum.Vacuum()
Molecule = SCEE.SCEE()
PCM1 = PCM1.PCM1()
PCM2 = PCM2.PCM2()
Analysis=Analysis.Analysis()

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

print('Does your system require separated .itp files? yes/no?')
#split=input()
split='yes'

print('What is your central and solvent molecues?')
print('Central molecules?')
#central=input()
central='Methylamine'

print('Its SMILE String?')
#mol = Chem.MolFromSmiles(input())
mol = 'CN'

print('Its residue name?')
#resname=input()
resname='MET'
solresnametop='Met'
print('Surrounding Solvent?')
#solvent=input()
solvent='Methylamine'

print('Its SMILE String?')
#solmol = Chem.MolFromSmiles(input())
solmol = 'CN'

print('Its residue name?')
#solresname=input()
solresname='MET'

if central == solvent:
    system_title=f'{central}'
else:
    system_title=f'{central}_in_{solvent}'

print('What is your solvent keyword for Gaussian?')
#sol_keyword=input()
sol_keyword='PropylAmine'

print('Does the central atom you are using require the application of a dummy atom?')
#dummy=input()
dummy='no'

print('What is the scaling ratio for your charges?')
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
gen_system='yes'


print('What is the experimental density of your?')
#density=input()
density=0.6562

print('What is the molar mass of your?')
#mol_mass=input()
mol_mass=31.1

print('What is the experimental dielectric constant of your molecule?')
#exp_diconst=input()
exp_diconst=4.9912

print('What is the experimental refractive index of your molecule?')
#ref_ind=input()
ref_ind=1.37

print('What is the experimental dipole moment of your molecule in the gas phase?')
#gas_dipole=input()
gas_dipole=1.31

print('How many replicas would you like to run?')
#reps=input()
replicas=range(1,2,1)
print('What temperatures (in K) would you like to run your simulations to run at? (please use 1, 2, ... format)')
#temps=input()
temps='298.0'
T_list=[temps]

print('What pressures (in K) would you like to run your simulations to run at? (please use 1, 2, ... format)')
#press=input()
press='1.0' 
P_list = [press]

f = open(f'{system_title}.top')
Topology = f.read()
f.close()

################################################################################
if split == 'yes':
    itp_list = Oniom_Generation.get_included_files(Topology)
    for filename in itp_list:
        if filename == f'{central}_AA.itp':
            f1 = open(f'{filename}')
            Solute = f1.read() 
            heavy_atoms=Oniom_Generation.Calc_Heavys(Solute,dummy,split)  
elif split == 'No':
    heavy_atoms=Oniom_Generation.Calc_Heavys(Solute,dummy,split)
cal_diconst=round(exp_diconst-(ref_ind**2)+1,3)
print(cal_diconst)

#Gro_Builder.AA_Structure(mol,system_title,resname)
#Gro_Builder.UA_Structure(solmol,system_title,solresname)

L=6.2
avo=6.022*10**23
N=(density*(((L*10**-7)**3))/(mol_mass))*avo
initial_molecules=int(np.ceil((N+0.05*N))-1)
print(initial_molecules)
#Pre_Eq_Simulations.Pre_Eq_Solute(central,L)
#Pre_Eq_Simulations.insert_Molecules(central,solvent,system_title,initial_molecules)
#Pre_Eq_Simulations.Write_to_top(system_title,initial_molecules,solresnametop)
#Pre_Eq_Simulations.Pre_Eq_System(system_title)
dipole_model=Pre_Eq_Simulations.get_dipole_model()
print(dipole_model)

HOMEDIR = os.getcwd()
################################################################################
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
Model_Dipole, Epsilon=Pre_Eq_Simulations.get_dipole_model_liquid(system_title)
ratio=Model_Dipole / qmax
num1=qr1*ratio*qmax
num2=qr2*ratio*qmax
num3=qr3*ratio*qmax
Molecule.natom=Total_Atoms
Vacuum.natom=Total_Atoms
PCM1.natom=Total_Atoms
PCM2.natom=Total_Atoms
################################################################################  
Vacuum.g09root='/home/zoe/Software/Gaussian/g09_pgi'
Vacuum.GAUSS_SCRDIR = '/home/zoe/Research/Gaussian/scratch'

Vacuum.rundir = 'opt-test_vacuum/'
Vacuum.natom=Total_Atoms

hostname = socket.gethostname()
print(f"Running on host: {hostname}")
if 'Tower3' in hostname:
    # Tower3: 48 cores, 125 GiB
    Vacuum.nproc = 8
    Vacuum.mem = '5GB'
    Vacuum.max_jobs = 4
else:
    # Home / office desktop: 24 cores, 15 GiB
    Vacuum.nproc = 6
    Vacuum.mem = '3GB'
    Vacuum.max_jobs = 3

#process *.gro files
#Vacuum.init_v0(system_title,sol_keyword)
#Vacuum.run_gaussian(step=0)
    
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
PCM1.g09root='/home/zoe/Software/Gaussian/g09_pgi/'
PCM1.GAUSS_SCRDIR = '/home/zoe/Research/Gaussian/scratch'
PCM1.rundir = 'opt-test-PCM1/'
PCM1.natom=Total_Atoms

hostname = socket.gethostname()
print(f"Running on host: {hostname}")
if 'Tower3' in hostname:
    # Tower3: 48 cores, 125 GiB
    PCM1.nproc = 8
    PCM1.mem = '5GB'
    PCM1.max_jobs = 4
else:
    # Home / office desktop: 24 cores, 15 GiB
    PCM1.nproc = 6
    PCM1.mem = '3GB'
    PCM1.max_jobs = 3
#PCM1
#PCM1.init_v0(system_title,sol_keyword,exp_diconst)
#PCM1.run_gaussian(step=0)
          
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
PCM2.natom=Total_Atoms

hostname = socket.gethostname()
print(f"Running on host: {hostname}")
if 'Tower3' in hostname:
    # Tower3: 48 cores, 125 GiB
    PCM2.nproc = 8
    PCM2.mem = '5GB'
    PCM2.max_jobs = 4

    # (Optionally match Vacuum/PCM1/PCM2 to something similar later)
else:
    # Home / office desktop: 24 cores, 15 GiB
    PCM2.nproc = 6
    PCM2.mem = '3GB'
    PCM2.max_jobs = 3
#PCM2 
#PCM2.init_v0(system_title,sol_keyword,cal_diconst)
#PCM2.run_gaussian(step=0)

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
dir_list = glob.glob('./Simulations/replica_1/*K/*.0Bar')
print(dir_list)

HOMEDIR = os.getcwd()
print(f'original directory: {HOMEDIR}')

Molecule.g09root='/home/zoe/Software/Gaussian/g09_pgi/'
Molecule.GAUSS_SCRDIR = '/home/zoe/Research/Gaussian/scratch'
Molecule.rundir = 'opt-test-SCEE/'
Molecule.natom=Total_Atoms
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
    Molecule.nproc = 6
    Molecule.mem = '3GB'
    Molecule.max_jobs = 3

count=0
plot={}
plot=pd.DataFrame(plot)
for rundir in dir_list:
    os.chdir('./' + rundir)
    cwd = os.getcwd()
    print(cwd)
    #process *.gro files
    shutil.copy2(HOMEDIR +'/' + 'oniom.inp', '.')
    shutil.copy2(HOMEDIR +'/' + 'shellO', '.')
    shutil.copy2(HOMEDIR +'/' + 'Shell_Oniom.f90', '.')
    
    #SCEE
    #Molecule.process_gro(exe=HOMEDIR +'/' + 'shellO',inp=HOMEDIR +'/' +'oniom.inp')
    #Molecule.init_v0()
    #Molecule.run_gaussian(step=0)
    #Molecule.init_v1()
    #Molecule.run_gaussian(step=1)
    df1 = Molecule.get_multipole_statistics()
    print(df1.describe())
    print(df1.head())
    muL_list = []
    for index, row in df1.iterrows():
        x = np.array([num1, num2, num3])
        y = np.array([row['dipole_l'], row['dipole_m'], row['dipole_h']])
        coeff = np.polyfit(x, y, 2)
        b = coeff[1] - 1
        fact = b*b-4*coeff[0]*coeff[2]
        if (coeff[0] >= 1.0e-6) and (fact > 0):
            muL = (-b - np.sqrt(fact)) / (2*coeff[0])
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
    os.chdir(HOMEDIR)  
    
    plot['muL_SCEE']=df1['muL_SCEE']
    plot['delta_mu'] = delta_mu_list
    plot['mu_liquid'] = mu_liquid_list

    count+=1
print(f'Plotting distribution for {central}')

plot_dipoledist(central)

collection={}
collection=pd.DataFrame(collection)
collection['Molecule']=f'{central}' #Input
collection['Temperatures']=T_list #Input
collection['Pressures']=P_list    #Input
collection['Replicas']=count  #counted in above iteration
collection['Number of Configurations']=len(plot['mu_liquid'])   #extracted from the total length of 

collection['Experimental Gas Dipole Moment']=gas_dipole
collection['Vacuum Dipole Moment Model']=dipole_model #Calculated from gromacs
collection['Liquid Dipole Moment Model']=Model_Dipole # calculated from gromacs
collection['Experimental Refractive Index']=ref_ind #Input

collection['Experimental Epsilon']=exp_diconst #Input
collection['Model Epsilon']=Epsilon #from gromacs at sane tune as liquid dipole moment
collection['cal_diconst']=cal_diconst #calculated
  
collection['mu_Vacuum']=mu_Vacuum #gaussian calc
collection['PCM1']=PCM1 #gaussian calc
collection['PCM2']=PCM2 #gaussian calc
collection['muL_SCEE']=plot['muL_SCEE'].mean()    #the configurations
collection['delta_mu']=plot['delta_mu'].mean()    
collection['mu_liquid']=plot['mu_liquid'].mean()  

collection['Dielectric Constant Correction'] = ref_ind**2 + (collection['mu_liquid'] / Model_Dipole) * (Epsilon + 1)
collection['STDEV']=np.std(plot['mu_liquid']) #Not sure if this is the correct approach for error calculations, need to figure our error by propogation for this instead

collection['SE']=collection['STDEV']/count  



################################################################################
