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
import argparse

pipe = Popen("/home/zoe/Software/gromacs-2024.2/bin/GMXRC.bash; env", stdout=PIPE, \
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
def read_data():
    csvfile_list = glob.glob(f'Simulations/replica_*/*/*/Dipole_Calculations.csv')
    df = pd.DataFrame()
    for i, csvfile in enumerate(csvfile_list, start=1):
        tmp = pd.read_csv(csvfile)
        df = pd.concat([df, tmp])
    return df
######################################################
def write_filtered_data(k):
    csvfile_list = glob.glob(f'Simulations/replica_*/*/*/Dipole_Calculations.csv')
    for csvfile in csvfile_list:
        df = pd.read_csv(csvfile)
        x = df['mu_liquid'].astype(float).to_numpy()
        finite_mask = np.isfinite(x)
        x_finite = x[finite_mask]
        if x_finite.size < 5:
            continue
        keep, out, (q1, q3, iqr, lo, hi) = _iqr_masks(x_finite, k=k)
        full_keep = np.zeros(len(x), dtype=bool)
        full_keep[finite_mask] = keep
        df_filtered = df.copy()
        cols_to_blank = ['muL_SCEE', 'delta_mu', 'mu_liquid']
        df_filtered.loc[~full_keep, cols_to_blank] = np.nan
        out_dir = os.path.dirname(csvfile)
        out_path = os.path.join(out_dir, 'Dipole_Calculations_filtered.csv')
        df_filtered.to_csv(out_path, index=False)
######################################################        
def _iqr_masks(x, k):
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
def plot_dipoledist(Molecule,k, bin_width=0.25, x_ceil=10.0):
    parser = argparse.ArgumentParser()
    parser.add_argument('--show', type=str, default='True', help='plot figures')
    args = parser.parse_args()
    show = (args.show == 'True')

    df = read_data()
    x_all = df['mu_liquid'].astype(float).to_numpy()
    x_all = x_all[np.isfinite(x_all)]
    if x_all.size < 5:
        raise ValueError("Not enough dipole samples to form a distribution.")

    keep, out, (q1, q3, iqr, lo, hi) = _iqr_masks(x_all, k=k)
    x_keep = x_all[keep]
    x_out  = x_all[out]

    mu    = x_keep.mean()
    sigma = x_keep.std(ddof=0)
    if sigma == 0.0:
        raise ValueError("Sigma is zero after outlier removal. Check your data.")

    # Axis limits
    xmin = 0.0
    largest_populated = x_all.max()  # will be overridden below
    # find the largest bin edge that actually has data within x_ceil
    x_within_ceil = x_all[x_all <= x_ceil]
    if x_within_ceil.size > 0:
        largest_populated = np.ceil(x_within_ceil.max() / bin_width) * bin_width
    xmax = min(x_all.max(), x_ceil)
    xmax = np.ceil(xmax / bin_width) * bin_width

    bin_edges = np.arange(xmin, xmax + bin_width, bin_width)
    bw = bin_edges[1] - bin_edges[0]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(x_keep, bins=bin_edges, alpha=0.9, color='purple', label="Data")
    if x_out.size:
        # only plot outliers that fall within the axis range
        x_out_visible = x_out[(x_out >= xmin) & (x_out <= xmax)]
        if x_out_visible.size:
            ax.hist(x_out_visible, bins=bin_edges, alpha=0.9,color='orange', label="Outliers (IQR rule)")

    x_line = np.linspace(xmin - bin_width, xmax + bin_width, 400)  # already correct
    y_line = _normal_pdf(x_line, mu, sigma) * x_keep.size * bw
    ax.plot(x_line, y_line, color='#008080', linewidth=2, label="Gaussian fit")
    ax.set_xlim(xmin, xmax)
    raw_spacing = (xmax - xmin) / 8                  # replace the old tick_spacing logic
    tick_spacing = round(raw_spacing * 2) / 2
    tick_spacing = max(tick_spacing, 0.5)
    ax.set_xticks(np.arange(xmin, xmax + tick_spacing, tick_spacing))

    ax.set_xlabel(r'$\mu_{liquid}$ / D', fontsize=16)
    ax.set_ylabel('Frequency', fontsize=16)
    ax.legend(frameon=False, fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.tight_layout()
    plt.savefig(f'../../../../../../../Research/SCEE_Project/Results/Aromatics_Chapter/figures/dipole_distribution_{Molecule}.png', bbox_inches="tight", dpi=300)
    plt.savefig(f'../../../../../../../Research/SCEE_Project/Results/Aromatics_Chapter/figures/dipole_distribution_{Molecule}.pdf', bbox_inches="tight", dpi=300)
    plt.close(fig)

    # Summary statistics table
    stats = {
        'Molecule':         [Molecule],
        'N_total':          [x_all.size],
        'N_kept':           [x_keep.size],
        'N_outliers':       [x_out.size],
        'mean_D':           [round(mu, 4)],
        'std_D':            [round(sigma, 4)],
        'median_D':         [round(np.median(x_keep), 4)],
        'IQR_lo':           [round(lo, 4)],
        'IQR_hi':           [round(hi, 4)],
        'q1':               [round(q1, 4)],
        'q3':               [round(q3, 4)],
    }
    stats_df = pd.DataFrame(stats)
    stats_df.to_csv(f'../../../../../../../Research/SCEE_Project/Results/Aromatics_Chapter/figures/vectors/dipole_stats_{Molecule}.csv', index=False)
    print(stats_df.to_string(index=False))
######################################################
def read_filtered_data():
    csvfile_list = glob.glob(f'Simulations/replica_*/*/*/Dipole_Calculations_filtered.csv')
    df = pd.DataFrame()
    for csvfile in csvfile_list:
        tmp = pd.read_csv(csvfile)
        df = pd.concat([df, tmp])
    return df
######################################################
def plot_dipoledist_filtered(Molecule, k, bin_width=0.25, x_ceil=10.0):
    parser = argparse.ArgumentParser()
    parser.add_argument('--show', type=str, default='True', help='plot figures')
    args = parser.parse_args()

    df_unfiltered = read_data()
    x_all = df_unfiltered['mu_liquid'].astype(float).to_numpy()
    x_all = x_all[np.isfinite(x_all)]
    if x_all.size < 5:
        raise ValueError("Not enough dipole samples to form a distribution.")

    keep, out, (q1, q3, iqr, lo, hi) = _iqr_masks(x_all, k=k)
    x_keep = x_all[keep]
    x_out  = x_all[out]

    mu    = x_keep.mean()
    sigma = x_keep.std(ddof=0)
    if sigma == 0.0:
        raise ValueError("Sigma is zero. Check your data.")

    # Gaussian threshold for axis limits
    gauss_threshold = 0.01
    x_test = np.linspace(0.0, x_ceil * 2, 10000)
    gauss_test = _normal_pdf(x_test, mu, sigma)
    gauss_peak = gauss_test.max()
    above_thresh = x_test[gauss_test >= gauss_threshold * gauss_peak]
    xmin_gauss = max(0.0, np.floor(above_thresh.min() / bin_width) * bin_width)
    xmax_gauss = min(x_ceil, np.ceil(above_thresh.max() / bin_width) * bin_width)

    # Gap detection to avoid long empty tails
    counts, edges = np.histogram(x_all[x_all <= x_ceil],
                                 bins=np.arange(0, x_ceil + bin_width, bin_width))
    nonzero_bins = np.where(counts > 0)[0]
    last_populated = edges[nonzero_bins[-1] + 1]
    gap_cut = None
    for i in range(len(nonzero_bins) - 1):
        gap = nonzero_bins[i+1] - nonzero_bins[i]
        if gap > 3:
            gap_cut = edges[nonzero_bins[i] + 1]
            break
    if gap_cut is not None:
        xmax_data = np.ceil((gap_cut + 0.5) / bin_width) * bin_width
    else:
        xmax_data = np.ceil((last_populated + 0.5) / bin_width) * bin_width

    xmax = min(x_ceil, max(xmax_gauss, xmax_data))

    # Extend to show mild outliers within 2 sigma
    mild_out = x_out[np.abs(x_out - mu) <= 2 * sigma]
    if mild_out.size:
        xmin_mild = np.floor(mild_out.min() / bin_width) * bin_width
        xmax_mild = np.ceil(mild_out.max() / bin_width) * bin_width
        xmin = max(0.0, min(xmin_gauss, xmin_mild))
        xmax = min(x_ceil, max(xmax, xmax_mild))
    else:
        xmin = xmin_gauss

    bin_edges = np.arange(xmin, xmax + bin_width, bin_width)
    bw = bin_edges[1] - bin_edges[0]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(x_keep, bins=bin_edges, alpha=0.9,color='purple', label="Data (outliers removed)")
    if x_out.size:
        x_out_visible = x_out[(x_out >= xmin) & (x_out <= xmax)]
        if x_out_visible.size:
            ax.hist(x_out_visible, bins=bin_edges, alpha=0.9,color='orange', label="Outliers (IQR rule)")

    x_line = np.linspace(xmin - bin_width, xmax + bin_width, 400)
    y_line = _normal_pdf(x_line, mu, sigma) * x_keep.size * bw
    ax.plot(x_line, y_line, color='#008080', linewidth=2, label="Gaussian fit")

    ax.set_xlim(xmin, xmax)
    raw_spacing = (xmax - xmin) / 8
    tick_spacing = round(raw_spacing * 2) / 2
    tick_spacing = max(tick_spacing, 0.5)
    ax.set_xticks(np.arange(xmin, xmax + tick_spacing, tick_spacing))
    ax.set_xlabel(r'$\mu_{liquid}$ / D', fontsize=16)
    ax.set_ylabel('Frequency', fontsize=16)
    ax.legend(frameon=False, fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.tight_layout()
    plt.savefig(f'../../../../../../../Research/SCEE_Project/Results/Aromatics_Chapter/figures/dipole_distribution_{Molecule}_filtered.png', bbox_inches="tight", dpi=300)
    plt.savefig(f'../../../../../../../Research/SCEE_Project/Results/Aromatics_Chapter/figures/dipole_distribution_{Molecule}_filtered.pdf', bbox_inches="tight", dpi=300)
    plt.close(fig)

    stats = {
        'Molecule':   [Molecule],
        'N_total':    [x_all.size],
        'N_kept':     [x_keep.size],
        'N_outliers': [x_out.size],
        'mean_D':     [round(mu, 4)],
        'std_D':      [round(sigma, 4)],
        'median_D':   [round(np.median(x_keep), 4)],
        'IQR_lo':     [round(lo, 4)],
        'IQR_hi':     [round(hi, 4)],
        'q1':         [round(q1, 4)],
        'q3':         [round(q3, 4)],
    }
    stats_df = pd.DataFrame(stats)
    stats_df.to_csv(f'../../../../../../../Research/SCEE_Project/Results/Aromatics_Chapter/figures/vectors/dipole_stats_{Molecule}_filtered.csv', index=False)
    print(stats_df.to_string(index=False))
######################################################
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
#Vacuum.init_v0(system_title,sol_keyword)
#Vacuum.run_gaussian(step=0)
    
df4 = Vacuum.get_multipole_statistics()
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
#PCM1.init_v0(system_title,sol_keyword,exp_diconst)
#PCM1.run_gaussian(step=0)
          
df2 = PCM1.get_multipole_statistics()
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
#PCM2.init_v0(system_title,sol_keyword,cal_diconst)
#PCM2.run_gaussian(step=0)

df3 = PCM2.get_multipole_statistics()
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


#Note this doesn't work anymore
#cmd=["f95 Shell_Oniom.f90 -o shellO"]
#subprocess.call(cmd)

count=0
plot = pd.DataFrame(columns=["muL_SCEE", "delta_mu", "mu_liquid"])
for rundir in dir_list:
    os.chdir('./' + rundir)
    cwd = os.getcwd()
    print(cwd)
    #process *.gro files
    #shutil.copy2(HOMEDIR +'/' + 'oniom.inp', '.')
    #shutil.copy2(HOMEDIR +'/' + 'shellO', '.')
    #shutil.copy2(HOMEDIR +'/' + 'Shell_Oniom.f90', '.')
    
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
    
    df['replica'] = rundir  # or count if you prefer an integer label
    
    plot = pd.concat([plot, df], ignore_index=True)

    count+=1


print(f'Plotting distribution for {central}')
IQR_K = 2.0
write_filtered_data(k=IQR_K)
plot_dipoledist(central,k=IQR_K)
plot_dipoledist_filtered(central,k=IQR_K)


plot_filtered = read_filtered_data()
x_filtered = plot_filtered['mu_liquid'].astype(float).to_numpy()
x_filtered = x_filtered[np.isfinite(x_filtered)]
n_filtered = len(x_filtered)

collection={}
collection = pd.DataFrame()
collection['Solvent']                          = [f'{central}']
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

collection['Experimental Gas Dipole Moment']   = [gas_dipole]
collection['Vacuum Dipole Moment Model']       = [dipole_model]
collection['Liquid Dipole Moment Model']       = [Model_Dipole]
collection['Experimental Refractive Index']    = [ref_ind]
collection['Experimental Epsilon']             = [exp_diconst]
collection['Model Epsilon']                    = [Epsilon]
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
keep, out, (q1, q3, iqr, lo, hi) = _iqr_masks(x_finite, k=IQR_K)
full_keep = np.zeros(len(x_all), dtype=bool)
full_keep[finite_mask] = keep
plot_filtered = plot[full_keep].copy()

replica_means = plot.groupby('replica')['mu_liquid'].mean()
mu_mean_unfiltered = replica_means.mean()
collection['STDEV'] = replica_means.std(ddof=1)  # ddof=1 for sample std
collection['SE']  = collection['STDEV'] / np.sqrt(count)
print(collection['mu_liquid'] > collection['Liquid Dipole Moment Model'])

if collection['mu_liquid'][0] > Model_Dipole:
    collection['Dielectric Constant Correction'] = (collection['Experimental Refractive Index']**2 +(collection['mu_liquid'] /collection['Liquid Dipole Moment Model']) * (collection['Model Epsilon'] + 1))
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

collection.to_csv(f'../../../../../../../Research/SCEE_Project/Results/Aromatics_Chapter/figures/vectors/Results_{central}.csv', index=False) 
print(f'Data saved in ~/Research/SCEE_Project/Results/Aromatics_Chapter/figures/vectors/Results_{central}.csv')
################################################################################  
