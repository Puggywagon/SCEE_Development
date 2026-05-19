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
import gromacs
import itertools
from contextlib import contextmanager

import Gro_Builder
import Gro_Simulations
import Simulations_Analysis
import Oniom_Generation
import Gaussian_Calculations
from Settings_Loader import load_settings

################################################################################
def dir_System(System_Title):
    dire=f'{System_Title}'
    os.chdir(dire)
################################################################################
@contextmanager
def cd(path):
    """Change directory for the duration of the block, restore on exit.
    Creates the directory if it doesn't exist.
    """
    prev = os.getcwd()
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev)
################################################################################
def confirm_structure():
    print('Please take a minute to look at the molecule you have just '
          'constructed using RD Kit. Confirm with the image generated if '
          'it resembles the molecule you are interested in studying. '
          'Also confirm your topology file matches the molecule gro files '
          'generated. Once ready to continue type "Yes", if you need to '
          'restart the process type "No".')
    # ready = input()
    ready = 'Yes'
    if ready == 'Yes':
        print('The script will continue on with the next steps of the system.')
    elif ready == 'No':
        print('Please restart this script after you have made your adjustments.')
        exit(0)
######################################################################################################
print('To start the process, you need to complete the Yaml script based on your system requirements. If you have already completed this please type Yes, If you need to go back and finish doing this please respond No.')
#ready=input()
ready='Yes'
if ready == 'Yes':
    print('This script will continue with the process.')
elif ready == 'No':
    print('Please restart this script once you are ready to start.')
    exit(0)


settings = load_settings()
Gro_Builder=Gro_Builder.Gro_Builder()
Gauss = Gaussian_Calculations.Gaussian_Calculations(settings)
MD = Gro_Simulations.Gro_Simulations()
Analysis=Simulations_Analysis.Simulations_Analysis()
Oniom_Generation=Oniom_Generation.Oniom_Generation() 

dir_System(settings.molecule.system_title)

MD.ntmpi = settings.software.gromacs.ntmpi
MD.ntomp = settings.software.gromacs.ntomp
MD.pin = settings.software.gromacs.pin
MD.pinoffset = settings.software.gromacs.pinoffset

if settings.build_solvent_gro:
    Gro_Builder.AA_Solvent(settings.molecule.solvent_smiles)
    if not settings.molecule.aa_surround:
        Gro_Builder.UA_Solvent(settings.molecule.solvent_smiles)
confirm_structure()
MD.run_md1(settings, role='pure')
MD.run_md2(settings, role='pure')        
######################################################################################################
#Making oniom here now
Oniom = 'oniom.inp'
Oniom_Generation.Gen_File(
    Oniom,
    settings.advanced.sampling.n_configurations,
    settings.molecule.system_title,
    settings.advanced.configuration.cluster_radius_nm,
)
Total_Atoms, qmax = Oniom_Generation.QM_Inputs(settings, role='pure', Oniom=Oniom)
Oniom_Generation.MM_Inputs(settings, Oniom=Oniom)
Oniom_Generation.Counting_Molecules(Oniom, settings.initial_molecules)
Gauss.natom = Total_Atoms
Gauss.qmax = qmax

hostname = socket.gethostname()
print(f"Running on host: {hostname}")
mu_Vacuum=Gauss.init0()
PCM1= Gauss.init1()
PCM2= Gauss.init2()

dipole_model=Analysis.get_dipole_model() #I have discovered that this step stops later steps from being able to use -v
gas_potential=Analysis.Pot_Gas() #I have discovered that this step stops later steps from being able to use -v
print(dipole_model)

HOMEDIR = os.getcwd()
replicas = range(1, settings.state.replicas + 1)
all_dipole_dfs = []
leaf_summary_rows = []
for replica, T, p in itertools.product(replicas, settings.state.temperatures, settings.state.pressures):
    leaf = f'Simulations/replica_{replica}/{T:g}K/{p:g}Bar'
    with cd(leaf):
        MD.run_md3(settings, role='pure', HOMEDIR=HOMEDIR, T=T, p=p)
        # Gaussian processing on the 200 configurations
        Gauss.process_gro(oniom_inp_path=f'{HOMEDIR}/oniom.inp')
        Gauss.init3()
        SCEE = Gauss.init4()
        
        dipole_df = Analysis.process_dipole(SCEE, mu_Vacuum, PCM1, PCM2, settings)
        sim_density = Analysis.Sim_Density()
        liq_dipole, model_eps = Analysis.get_dipole_model_liquid()
        liq_potential = Analysis.Pot_Liq()
        enth_of_vap = Analysis.Sim_Enth(gas_potential,liq_potential, T,N=settings.initial_molecules)
        
        
    dipole_df['replica'] = replica
    dipole_df['T'] = T
    dipole_df['p'] = p
    all_dipole_dfs.append(dipole_df)
    leaf_summary_rows.append({
        'replica': replica, 
        'T': T, 
        'p': p,
        'Simulation Density': sim_density,
        'Liquid Dipole Moment Model': liq_dipole,
        'Model Epsilon': model_eps,
        'Liquid Potential': liq_potential,
        'Enthalpy of Vaporisation': enth_of_vap,    
    })    
    

dipole_results = pd.concat(all_dipole_dfs, ignore_index=True)
leaf_summary = pd.DataFrame(leaf_summary_rows)    

print(f'Filtering and plotting per-(T,p) distributions for {settings.molecule.system_title}')
IQR_K = 2.0

# Delete dipole_stats.csv so we start fresh on each run
if os.path.exists('dipole_stats.csv'):
    os.remove('dipole_stats.csv')

keep_masks = {}  # (T, p) -> bool array into the pool at that (T,p)
for T, p in itertools.product(settings.state.temperatures, settings.state.pressures):
    df_pool = dipole_results[
        (dipole_results['T'] == T) & (dipole_results['p'] == p)
    ].reset_index(drop=True)
    keep_masks[(T, p)] = Analysis.plot_state_point(
        df_pool, T, p, mu_Vacuum=mu_Vacuum, k=IQR_K,
    )

# Write per-leaf Dipole_Calculations_filtered.csv files using pooled keep decisions
Analysis.write_filtered_per_leaf(dipole_results, keep_masks)

per_replica_df, system_df = Analysis.build_collection(
    dipole_results=dipole_results,
    leaf_summary=leaf_summary,
    keep_masks=keep_masks,
    settings=settings,
    mu_Vacuum=mu_Vacuum,
    PCM1=PCM1,
    PCM2=PCM2,
    dipole_model=dipole_model,
    gas_potential=gas_potential,
)

per_replica_df.to_csv('Per_Replica.csv', index=False)
system_df.to_csv('Results.csv', index=False)
print(f'Per-replica data saved in Per_Replica.csv')
print(f'System data saved in Results.csv')
                
if settings.is_mixture:
    if settings.build_solute_gro:
        Gro_Builder.AA_Solute(settings.molecule.solute_smiles)

    confirm_structure()
    
    MD.run_md1(settings, role='mixture')
    MD.run_md2(settings, role='mixture')    
    
    replicas = range(1, settings.state.replicas + 1)
    all_results = []
    for replica, T, p in itertools.product(replicas, settings.state.temperatures, settings.state.pressures):
        leaf = f'Simulations/replica_{replica}/{T}K/{p}Bar'
        with cd(leaf):
            MD.run_md3(settings, role='mixture', HOMEDIR=HOMEDIR, T=T, p=p)
