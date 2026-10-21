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
from datetime import datetime

import Gro_Builder
import Gro_Simulations
import Simulations_Analysis
import Oniom_Generation
import Gaussian_Calculations
from Settings_Loader import load_settings

################################################################################
def dir_System(System_Title):
    dire = f'{System_Title}'
    os.chdir(dire)
################################################################################
def save_settings_snapshot(settings, source_path):
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    dest = f'./Results/settings_{settings.molecule.system_title}_{timestamp}.yml'
    shutil.copy(source_path, dest)
    print(f'Settings snapshot saved as {dest}')
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
print('To start the process, you need to complete the Yaml script based on your system requirements. If you have already completed this please type Yes, If you need to go back and finish doing this please respond No.')
#ready=input()
ready='Yes'
if ready == 'Yes':
    print('This script will continue with the process.')
elif ready == 'No':
    print('Please restart this script once you are ready to start.')
    exit(0)


settings = load_settings()
SETTINGS_PATH = os.path.abspath('Settings.yml')

# Which stages to run: set in the Run_Control block of Settings.yml
RUN_BUILD      = settings.run.build
RUN_MD_SETUP   = settings.run.md_setup
RUN_PRODUCTION = settings.run.production
RUN_SCEE_QM    = settings.run.scee_qm
RUN_GAS_QM     = settings.run.gas_qm
print(f'Run control: build={RUN_BUILD}, md_setup={RUN_MD_SETUP}, '
      f'production={RUN_PRODUCTION}, scee_qm={RUN_SCEE_QM}, gas_qm={RUN_GAS_QM}')

Gro_Builder=Gro_Builder.Gro_Builder()
Gauss = Gaussian_Calculations.Gaussian_Calculations(settings)
MD = Gro_Simulations.Gro_Simulations()
Analysis=Simulations_Analysis.Simulations_Analysis()
Oniom_Generation=Oniom_Generation.Oniom_Generation() 

dir_System(settings.molecule.system_title)
save_settings_snapshot(settings, SETTINGS_PATH)

MD.ntmpi = settings.software.gromacs.ntmpi
MD.ntomp = settings.software.gromacs.ntomp
MD.pin = settings.software.gromacs.pin
MD.pinoffset = settings.software.gromacs.pinoffset

if RUN_BUILD and settings.build_solvent_gro:
    Gro_Builder.AA_Solvent(settings.molecule.solvent_smiles)
    if not settings.molecule.aa_surround:
        Gro_Builder.UA_Solvent(settings.molecule.solvent_smiles)
    confirm_structure()
if RUN_MD_SETUP:
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

# Gas-phase QM: Vacuum, PCM1 and PCM2 on configurations from the vacuum NVT run,
# each optimised/calculated at v0 then single-pointed at v1, then averaged.
# With RUN_GAS_QM = False, existing gas-phase .out files are re-read instead.
run_gas_qm = RUN_GAS_QM
if RUN_GAS_QM:
    gas_gros = MD.extract_vacuum_configurations(settings, role='pure')
else:
    # Re-use the configurations from the previous extraction (or the final frame)
    gas_gros = sorted(glob.glob('Gas_Configs/gas_conf_*.gro'),
                      key=lambda s: int(s.split('_')[-1].replace('.gro', ''))) or ['nvt_vacuum2.gro']
print(f'Using {len(gas_gros)} gas-phase configuration(s)')
mu_Vacuum = Gauss.init0(gas_gros, run=run_gas_qm)
PCM1 = Gauss.init1(gas_gros, run=run_gas_qm)
PCM2 = Gauss.init2(gas_gros, run=run_gas_qm)

dipole_model=Analysis.get_dipole_model() #I have discovered that this step stops later steps from being able to use -v
gas_potential=Analysis.Pot_Gas() #I have discovered that this step stops later steps from being able to use -v

HOMEDIR = os.getcwd()
replicas = range(1, settings.state.replicas + 1)
all_dipole_dfs = []
leaf_summary_rows = []
for replica, T, p in itertools.product(replicas, settings.state.temperatures, settings.state.pressures):
    leaf = Analysis.leaf_path(replica, T, p)
    if not RUN_PRODUCTION and not os.path.isdir(leaf):
        raise FileNotFoundError(f'Expected existing run folder {leaf} (reanalysis mode)')
    with cd(leaf):
        if RUN_PRODUCTION:
            MD.run_md3(settings, role='pure', HOMEDIR=HOMEDIR, T=T, p=p)
        if RUN_SCEE_QM:
            # Build the 3 x 200 ONIOM inputs and run the v0 optimisations
            Gauss.process_gro(oniom_inp_path=f'{HOMEDIR}/oniom.inp')
            Gauss.init3(run=True)
        # v1 single points (run only if RUN_SCEE_QM), then read all results
        scee_df = Gauss.init4(run=RUN_SCEE_QM)
        
        dipole_df = Analysis.process_dipole(scee_df, mu_Vacuum, PCM1, PCM2, settings)
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

# Delete Dipole_Stats.csv so we start fresh on each run
if os.path.exists(Analysis.STATS_PATH):
    os.remove(Analysis.STATS_PATH)

filter_results = {}
for T, p in itertools.product(settings.state.temperatures, settings.state.pressures):
    df_pool = dipole_results[
        (dipole_results['T'] == T) & (dipole_results['p'] == p)
    ].reset_index(drop=True)
    filter_results[(T, p)] = Analysis.plot_state_point(
        df_pool, T, p, mu_Vacuum=mu_Vacuum,    
        mu_gas_exp=settings.molecule.gas_dipole,
        k=settings.advanced.filtering.iqr_k,
        max_iqr_passes=settings.advanced.filtering.max_iqr_passes,
        min_pass_outliers=settings.advanced.filtering.min_pass_outliers,
        min_kept_configs=settings.advanced.filtering.min_kept_configs,
    )

keep_masks = {tp: r['keep_mask'] for tp, r in filter_results.items()}

Analysis.write_filtered_per_leaf(dipole_results, keep_masks)

per_replica_df, system_df = Analysis.build_collection(
    dipole_results=dipole_results,
    leaf_summary=leaf_summary,
    keep_masks=keep_masks,
    filter_results=filter_results,
    settings=settings,
    mu_Vacuum=mu_Vacuum,
    dipole_model=dipole_model,
    PCM1=PCM1,
    PCM2=PCM2,
    gas_stats=Gauss.gas_stats,
)

per_replica_df.to_csv('./Results/Per_Replica.csv', index=False)
system_df.to_csv('./Results/Results.csv', index=False)
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
        leaf = f'Simulations/replica_{replica}/{T:g}K/{p:g}Bar'
        with cd(leaf):
            MD.run_md3(settings, role='mixture', HOMEDIR=HOMEDIR, T=T, p=p)
