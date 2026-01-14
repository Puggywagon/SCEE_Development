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
import Oniom_Generation

Gro_Builder=Gro_Builder.Gro_Builder()
Oniom_Generation=Oniom_Generation.Oniom_Generation()

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
# Here we calculate the model dipole moment of the production run box and use this in calculating the scaled charges
Model_Dipole=Pre_Eq_Simulations.get_dipole_model_liquid(system_title) # We could use this to extract the liquid dielectric constant
ratio=Model_Dipole / qmax
num1=qr1*ratio*qmax
num2=qr2*ratio*qmax
num3=qr3*ratio*qmax

