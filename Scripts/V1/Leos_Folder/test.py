#!/usr/bin/python3

import argparse
import yaml

import SCEE_GROMACS
import SCEE_Gaussian

import Gro_Builder
import Gro_Simulations
import Simulations_Analysis
import Oniom_Generation
import Gaussian_Calculations


################################################################################
################################################################################
################################################################################




with open("Settings_LL.yml", "r") as f:
    config = yaml.safe_load(f)


################################################################################
################################################################################
################################################################################
Gro_Builder=Gro_Builder.Gro_Builder()
MD=Gro_Simulations.Gro_Simulations()
Analysis=Simulations_Analysis.Simulations_Analysis() #Make this into an analysis script rather than what I have it now.
Oniom_Generation=Oniom_Generation.Oniom_Generation()
Gauss=Gaussian_Calculations.Gaussian_Calculations()

md = SCEE_GROMACS.SCEE_GROMACS()



print('Building *.gro file')
Build_Gro=config['User']["Build_Gro"]
Topology_File=Build_Gro["Solvent_Topology_File"]
solvent=Build_Gro["solvent"]
density=Build_Gro["density"]
mol_mass=Build_Gro["mol_mass"]
solresnametop=Build_Gro["solresnametop"]
solresname=Build_Gro["solresname"]
solmol=Build_Gro["solmol"]
Box_Build='Yes'
print('before')
Gro_file = md.AA_Structure(solmol,solvent,solresname)
Gro_File=Gro_Builder.AA_Structure(solmol,solvent,solresname) #solute #Rename this to handle UA and AA models rather than how I have been doing things..., get this to return the gro file name #Need to figure out how to maintain consistency between gro structure and top file.


################################################################################
################################################################################
qm = SCEE_Gaussian.SCEE_Gaussian()

itp = 'EthaneDiol.itp'
ref_itp = 'ffnonbonded.itp'
qm.create_atom_list(itp, ref_itp)
qm.write_atom_list()

#for atom in qm.atom_list:
#    print(atom['nr'], atom['type'], atom['ptype'], atom['charge'], atom['sigma'], atom['epsilon'])
    


qm['max_jobs'] = config['User']['Gaussian']['max_jobs']
qm['g09root'] = config['User']['Gaussian']['g09root']
qm['GAUSS_SCRDIR'] = config['User']['Gaussian']['GAUSS_SCRDIR']
qm['nproc'] = 8
qm['mem'] = '5GB'





################################################################################
################################################################################
################################################################################

