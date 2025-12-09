#!/usr/bin/python3
import subprocess
import numpy as np
import shutil
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd
from statistics import mean
import math 
import csv
import gromacs
import argparse

import Molecule
import Particle
import SimBox
#########################################################################################################
class Analysis(object):
    def __init__(self):
        pass
#########################################################################################################
    def get_density(self,system_title):
        gromacs.environment.flags['capture_output'] = True
        input_str = 'Density 0\n'
        tmp = gromacs.tools.Energy(f=f'{system_title}_QMMM_md3',
                               s=f'{system_title}_QMMM_md3',
                               input=input_str)
        chk, density_output, stderr = tmp.run()

        f = open(f'Density.txt', 'w')
        f.write(density_output)
        f.close()

        density_lines = density_output.splitlines()
        density_line = density_lines[9]  # Assuming density is on the 5th line
        density = float(density_line.split()[1])
        
        return density
###############################################################################################
    def Averagre_dipole(self):
        df = pd.read_csv(f'./Dipole.csv')
        df.columns = ['config', 'dipole_l', 'dipole_m', 'dipole_h','muL']
        dipole_moments=df['muL']
        dipole=np.mean(dipole_moments)
        std_dipole=np.std(dipole_moments)
        
        return dipole_moments,dipole,std_dipole
###############################################################################################
    def make_ndx(self,system_title):
        input_1='del 2 \n del 1 \n del 0 \n a O* \n name 0 Oxygen \n a H* \n name 1 Hydrogen \n q'
        tmp = gromacs.tools.Make_ndx(f=f'{system_title}_QMMM_md3',
                                o='index.ndx',
                                input=input_1)
        chk, makendx_output, stderr = tmp.run()
################################################################################
    def rdf_OO(self,system_title):
        input_str='0\n0\n'
        tmp = gromacs.tools.Rdf(f=f'{system_title}_QMMM_md3', 
                                s=f'{system_title}_QMMM_md3',
                                n='index.ndx',
                                excl=[],
                                bin=0.004,
                                o='OO_RDF.xvg',
                                input=input_str)
        chk, rdf_output, stderr = tmp.run()
################################################################################
    def rdf_OH(self,system_title):
        input_str='0\n 1\n'
        tmp = gromacs.tools.Rdf(f=f'{system_title}_QMMM_md3', 
                                s=f'{system_title}_QMMM_md3',
                                n='index.ndx',
                                bin=0.004,
                                excl=[],
                                o='OH_RDF.xvg',
                                input=input_str)
        chk, rdf_output, stderr = tmp.run()
################################################################################
    def rdf_HH(self,system_title):
        input_str='1\n 1\n'
        tmp = gromacs.tools.Rdf(f=f'{system_title}_QMMM_md3', 
                                s=f'{system_title}_QMMM_md3',
                                n='index.ndx',
                                bin=0.004,
                                excl=[],
                                o='HH_RDF.xvg',
                                input=input_str)
        chk, rdf_output, stderr = tmp.run()
##############################################################################################
    def read_conf(self,filename):
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()
    
        lines.pop(0)
        nmol = lines.pop(0)
        natoms = 4
#        print(nmol)
#        print(lines[-1])
        tmp = lines.pop(-1)
        data = tmp.split()
        x = float(data[0])
        y = float(data[1])
        z = float(data[2])
        box = SimBox.SimBox(x, y, z)
#        print(tmp)

        molecule_dict = {}
        for line in lines:
            data = line.split()
            m_name = data[0]
            if (m_name not in molecule_dict):
                m = Molecule.Molecule(m_name)
                m.box = box
                molecule_dict[m_name] = m
                
            p_name = data[1]
            p = Particle.Particle(p_name)
            x = float(data[3])
            y = float(data[4])
            z = float(data[5])
            p.position = np.array([x, y, z])
        
            molecule_dict[m_name].add_particle(p)

        for name, m in molecule_dict.items():
            m.add_acceptor(m.particle_list[0])
            m.add_donor(m.particle_list[1], m.particle_list[0])
            m.add_donor(m.particle_list[2], m.particle_list[0])

        return molecule_dict
################################################################################
    def hbonding_anlaysis(self):
        df = pd.read_csv('Dipole.csv')
        df.sort_values(by='config', inplace=True)
        print(df.head())

        file_list = glob.glob('*_*.gro')

        #filename = 'conf_1.gro'
        data_dict = {'config': [], 'donors': [], 'acceptors': [], 'nHB': []}
        for filename in file_list:
            tmp = filename.split('_')[-1]
            config = tmp.replace('.gro', '')
            molecule_dict = self.read_conf(filename)
    
            molecule_list = list(molecule_dict.values())
            m0 = molecule_list.pop(0)
    
            donors = 0
            acceptors = 0
            for m in molecule_list:
                if (m0.is_acceptor(m)):
                    acceptors += 1
                if (m0.is_donor(m)):
                    #print(m.ID)
                    donors += 1
    
            nHB = acceptors + donors        
            print(f'{config}:  nHB = {nHB}, donors={donors}, acceptors={acceptors}')
            data_dict['config'].append(config)
            data_dict['donors'].append(donors)
            data_dict['acceptors'].append(acceptors)
            data_dict['nHB'].append(nHB)
        
        tmp = pd.DataFrame.from_dict(data_dict)
        tmp.sort_values(by='config', inplace=True)

        df['donors'] = tmp['donors']
        df['acceptors'] = tmp['acceptors']
        df['nHB'] = tmp['nHB']
        
        No_Hbonds= df['nHB']
        No_Donors=df['donors']
        No_Acceptors=df['acceptors']
    
        df.to_csv('hbanalysis.csv', index=False)
        
        return No_Hbonds, No_Donors, No_Acceptors
################################################################################
    def build_CSV(self,density,dipole_moments,dipole,std_dipole,No_Hbonds,No_Donors,No_Acceptors):
        data={
        'Replicas':i,
        'Temperature':T,
        'Pressure':P,
        'Density':Density,
        'Dipole Distribution':dipole_moments,
        'Average Dipole Moments':dipole,
        'STD Of Dipole Moments':std_dipole,
        'Number of H_bonds':No_Hbonds,
        'Number of H-Donors':No_Donors,
        'Number of H-Acceptors':No_Acceptors
        }
        df=pd.DataFrame(data)
        with open(f"../../../../Results.csv", "a", newline="") as file:  # Open in append mode
            df.to_csv(file, header=False, index=False)
            return df
################################################################################

