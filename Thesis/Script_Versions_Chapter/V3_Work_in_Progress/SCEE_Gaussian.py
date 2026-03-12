#!/usr/bin/python3

import io
import itertools
import os
import pandas as pd
import subprocess

import Atom


################################################################################
################################################################################
################################################################################
class SCEE_Gaussian(object):

    def __init__(self):

        self.param = {}
        self.atom_list = []
        
        self.nproc = 8
        self.basis = 'aug-cc-pvtz'
        self.method = 'b3lyp'

        self.rundir = './opt_b3lypaug-cc-pvtz_spb3lyp/'
        self.g09root='/home/mjksill/Research/PhD/MacPherson/SCEE/Gaussian/g09_pgi/'
        self.GAUSS_SCRDIR = '/home/mjksill/idrive/ProcessModelling/mjksill/scratch/'


################################################################################
################################################################################
    def __setitem__(self, key, value):
        self.param[key] = value

        
    def __getitem__(self, key):
        return self.param[key]


################################################################################
################################################################################
    def create_atom_list(self, itp, ref_itp):

        header = ['type', 'name', 'at.num', 'mass', 'charge', 'ptype', 'sigma', 'epsilon']
        df_opls = pd.read_csv(ref_itp, skiprows=3, sep=r'\s+', comment=';',
                              header=None, names=header,
                              index_col='type')
#        df_opls['at.num'] = df_opls['at.num'].astype(int)
#        print(df_opls['at.num'])

        with open(itp, 'r') as f:
            header_to_end = itertools.dropwhile(lambda x: x.strip() != '[ atoms ]', f)
            atoms_section = itertools.takewhile(lambda x: x.strip() != '', header_to_end)
            data_str = ''.join(atoms_section)
        header = ['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass']
        df = pd.read_csv(io.StringIO(data_str), skiprows=[0],
                         sep=r'\s+', header=None, names=header, comment=';')
#        print(df)

        self.atom_list = []
        for index, row in df.iterrows():
            atom = Atom.Atom()
            type_str = row['type']
            if (type_str in df_opls.index):
                row_opls = df_opls.loc[type_str]
                for key, value in row_opls.items():
                    atom[key] = value
            for key, value in row.items():
#                print(key, value)
                atom[key] = value
            self.atom_list.append(atom)


################################################################################
################################################################################
    def write_atom_list(self):

        text  = f"{'nr':>5s}"
        text += f"{'name':>8s}"
        text += f"{'type':>15s}"
        text += f"{'at.num':>8s}"
        text += f"{'ptype':>8s}"
        text += f"{'mass':>10s}"
        text += f"{'charge':>8s}"
        text += f"{'sigma':>10s}"
        text += f"{'epsilon':>10s}"
        print(text)
        for atom in self.atom_list:
            text  = f"{atom['nr']:5d}"
            text += f"{atom['name']:>8s}"
            text += f"{atom['type']:>15s}"
            text += f"{atom['at.num']:>8.0f}"
            text += f"{atom['ptype']:>8s}"
            text += f"{atom['mass']:>10.5f}"
            text += f"{atom['charge']:>8.3f}"
            text += f"{atom['sigma']:>10.5f}"
            text += f"{atom['epsilon']:>10.5f}"
            print(text)



################################################################################
################################################################################
################################################################################

