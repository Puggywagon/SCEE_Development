#!/usr/bin/python3

import io
import itertools
import os
import pandas as pd
import subprocess

import MDAnalysis

import Atom


################################################################################
################################################################################
################################################################################
class SCEE_Gaussian(object):

    def __init__(self, config):

        self.param = {}
        self.atom_list = []
        
        self.nproc = 8
        self.param['nproc'] = 8
        self.param['rundir'] = './opt_b3lypaug-cc-pvtz_spb3lyp/'
        self.param['g09root'] ='/home/mjksill/Research/PhD/MacPherson/SCEE/Gaussian/g09_pgi/'
        self.param['GAUSS_SCRDIR'] = '/home/mjksill/idrive/ProcessModelling/mjksill/scratch/'

        self.param['basis_vacuum'] = '6-31+G(d,p)'
        self.param['method_vacuum'] = 'b3lyp'
        self.param['basis_v1'] = 'aug-cc-pvtz'
        self.param['method_v1'] = 'blyp'


        for key, value in config.items():
            self.param[key] = value
#            print(key, value)
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
    def run_gaussian(self, run_str, step=0):

        step_str = ''
        if (step == 1):
            step_str = '_v1'
        
        scratch_list_bash = ' '.join(self.scratch_dirs)

        text  = '#!/bin/bash' + '\n'
        text += f'cd {self.rundir}' + '\n'
        text += f'export g09root={self.param["g09root"]}' + '\n'
        text += 'source $g09root/g09/bsd/g09.profile' + '\n'

        # scratch directories as a bash array (can be one or many)
        scratch_list = ' '.join(self.scratch_dirs)
        text += f'SCRATCH_DIRS=({scratch_list})' + '\n'

        # global maximum number of Gaussian jobs allowed on this machine
        text += f'MAXJOBS={self.max_jobs}' + '\n'
        text += 'i=0' + '\n'

        # loop over all .dat files for this step (v0 or v1)
        text += f'for datfile in $(ls {run_str}{step_str}.dat); do' + '\n'
        text += '  outfile=${datfile%.dat}.out' + '\n'
        text += '  idx=$(( i % ${#SCRATCH_DIRS[@]} ))' + '\n'
        text += '  GAUSS_SCRDIR=${SCRATCH_DIRS[$idx]}' + '\n'
        text += '  echo "Preparing $datfile -> $outfile on $GAUSS_SCRDIR"' + '\n'

        # --- GLOBAL THROTTLE: never exceed MAXJOBS g09 processes on this machine ---
        text += '  while true; do' + '\n'
        text += '      g09c=$(pgrep -c g09 2>/dev/null)\n'
        text += '      mdc=$(pgrep -c mdrun 2>/dev/null)\n'
        text += '      g09c=${g09c:-0}\n'
        text += '      mdc=${mdc:-0}\n'
        text += '      if [ $((g09c + mdc)) -lt "$MAXJOBS" ]; then\n'
        text += '          break' + '\n'
        text += '      fi' + '\n'
        text += '      sleep 5' + '\n'
        text += '  done' + '\n'

        # launch Gaussian job in the background, with its own GAUSS_SCRDIR
        text += '  GAUSS_SCRDIR=$GAUSS_SCRDIR g09 < "${datfile}" > "${outfile}" &' + '\n'
        text += '  i=$((i+1))' + '\n'
        text += 'done' + '\n'
        text += 'wait' + '\n'

        f = open(f'run_{step_str}.sh', 'w')
        f.write(text)
        f.close()

        # make sure the script is executable (chmod 755)
        os.chmod(f'run_{step_str}.sh', 0o755)

        cmd = ['bash', f'run_{step_str}.sh']
        check_call(cmd)

            
################################################################################
################################################################################
    def init_vacuum(self, workdir='./'):
        sol_keyword = 'Ethanol'
        rundir = 'SCEE_Gaussian_vacuum/'
        text = f'%chk=Vacuum.chk\n'
        text += f'%nprocshared={self.param["nproc"]}\n'
        text += f'%mem={self.param["mem"]}\n'
        text += f'#p {self.param["method_vacuum"]}/{self.param["basis_vacuum"]} gfprint fopt=tight' + '\n'
        text += '# nosymm pop=full scf=(verytight) density=current\n'
        text += '# integral=(ultrafine,NoXCTest)\n\n'
        text += f'{sol_keyword}\n\n'
        text += '0 1\n'

        gro_data = MDAnalysis.Universe('./Ethanol_UA/nvt_vacuum2.gro')
        for atom in gro_data.atoms:
            x = atom.position[0]
            y = atom.position[1]
            z = atom.position[2]
            text += f'{atom.resname:s} {x:9.4f} {y:8.4f} {z:8.4f}\n'
        text += '\n'  

#        if (not os.path.isdir(workdir + rundir)):
#            os.mkdir(workdir + rundir)
#            
#        f = open(workdir + rundir+'Vacuum.dat', 'w')
#        f.write(text)
#        f.close()
#        run_str='vacuum'
#        self.run_gaussian(run_str,step=0)
#        dips_out=run_str

#        df4 = self.get_multipole_statistics(dips_out,rundir)
#        muG_list=df4.mean()
#        df4['muG'] = muG_list
#        df4.to_csv('Dipole_Vacuum.csv', index=False)
#        print(muG_list)
#        
#        mu_Vacuum=df4['muG']
#        return mu_Vacuum

        
        print(text)
        
        
################################################################################
################################################################################
################################################################################

