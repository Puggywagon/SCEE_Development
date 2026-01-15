#!/usr/bin/python3

import glob
import gromacs
import pandas as pd
import re
from subprocess import check_call
import os
from dat_generation import Dat_Generation

import Atoms

Atoms_Dict=Atoms.Atoms()
################################################################################
################################################################################
################################################################################
class PCM2(object):

    def __init__(self):
        self.runname = 'TIP4P2005'

        self.nproc = 4
        self.mem = '3GB'
        self.max_jobs = 3

        #self.basis_v0 = 'aug-cc-pvtz'
        #self.method_v0 = 'b3lyp'
        #self.basis_v1 = 'aug-cc-pvtz'
        #self.method_v1 = 'b3lyp'
        self.basis_v0 = '6-31+G(d,p)'
        self.method_v0 = 'b3lyp'
        self.basis_v1 = 'aug-cc-pvtz'
        self.method_v1 = 'blyp'
        self.natom = 6

        self.rundir = './opt_b3lypaug-cc-pvtz_spb3lyp/'
        self.g09root='/home/zoe/Software/Gaussian/g09_pgi'
        self.GAUSS_SCRDIR = '/home/zoe/Research/Gaussian/scratch' 
        self.dat_generation=Dat_Generation(filename='nvt_vacuum2.gro')
################################################################################
################################################################################
    def init(self,system_title,sol_keyword,cal_diconst, workdir='./'):
        dat_atoms,dat_xs,dat_ys,dat_zs=self.dat_generation.gro_to_dat()

        
        text = f'%chk={system_title}_pcm2.chk' + '\n'        
        text += f'%nprocshared={self.nproc}\n'
        text += f'%mem={self.mem}\n'
        text += f'#p b3lyp/cc-pvtz gfprint scrf=(pcm,solvent={sol_keyword},read)' + '\n'
        text += '# nosymm pop=full scf=(verytight) density=current\n'
        text += '# integral=(ultrafine,NoXCTest)\n\n'
        text += f'{sol_keyword}\n\n'
        text += '0 1 \n'
        for dat_atom,dat_x,dat_y,dat_z in zip(dat_atoms,dat_xs,dat_ys,dat_zs):
            text += f'{dat_atom:s} {dat_x:9.4f} {dat_y:8.4f} {dat_z:8.4f}\n'
        text += '\n'
        text += f'eps={cal_diconst} \n\n'
        
        text  += f'\n--Link1--\n'                    
        text += f'%chk={system_title}_pcm2.chk' + '\n'   
        text += f'%nprocshared={self.nproc}\n'
        text += f'%mem={self.mem}\n'
        text += f'#p b3lyp/aug-cc-pvtz gfprint scrf=(pcm,solvent={sol_keyword},read)' + '\n'
        text += '# nosymm pop=full scf=(verytight) density=current' + '\n'
        text += '# integral=(ultrafine,NoXCTest) geom=checkpoint\n\n'
        text += f'{sol_keyword}\n\n'
        text += '0 1 \n\n'
        text += f'eps={cal_diconst}\n\n\n\n' 
               
        if (not os.path.isdir(workdir + self.rundir)):
            os.mkdir(workdir + self.rundir)
        
        f = open(workdir + self.rundir+'PCM2.dat', 'w')
        f.write(text)
        f.close()
################################################################################
    def run_gaussian(self, step=0):

        step_str = ''

        text  = '#!/bin/bash\n'
        text += f'cd {self.rundir}\n'
        text += f'export g09root={self.g09root}\n'
        text += 'source $g09root/g09/bsd/g09.profile\n'
        text += f'export GAUSS_SCRDIR={self.GAUSS_SCRDIR}\n'
        text += f'MAXJOBS={self.max_jobs}\n'

        text += 'for datfile in $(ls *.dat); do\n'
        text += '  outfile=${datfile%.dat}.out\n'
        text += '  echo "Preparing $datfile -> $outfile"\n'
        text += '  while true; do\n'
        text += '      g09c=$(pgrep -c g09 2>/dev/null)\n'
        text += '      mdc=$(pgrep -c mdrun 2>/dev/null)\n'
        text += '      g09c=${g09c:-0}\n'
        text += '      mdc=${mdc:-0}\n'
        text += '      if [ $((g09c + mdc)) -lt "$MAXJOBS" ]; then\n'
        text += '          break\n'
        text += '      fi\n'
        text += '      sleep 5\n'
        text += '  done\n'
        text += '  g09 < "${datfile}" > "${outfile}" &\n'
        text += 'done\n'
        text += 'wait\n'

        f = open(f'run_{step_str}.sh', 'w')
        f.write(text)
        f.close()

        os.chmod(f'run_{step_str}.sh', 0o755)
        cmd = ['bash', f'run_{step_str}.sh']
        check_call(cmd)
        

################################################################################
    def read_multipoles(self, filename):
#        outfile = '../TIP4P/tmp/opt_b3lypaug-cc-pvtz_spb3lyp/TIP4P2005_c1_q1_v1.out'
#        outfile = '../tmp/opt-test/wat_SPCE_oniom_c1_q1_v1.out'
#        filename = './opt-orig/TIP4P2005_c14_q1_v1.out'
        f = open(filename, 'r')
        text_file = f.read()
        f.close()
    
#        re_str = ' Dipole moment \(field-independent basis, Debye\):(?:.*\n){17}'
#        raw_data = re.findall(re_str, text_file, re.M)
#        data = raw_data[-1].split('\n')
#        for line in data:
#            print(line)

        multipole = {}
            
        re_dipole = ' Dipole moment \(field-independent basis, Debye\):(?:.*\n){2}'
        raw_data = re.findall(re_dipole, text_file, re.M)
        data = raw_data[-1].split('\n')
#        print(data)
        dipole_data = data[1].split()
        multipole['x'] = float(dipole_data[1])
        multipole['y'] = float(dipole_data[3])
        multipole['z'] = float(dipole_data[5])
        multipole['total dipole'] = float(dipole_data[7])
#        print(multipole['total dipole'])

        re_quad = ' Quadrupole moment \(field-independent basis, Debye-Ang\):(?:.*\n){3}'
        raw_data = re.findall(re_quad, text_file, re.M)
        data = raw_data[-1].split('\n')
        quad_data = data[1].split()
        multipole['xx'] = float(quad_data[1])
        multipole['yy'] = float(quad_data[3])
        multipole['zz'] = float(quad_data[5])
        quad_data = data[2].split()
        multipole['xy'] = float(quad_data[1])
        multipole['xz'] = float(quad_data[3])
        multipole['yz'] = float(quad_data[5])

        re_oct = ' Octapole moment \(field-independent basis, Debye-Ang\*\*2\):(?:.*\n){4}'
        raw_data = re.findall(re_oct, text_file, re.M)
        data = raw_data[-1].split('\n')
        quad_data = data[1].split()
        multipole['xxx'] = float(quad_data[1])
        multipole['yyy'] = float(quad_data[3])
        multipole['zzz'] = float(quad_data[5])
        multipole['xyy'] = float(quad_data[7])
        quad_data = data[2].split()
        multipole['xxy'] = float(quad_data[1])
        multipole['xxz'] = float(quad_data[3])
        multipole['xzz'] = float(quad_data[5])
        multipole['yzz'] = float(quad_data[7])
        quad_data = data[3].split()
        multipole['yyz'] = float(quad_data[1])
        multipole['xyz'] = float(quad_data[3])
        
        return multipole
    
            
################################################################################
################################################################################
    def get_multipole_statistics(self):

        filename = self.rundir + f'PCM2.out'
        multipole = self.read_multipoles(filename)
        dipole=multipole['total dipole']
    
        print(multipole['total dipole'])
        data_dict={'dipole_l': [dipole],'dipole_m': [dipole],'dipole_h': [dipole]} 
        df3 = pd.DataFrame(data_dict)
        print(df3)
        return df3

            
################################################################################
################################################################################
################################################################################
        
