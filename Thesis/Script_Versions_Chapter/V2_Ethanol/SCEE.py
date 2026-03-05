#!/usr/bin/python3

import glob
import pandas as pd
import re
from subprocess import check_call
import os
import sys
from shutil import which
import numpy as np

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


################################################################################
################################################################################
################################################################################
class SCEE(object):

    def __init__(self):
        self.runname = 'TIP4P2005'
        self.atom_dict = {1: 'H', 6: 'C', 7:'N', 8: 'O'}
        self.nproc = 8
        self.basis_v0 = 'aug-cc-pvtz'
        self.method_v0 = 'b3lyp'
        self.basis_v1 = 'aug-cc-pvtz'
        self.method_v1 = 'b3lyp'
        #self.basis_v0 = '6-31+G(d,p)'
        #self.method_v0 = 'b3lyp'
        #self.basis_v1 = 'aug-cc-pvtz'
        #self.method_v1 = 'blyp'
        self.natom = 3

        self.rundir = './opt_b3lypaug-cc-pvtz_spb3lyp/'
        self.g09root='/home/zoe/Software/Gaussian/g09_pgi'
        self.GAUSS_SCRDIR = '/home/zoe/Research/Gaussian/scratch'
        
        

################################################################################
################################################################################
    def process_gro(self, exe, inp):
        print('processing conf_*.gro files')
        logfile = open('junk.log', 'w')
        cmd = [exe]
        check_call(cmd, stdout=logfile, stderr=logfile)
        

################################################################################
################################################################################
    def init_v0(self, workdir='./'):
        
        text = f'%nproc={self.nproc}' + '\n'
        text +='%mem=50MW'+'\n'
        text += f'#p oniom({self.method_v0}/{self.basis_v0}:amber=(print,SoftFirst,LastEquiv))=embedcharge gfprint' + '\n'
        text += '# opt=quadmacro iop(1/98=66,1/19=7)\n'
        text += '# nosymm pop=full scf=(verytight) density=current\n'
        text += '# integral=(ultrafine,NoXCTest) geom=connectivity\n\n'
        text += 'qteste\n\n'
        text += '0 1 0 1 0 1\n'

        if (not os.path.isdir(workdir + self.rundir)):
            os.mkdir(workdir + self.rundir)
        filelist = glob.glob(workdir + '*_c*_q?.inp')
        print(filelist)
        for inpfile in filelist:
            text_str = text
            
            file0_str = inpfile.replace('.inp', '')
            print(file0_str)
            file0_str = file0_str.split('/')[-1]
            datfile = workdir + self.rundir + file0_str + '.dat'
            print(datfile)
        
            f = open(inpfile, 'r')
            text_str += f.read()
            f.close()
#            print(text_str)

            f = open(datfile, 'w')
            f.write(text_str)
            f.close()
        

################################################################################
    def init_v1(self, workdir='./'):

        text  = f'%nproc={self.nproc}' + '\n'
        text +='%mem=50MW'+'\n'
        text += f'#p {self.method_v1}/{self.basis_v1} gfprint charge' + '\n'
        text += '# nosymm pop=full scf=(verytight) density=current' + '\n'
        text += '# integral=(ultrafine,NoXCTest)\n\n'
        text += 'teste\n\n'
        text += '0 1\n'
        
        re_str = 'Z-Matrix orientation:(?:.*\n){' + str(5+self.natom) + '}'
        
        filelist = glob.glob(workdir + '*_c*_q?_chg.inp')
        print(filelist)
        for inpfile in filelist:
            if not os.path.exists(inpfile) or os.path.getsize(inpfile) == 0:
                print(f"Skipping {inpfile} — missing or empty _chg.inp (v0 must have failed)")
                continue

            file0_str = inpfile.replace('_chg.inp', '')
            file0_str = file0_str.split('/')[-1]
    
            text_str = f'%chk={file0_str}.chk' + '\n' + text
    
            outfile = workdir + self.rundir + file0_str + '.out'
            datfile = workdir + self.rundir + file0_str + '_v1.dat'
            print(outfile, datfile)

            if not os.path.exists(outfile) or os.path.getsize(outfile) == 0:
                print(f"Skipping {outfile} — missing or empty Gaussian output")
                continue

            with open(outfile, 'r') as f:
                text_file = f.read()

            if "Z-Matrix orientation" not in text_file:
                print(f"Skipping {outfile} — no Z-Matrix found (Gaussian likely failed)")
                continue

            raw_data = re.findall(re_str, text_file, re.M)
            if not raw_data:
                print(f"Skipping {outfile} — regex match failed on Z-Matrix block")
                continue

            data = raw_data[-1].split('\n')
            for line in data[5:5+self.natom]:
                row = line.split()
                atom = self.atom_dict[int(row[1])]
                x, y, z = row[3], row[4], row[5]
                text_str += f'{atom} {x} {y} {z}' + '\n'

            text_str += '\n'
            with open(inpfile, 'r') as f:
                text_str += f.read()

            with open(datfile, 'w') as f:
                f.write(text_str)


################################################################################
    def run_gaussian(self, step=0):

        step_str = ''
        if (step == 1):
            step_str = '_v1'
        
        text  = '#!/bin/bash' + '\n'
        text += f'cd {self.rundir}' + '\n'
        text += f'export g09root={self.g09root}' + '\n'
        text += 'source $g09root/g09/bsd/g09.profile' + '\n'
        text += f'export GAUSS_SCRDIR={self.GAUSS_SCRDIR}' + '\n'
        text += f'for datfile in $(ls *_c*_q?{step_str}.dat); do' + '\n'    
        text += '  outfile=${datfile%.dat}.out' + '\n'    
        text += '  echo "g09 < ${datfile} > ${outfile}"' + '\n'    
        text += '  g09 < ${datfile} > ${outfile}' + '\n'    
        text += 'done' + '\n'

        f = open(f'run_{step_str}.sh', 'w')
        f.write(text)
        f.close()

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
        
        if raw_data:
            data = raw_data[-1].split('\n')
            dipole_data = data[1].split()
            multipole['x'] = float(dipole_data[1])
            multipole['y'] = float(dipole_data[3])
            multipole['z'] = float(dipole_data[5])
            multipole['total dipole'] = float(dipole_data[7])        
        else:
            multipole['x'] = np.nan
            multipole['y'] = np.nan
            multipole['z'] = np.nan
            multipole['total dipole'] = np.nan
        
        re_quad = ' Quadrupole moment \(field-independent basis, Debye-Ang\):(?:.*\n){3}'
        raw_data = re.findall(re_quad, text_file, re.M)
        if raw_data:        
            data = raw_data[-1].split('\n')
            quad_data = data[1].split()
            multipole['xx'] = float(quad_data[1])
            multipole['yy'] = float(quad_data[3])
            multipole['zz'] = float(quad_data[5])
            quad_data = data[2].split()
            multipole['xy'] = float(quad_data[1])
            multipole['xz'] = float(quad_data[3])
            multipole['yz'] = float(quad_data[5])
        else:
            multipole['xx'] = np.nan
            multipole['yy'] = np.nan
            multipole['zz'] = np.nan
            multipole['xy'] = np.nan
            multipole['xz'] = np.nan
            multipole['yz'] = np.nan

        re_oct = ' Octapole moment \(field-independent basis, Debye-Ang\*\*2\):(?:.*\n){4}'
        raw_data = re.findall(re_oct, text_file, re.M)
        if raw_data:        
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
        else:
            multipole['xxx'] = np.nan
            multipole['yyy'] = np.nan
            multipole['zzz'] = np.nan
            multipole['xyy'] = np.nan
            multipole['xxy'] = np.nan
            multipole['xxz'] = np.nan
            multipole['xzz'] = np.nan
            multipole['yzz'] = np.nan
            multipole['yyz'] = np.nan
            multipole['xyz'] = np.nan
        
        return multipole
    
            
################################################################################
################################################################################
    def get_multipole_statistics(self):
        raw_dict = {}
        case_dict = {1: 'dipole_l', 2: 'dipole_m', 3: 'dipole_h'}        
        for Q in range(1,4):
#            print(Q)
            filelist = glob.glob(self.rundir + f'*_c*_q{Q}_v1.out')
            for filename in filelist:
                multipole = self.read_multipoles(filename)
                config = filename.split('_')[-3].replace('c', '')
                if config not in raw_dict:
                    raw_dict[config] = {}
                raw_dict[config][case_dict[Q]] = multipole['total dipole']
    
                #print(filename, config, case_dict[Q], multipole['total dipole'])
                
        data_dict = {'config': [], 'dipole_l': [], 'dipole_m': [], 'dipole_h': []}    
        for config, value in raw_dict.items():
            dipole_l=value.get('dipole_l')
            dipole_m=value.get('dipole_m')
            dipole_h=value.get('dipole_h')
            if not any(np.isnan(x) for x in [dipole_l,dipole_m,dipole_h]):
                data_dict['config'].append(config)
                data_dict['dipole_l'].append(dipole_l)
                data_dict['dipole_m'].append(dipole_m)
                data_dict['dipole_h'].append(dipole_h)
                #print(config)

        df = pd.DataFrame.from_dict(data_dict)
        return df           
################################################################################
################################################################################
################################################################################
        
