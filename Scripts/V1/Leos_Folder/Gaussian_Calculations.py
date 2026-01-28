#!/usr/bin/python3

import glob
import gromacs
import pandas as pd
import re
from subprocess import check_call
import os
from dat_generation import Dat_Generation
import numpy as np


import Atoms

Atoms_Dict=Atoms.Atoms()

################################################################################
class Gaussian_Calculations(object):
    def __init__(self):
        self.runname = 'TIP4P2005'       

        #self.basis_v0 = 'aug-cc-pvtz'
        #self.method_v0 = 'b3lyp'
        #self.basis_v1 = 'aug-cc-pvtz'
        #self.method_v1 = 'b3lyp'
        self.basis_v0 = '6-31+G(d,p)'
        self.method_v0 = 'b3lyp'
        self.basis_v1 = 'aug-cc-pvtz'
        self.method_v1 = 'blyp'
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
    def get_multipole_statistics(self,rundir):
        raw_dict = {}
        case_dict = {1: 'dipole_l', 2: 'dipole_m', 3: 'dipole_h'}        
        for Q in range(1,4):
#            print(Q)
            filelist = glob.glob(rundir + f'*_c*_q{Q}_v1.out')
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
    def process_gro(self, exe, inp):
    
        shutil.copy2(HOMEDIR +'/' + 'oniom.inp', '.')
        shutil.copy2(HOMEDIR +'/' + 'shellO', '.')
        shutil.copy2(HOMEDIR +'/' + 'Shell_Oniom.f90', '.')
        print('processing conf_*.gro files')
        logfile = open('junk.log', 'w')
        cmd = [exe]
        check_call(cmd, stdout=logfile, stderr=logfile)

################################################################################
    def run_gaussian(self, step=0):

        step_str = ''
        if (step == 1):
            step_str = '_v1'
        
        scratch_list_bash = ' '.join(self.scratch_dirs)

        text  = '#!/bin/bash' + '\n'
        text += f'cd {self.rundir}' + '\n'
        text += f'export g09root={self.g09root}' + '\n'
        text += 'source $g09root/g09/bsd/g09.profile' + '\n'

        # scratch directories as a bash array (can be one or many)
        scratch_list = ' '.join(self.scratch_dirs)
        text += f'SCRATCH_DIRS=({scratch_list})' + '\n'

        # global maximum number of Gaussian jobs allowed on this machine
        text += f'MAXJOBS={self.max_jobs}' + '\n'
        text += 'i=0' + '\n'

        # loop over all .dat files for this step (v0 or v1)
        text += f'for datfile in $(ls *_c*_q?{step_str}.dat); do' + '\n'
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
    def init(self,sol_keyword,Gaus='Vacuum', workdir='./'):
        rundir = 'Vacuum/'
        dat_atoms,dat_xs,dat_ys,dat_zs=self.dat_generation.gro_to_dat()
        
        text = f'%chk=Vacuum.chk\n'
        text += f'%nprocshared={self.nproc}\n'
        text += f'%mem={self.mem}\n'
        text += f'#p {self.method_v0}/{self.basis_v0} gfprint fopt=tight' + '\n'
        text += '# nosymm pop=full scf=(verytight) density=current\n'
        text += '# integral=(ultrafine,NoXCTest)\n\n'
        text += f'{sol_keyword}\n\n'
        text += '0 1\n'
        for dat_atom,dat_x,dat_y,dat_z in zip(dat_atoms,dat_xs,dat_ys,dat_zs):
            text += f'{dat_atom:s} {dat_x:9.4f} {dat_y:8.4f} {dat_z:8.4f}\n'
        text += '\n'  
        
        if (not os.path.isdir(workdir + rundir)):
            os.mkdir(workdir + rundir)
            
        f = open(workdir + rundir+'Vacuum.dat', 'w')
        f.write(text)
        f.close()
        self.run_gaussian(step=0)
    
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
        print(muG_list)
        
        mu_Vacuum=df4['muG']
        return mu_Vacuum
################################################################################
    def init(self,sol_keyword,exp_diconst,Gaus='PCM1', workdir='./'):
        rundir = 'PCM1/'
        dat_atoms,dat_xs,dat_ys,dat_zs=self.dat_generation.gro_to_dat()
                
        text = f'%chk=PCM1.chk\n'
        text += f'%nprocshared={self.nproc}\n'
        text += f'%mem={self.mem}\n'
        text += f'#p b3lyp/cc-pvtz gfprint scrf=(pcm,solvent={sol_keyword},read)' + '\n'
        text += '# nosymm pop=full density=current scf=(verytight)\n'
        text += '# integral=(ultrafine,NoXCTest)\n\n'
        text += f'{sol_keyword}\n\n'
        text += '0 1 \n'
        
        for dat_atom,dat_x,dat_y,dat_z in zip(dat_atoms,dat_xs,dat_ys,dat_zs):
            text += f'{dat_atom:s} {dat_x:9.4f} {dat_y:8.4f} {dat_z:8.4f}\n'
        text +='\n'
        text += f'eps={exp_diconst}\n\n'

        text  += f'--Link1--\n'
        text += f'%chk=PCM1.chk\n'
        text += f'%nprocshared={self.nproc}\n'
        text += f'%mem={self.mem}\n'
        text += f'#p b3lyp/aug-cc-pvtz gfprint scrf=(pcm,solvent={sol_keyword},read)' + '\n'
        text += '# nosymm pop=full scf=(verytight) density=current' + '\n'
        text += '# integral=(ultrafine,NoXCTest) geom=checkpoint\n\n'
        text += f'{sol_keyword}\n\n'
        text += '0 1 \n\n'
        text += f'eps={exp_diconst}\n\n\n\n'
        
        if (not os.path.isdir(workdir + rundir)):
            os.mkdir(workdir + rundir)
        
        f = open(workdir + rundir+'PCM1.dat', 'w')
        f.write(text)
        f.close()
        
        
        self.run_gaussian(step=0)
          
        df2 = PCM1.get_multipole_statistics(rundir)
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
        return PCM1
################################################################################
    def init(self,sol_keyword,cal_diconst,Gaus='PCM2', workdir='./'):
        rundir = 'PCM2/'
        dat_atoms,dat_xs,dat_ys,dat_zs=self.dat_generation.gro_to_dat()

        
        text = f'%chk=PCM2.chk' + '\n'        
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
        text += f'%chk=PCM2.chk' + '\n'   
        text += f'%nprocshared={self.nproc}\n'
        text += f'%mem={self.mem}\n'
        text += f'#p b3lyp/aug-cc-pvtz gfprint scrf=(pcm,solvent={sol_keyword},read)' + '\n'
        text += '# nosymm pop=full scf=(verytight) density=current' + '\n'
        text += '# integral=(ultrafine,NoXCTest) geom=checkpoint\n\n'
        text += f'{sol_keyword}\n\n'
        text += '0 1 \n\n'
        text += f'eps={cal_diconst}\n\n\n\n' 
               
        if (not os.path.isdir(workdir + rundir)):
            os.mkdir(workdir + rundir)
        
        f = open(workdir + rundir+'PCM2.dat', 'w')
        f.write(text)
        f.close()
        
        PCM2.run_gaussian(step=0)

        df3 = PCM2.get_multipole_statistics(rundir)
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
        PCM2=df3['muL_PCM2']
        
        return PCM2
################################################################################
    def init(self,Gaus='SCEE_V0', workdir='./'):
    
        rundir = 'SCEE/'
        
        text  = f'%nprocshared={self.nproc}\n'
        text += f'%mem={self.mem}\n'
        text += f'#p oniom({self.method_v0}/{self.basis_v0}:amber=(print,SoftFirst,LastEquiv))=embedcharge gfprint' + '\n'
        text += '# opt=quadmacro iop(1/98=66,1/19=7)\n'
        text += '# nosymm pop=full scf=(verytight) density=current\n'
        text += '# integral=(ultrafine,NoXCTest) geom=connectivity\n\n'
        text += 'qteste\n\n'
        text += '0 1 0 1 0 1\n'

        if (not os.path.isdir(workdir + rundir)):
            os.mkdir(workdir + rundir)
        filelist = glob.glob(workdir + '*_c*_q?.inp')
        print(filelist)
        for inpfile in filelist:
            text_str = text
            
            file0_str = inpfile.replace('.inp', '')
            print(file0_str)
            file0_str = file0_str.split('/')[-1]
            datfile = workdir + rundir + file0_str + '.dat'
            print(datfile)
        
            f = open(inpfile, 'r')
            text_str += f.read()
            f.close()
#            print(text_str)

            f = open(datfile, 'w')
            f.write(text_str)
            f.close()
        self.run_gaussian(step=0)
        
################################################################################
    def init(self,Gaus='SCEE_V1', workdir='./'):
        rundir = 'SCEE/'

        text  = f'%nprocshared={self.nproc}\n'
        text += f'%mem={self.mem}\n'
        text += f'#p {self.method_v1}/{self.basis_v1} gfprint charge' + '\n'
        text += '# nosymm pop=full scf=(verytight) density=current' + '\n'
        text += '# integral=(ultrafine,NoXCTest)\n\n'
        text += 'teste\n\n'
        text += '0 1\n'
        
        re_str = 'Z-Matrix orientation:(?:.*\n){' + str(5+self.natom) + '}'
        
        filelist = glob.glob(workdir + '*_c*_q?_chg.inp')
        print(filelist)
        for inpfile in filelist:

            file0_str = inpfile.replace('_chg.inp', '')
            file0_str = file0_str.split('/')[-1]
            
            text_str = f'%chk={file0_str}.chk' + '\n' + text
            
            outfile = workdir + rundir + file0_str + '.out'
            datfile = workdir + rundir + file0_str + '_v1.dat'
            print(outfile, datfile)
            f = open(outfile, 'r')
            text_file = f.read()
            raw_data = re.findall(re_str, text_file, re.M)
            data = raw_data[-1].split('\n')
            for line in data[5:5+self.natom]:
                row = line.split()
                atom = Atoms.Atom_Symbol(Atom=int(row[1]))
                x, y, z = row[3], row[4], row[5]
                text_str += f'{atom} {x} {y} {z}' + '\n'
            f.close()
        
            text_str += '\n'
            f = open(inpfile, 'r')
            text_str += f.read()
            f.close()

            f = open(datfile, 'w')
            f.write(text_str)
            f.close()
            
        self.run_gaussian(step=1)
        df1 = Molecule.get_multipole_statistics(rundir)
        print(df1.describe())
        print(df1.head())
        muL_list = []
        for index, row in df1.iterrows():
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
             muL_list.append(muL)
        df1['muL_SCEE'] = muL_list
        df1.to_csv('Dipole_SCEE.csv', index=False)
        SCEE=df1['muL_SCEE']
        return SCEE
################################################################################
