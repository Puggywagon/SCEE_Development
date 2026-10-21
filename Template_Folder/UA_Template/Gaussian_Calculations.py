#!/usr/bin/python3

import glob
import gromacs
import pandas as pd
import re
from subprocess import check_call
import os
import numpy as np


import Atoms
import Simulations_Analysis
from Shell_Oniom import ShellOniom

Atoms_Dict=Atoms.Atoms()
Analysis=Simulations_Analysis.Simulations_Analysis() 
################################################################################
class Gaussian_Calculations(object):
    def __init__(self, settings):
        self.settings = settings
        
        # Gaussian runtime config
        gauss = settings.software.gaussian
        self.nproc = gauss.nproc
        self.mem = gauss.mem
        self.max_jobs = gauss.max_jobs
        
        if gauss.g09root_override is None:
            raise NotImplementedError(
                "Auto-detection of g09root not yet implemented. "
                "Set Software.Gaussian.g09root_override in Settings.yml."
            )
        self.g09root = gauss.g09root_override
        self.scratch_dirs = [gauss.scratch_dir]
        
        # Electronic structure: v0 = optimisation step, v1 = single point
        es = settings.advanced.electronic_structure
        self.method_v0 = es.v0.method
        self.basis_v0 = es.v0.basis
        self.method_v1 = es.v1.method
        self.basis_v1 = es.v1.basis
#################################################################################      
    def gro_to_dat(self):
        f=open('nvt_vacuum2.gro')
        gro_file = f.readlines()        
        f.close()
        gro_file = gro_file[2:-1]
        gro_dict=self.get_gro(gro_file)  
        gro_df=pd.DataFrame(gro_dict)
        gro_df.columns=['Residue','gro_atom','atom_number','x','y','z','Vel x','Vel y','Vel z']
        dat_atoms=[]
        dat_xs=[]
        x_mean=gro_df['x'].mean()
        dat_ys=[]
        y_mean=gro_df['y'].mean()
        dat_zs=[]
        z_mean=gro_df['z'].mean()
        for index, row in gro_df.iterrows():
            Gros=gro_df.iloc[index]['gro_atom']
            dat_atom = Atoms_Dict.gaussian_symbol_from_gro_name(row['gro_atom'])
            dat_x=gro_df.iloc[index]['x']
            dat_x=(dat_x-x_mean)*10
            dat_y=gro_df.iloc[index]['y']
            dat_y=(dat_y-y_mean)*10
            dat_z=gro_df.iloc[index]['z']
            dat_z=(dat_z-z_mean)*10
            dat_atoms.append(dat_atom)
            dat_xs.append(dat_x)
            dat_ys.append(dat_y)
            dat_zs.append(dat_z)
        return dat_atoms,dat_xs,dat_ys,dat_zs
#################################################################################            
    def get_gro(self, gro_file):
        gro_dict = []
        for line in gro_file:
            data = line.split()
            gro_dicts = {
                'Residue': data[0],
                'gro_atoms': data[1],
                'atom_number': data[2],
                'x': float(data[3]),
                'y': float(data[4]),
                'z': float(data[5]),
            }
            # Velocities are optional in .gro format
            if len(data) >= 9:
                gro_dicts['Vel x'] = float(data[6])
                gro_dicts['Vel y'] = float(data[7])
                gro_dicts['Vel z'] = float(data[8])
            gro_dict.append(gro_dicts)
        return gro_dict
################################################################################
    def read_multipoles(self, filename):
    
        f = open(filename, 'r')
        text_file = f.read()
        f.close()
        
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
    def get_multipole_statistics_scee(self, rundir):
        raw_dict = {}
        case_dict = {1: 'dipole_l', 2: 'dipole_m', 3: 'dipole_h'}
    
        for Q in range(1, 4):
            filelist = glob.glob(rundir + f'*_c*_q{Q}_v1.out')
            for filename in filelist:
                multipole = self.read_multipoles(filename)
                config = filename.split('_')[-3].replace('c', '')
                if config not in raw_dict:
                    raw_dict[config] = {}
                raw_dict[config][case_dict[Q]] = multipole['total dipole']
        
        data_dict = {'config': [], 'dipole_l': [], 'dipole_m': [], 'dipole_h': []}
        for config, value in raw_dict.items():
            dipole_l = value.get('dipole_l')
            dipole_m = value.get('dipole_m')
            dipole_h = value.get('dipole_h')
            if not any(np.isnan(x) for x in [dipole_l, dipole_m, dipole_h]):
                data_dict['config'].append(config)
        case_dict = {1: 'dipole_l', 2: 'dipole_m', 3: 'dipole_h'}

        for Q in range(1, 4):
            filelist = glob.glob(rundir + f'*_c*_q{Q}_v1.out')
            for filename in filelist:
                multipole = self.read_multipoles(filename)
                config = filename.split('_')[-3].replace('c', '')
                if config not in raw_dict:
                    raw_dict[config] = {}
                raw_dict[config][case_dict[Q]] = multipole['total dipole']

        data_dict = {'config': [], 'dipole_l': [], 'dipole_m': [], 'dipole_h': []}
        for config, value in raw_dict.items():
            dipole_l = value.get('dipole_l')
            dipole_m = value.get('dipole_m')
            dipole_h = value.get('dipole_h')
            if not any(np.isnan(x) for x in [dipole_l, dipole_m, dipole_h]):
                data_dict['config'].append(config)
                data_dict['dipole_l'].append(dipole_l)
                data_dict['dipole_m'].append(dipole_m)
                data_dict['dipole_h'].append(dipole_h)

        return pd.DataFrame.from_dict(data_dict)
################################################################################
    def get_single_dipole(self, rundir, name):
        filename = rundir + f'{name}.out'
        multipole = self.read_multipoles(filename)
        return multipole['total dipole']  
################################################################################
    def process_gro(self, oniom_inp_path):
        shell = ShellOniom()
        shell.process(oniom_inp_path)

################################################################################
    def run_gaussian(self, rundir, run_str, step=0):

        step_str = ''
        if (step == 1):
            step_str = '_v1'
        
        scratch_list_bash = ' '.join(self.scratch_dirs)

        text  = '#!/bin/bash' + '\n'
        text += f'cd {rundir}' + '\n'
        text += f'export g09root={self.g09root}' + '\n'
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
    def init0(self):
        rundir = 'Vacuum/'
        sol_keyword = self.settings.molecule.sol_keyword
        dat_atoms,dat_xs,dat_ys,dat_zs=self.gro_to_dat()
    
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
                        
        text  += f'--Link1--\n'              
        text += f'%chk=Vacuum.chk' + '\n'        
        text += f'%nprocshared={self.nproc}' + '\n'
        text += f'#p {self.method_v0}/{self.basis_v0} gfprint fopt=tight' + '\n'
        text += '# nosymm pop=full scf=(verytight) density=current\n'
        text += '# integral=(ultrafine,NoXCTest) geom=checkpoint\n\n'
        text += f'{sol_keyword}\n\n'
        text += '0 1 \n\n\n\n'
        
        
        if (not os.path.isdir(rundir)):
            os.mkdir(rundir)
            
        f = open(rundir+'Vacuum.dat', 'w')
        f.write(text)
        f.close()
        run_str='Vacuum'
        self.run_gaussian(rundir,run_str,step=0)
        
        mu_Vacuum = self.get_single_dipole(rundir, 'Vacuum')
        pd.DataFrame({'Vacuum': [mu_Vacuum]}).to_csv('Dipole_Vacuum.csv', index=False)
        return mu_Vacuum
################################################################################
    def init1(self):
        sol_keyword = self.settings.molecule.sol_keyword
        exp_diconst = self.settings.molecule.dielectric_constant
    
        rundir = 'PCM1/'
        dat_atoms,dat_xs,dat_ys,dat_zs=self.gro_to_dat()
                
        text = f'%chk=PCM1.chk\n'
        text += f'%nprocshared={self.nproc}\n'
        text += f'%mem={self.mem}\n'
        text += f'#p {self.method_v0}/{self.basis_v0} gfprint scrf=(pcm,solvent={sol_keyword},read)' + '\n'
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
        text += f'#p {self.method_v1}/{self.basis_v1} gfprint scrf=(pcm,solvent={sol_keyword},read)' + '\n'
        text += '# nosymm pop=full scf=(verytight) density=current' + '\n'
        text += '# integral=(ultrafine,NoXCTest) geom=checkpoint\n\n'
        text += f'{sol_keyword}\n\n'
        text += '0 1 \n\n'
        text += f'eps={exp_diconst}\n\n\n\n'
        
        if (not os.path.isdir(rundir)):
            os.mkdir(rundir)
        
        f = open(rundir+'PCM1.dat', 'w')
        f.write(text)
        f.close()
                
        run_str='PCM1'
        self.run_gaussian(rundir,run_str,step=0)
        PCM1 = self.get_single_dipole(rundir, 'PCM1')
        pd.DataFrame({'PCM1': [PCM1]}).to_csv('Dipole_PCM1.csv', index=False)
        return PCM1
################################################################################
    def init2(self):
        sol_keyword = self.settings.molecule.sol_keyword
        cal_diconst = self.settings.molecule.calculated_dielectric
    
        rundir = 'PCM2/'
        dat_atoms,dat_xs,dat_ys,dat_zs=self.gro_to_dat()

        
        text = f'%chk=PCM2.chk' + '\n'        
        text += f'%nprocshared={self.nproc}\n'
        text += f'%mem={self.mem}\n'
        text += f'#p {self.method_v0}/{self.basis_v0} gfprint scrf=(pcm,solvent={sol_keyword},read)' + '\n'
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
        text += f'#p {self.method_v1}/{self.basis_v1} gfprint scrf=(pcm,solvent={sol_keyword},read)' + '\n'
        text += '# nosymm pop=full scf=(verytight) density=current' + '\n'
        text += '# integral=(ultrafine,NoXCTest) geom=checkpoint\n\n'
        text += f'{sol_keyword}\n\n'
        text += '0 1 \n\n'
        text += f'eps={cal_diconst}\n\n\n\n' 
               
        if (not os.path.isdir(rundir)):
            os.mkdir(rundir)
        
        f = open(rundir+'PCM2.dat', 'w')
        f.write(text)
        f.close()
                
        run_str='PCM2'
        self.run_gaussian(rundir,run_str,step=0)

        PCM2 = self.get_single_dipole(rundir, 'PCM2')
        pd.DataFrame({'PCM2': [PCM2]}).to_csv('Dipole_PCM2.csv', index=False)
        
        return PCM2
################################################################################
    def init3(self):
    
        rundir = 'SCEE/'
        
        text  = f'%nprocshared={self.nproc}\n'
        text += f'%mem={self.mem}\n'
        text += f'#p oniom({self.method_v0}/{self.basis_v0}:amber=(print,SoftFirst,LastEquiv))=embedcharge gfprint' + '\n'
        text += '# opt=quadmacro iop(1/98=66,1/19=7)\n'
        text += '# nosymm pop=full scf=(verytight) density=current\n'
        text += '# integral=(ultrafine,NoXCTest) geom=connectivity\n\n'
        text += 'qteste\n\n'
        text += '0 1 0 1 0 1\n'

        if (not os.path.isdir(rundir)):
            os.mkdir(rundir)
        filelist = glob.glob('*_c*_q?.inp')
        print(filelist)
        for inpfile in filelist:
            text_str = text
            
            file0_str = inpfile.replace('.inp', '')
            print(file0_str)
            file0_str = file0_str.split('/')[-1]
            datfile = rundir + file0_str + '.dat'
            print(datfile)
        
            f = open(inpfile, 'r')
            text_str += f.read()
            f.close()

            f = open(datfile, 'w')
            f.write(text_str)
            f.close()
                
        run_str='*_c*_q?'
        self.run_gaussian(rundir,run_str,step=0)
        
################################################################################
    def init4(self):
        rundir = 'SCEE/'

        text  = f'%nprocshared={self.nproc}\n'
        text += f'%mem={self.mem}\n'
        text += f'#p {self.method_v1}/{self.basis_v1} gfprint charge' + '\n'
        text += '# nosymm pop=full scf=(verytight) density=current' + '\n'
        text += '# integral=(ultrafine,NoXCTest)\n\n'
        text += 'teste\n\n'
        text += '0 1\n'
        
        re_str = 'Z-Matrix orientation:(?:.*\n){' + str(5+self.natom) + '}'
        
        filelist = glob.glob('*_c*_q?_chg.inp')
        print(filelist)
        
        
        Model_Dipole, _ =Analysis.get_dipole_model_liquid() 
        ratio=Model_Dipole / self.qmax
        num1=self.settings.advanced.charge_scaling.qr1*ratio*self.qmax
        num2=self.settings.advanced.charge_scaling.qr2*ratio*self.qmax
        num3=self.settings.advanced.charge_scaling.qr3*ratio*self.qmax
        for inpfile in filelist:

            file0_str = inpfile.replace('_chg.inp', '')
            file0_str = file0_str.split('/')[-1]
            
            text_str = f'%chk={file0_str}.chk' + '\n' + text
            
            outfile = rundir + file0_str + '.out'
            datfile = rundir + file0_str + '_v1.dat'
            print(outfile, datfile)
            f = open(outfile, 'r')
            text_file = f.read()
            raw_data = re.findall(re_str, text_file, re.M)
            data = raw_data[-1].split('\n')
            for line in data[5:5+self.natom]:
                row = line.split()
                atom = Atoms_Dict.Atom_Symbol(int(row[1]))
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
                
        run_str='*_c*_q?'
        self.run_gaussian(rundir,run_str,step=1)
        scee_df = self.get_multipole_statistics_scee(rundir)
        print(scee_df.describe())
        print(scee_df.head())
        muL_list = []
        for index, row in scee_df.iterrows():
             x = np.array([num1, num2, num3])
             y = np.array([row['dipole_l'], row['dipole_m'], row['dipole_h']])
             coeff = np.polyfit(x, y, 2)
             b = coeff[1] - 1
             fact = b*b-4*coeff[0]*coeff[2]
             if (coeff[0] >= 1.0e-6) and (fact > 0):
                  muL = (-b - np.sqrt(fact)) / (2*coeff[0])
             else:
                 coeff = np.polyfit(x, y, 1)
                 b = coeff[0] - 1
                 muL = - coeff[1]/b
             muL_list.append(muL)
        scee_df['muL_SCEE'] = muL_list
        scee_df.to_csv('Dipole_SCEE.csv', index=False)
        
        return scee_df
################################################################################
