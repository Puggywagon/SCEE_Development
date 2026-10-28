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

        # Filled by run_gas_phase: {'Vacuum': {...}, 'PCM1': {...}, 'PCM2': {...}}
        self.gas_stats = {}
#################################################################################      
    def gro_to_dat(self, gro_path='nvt_vacuum2.gro'):
        """Read a single-molecule .gro and return element symbols and coordinates
        (Angstrom, centred on the geometric centre).

        Fixed-width columns are used (as in Shell_Oniom.parse_gro), so frames
        with or without velocities both work.
        """
        with open(gro_path) as f:
            lines = f.readlines()
        n_atoms = int(lines[1].strip())

        names, xs, ys, zs = [], [], [], []
        for line in lines[2:2 + n_atoms]:
            names.append(line[10:15].strip())
            xs.append(float(line[20:28]))
            ys.append(float(line[28:36]))
            zs.append(float(line[36:44]))

        x_mean, y_mean, z_mean = np.mean(xs), np.mean(ys), np.mean(zs)
        dat_atoms = [Atoms_Dict.gaussian_symbol_from_gro_name(n) for n in names]
        dat_xs = [(x - x_mean) * 10 for x in xs]
        dat_ys = [(y - y_mean) * 10 for y in ys]
        dat_zs = [(z - z_mean) * 10 for z in zs]
        return dat_atoms, dat_xs, dat_ys, dat_zs
################################################################################
    def read_multipoles(self, filename):
    
        f = open(filename, 'r')
        text_file = f.read()
        f.close()
        
        multipole = {}
            
        re_dipole = r' Dipole moment \(field-independent basis, Debye\):(?:.*\n){2}'
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
        
        re_quad = r' Quadrupole moment \(field-independent basis, Debye-Ang\):(?:.*\n){3}'
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

        re_oct = r' Octapole moment \(field-independent basis, Debye-Ang\*\*2\):(?:.*\n){4}'
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
            # Skip configs with a missing (None) or failed (NaN) charge-set result
            if not any(x is None or np.isnan(x) for x in [dipole_l, dipole_m, dipole_h]):
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
    def run_gaussian(self, rundir, run_str, step=0, datfiles=None):
        # datfiles: optional explicit list of .dat names (inside rundir) to run
        # instead of the run_str glob

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
        if datfiles is not None:
            text += f'for datfile in {" ".join(datfiles)}; do' + '\n'
        else:
            text += f'for datfile in $(ls {run_str}{step_str}.dat); do' + '\n'
        text += '  outfile=${datfile%.dat}.out' + '\n'
        text += '  idx=$(( i % ${#SCRATCH_DIRS[@]} ))' + '\n'
        text += '  GAUSS_SCRDIR=${SCRATCH_DIRS[$idx]}' + '\n'
        text += '  echo "Preparing $datfile -> $outfile on $GAUSS_SCRDIR"' + '\n'

        # --- GLOBAL THROTTLE: never exceed MAXJOBS g09 processes on this machine ---
        text += '  while true; do' + '\n'
        text += '      g09c=$(pgrep -c g09 2>/dev/null)\n'
        text += '      mdc=$(pgrep -fc "gmx mdrun" 2>/dev/null)\n'
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
    def _gas_phase_input(self, kind, chk_name, dat_atoms, dat_xs, dat_ys, dat_zs):
        """Build the two-step Gaussian input used for Vacuum, PCM1 and PCM2.

        Step 1 (v0): Vacuum -> geometry optimisation (fopt=tight);
                     PCM1/PCM2 -> PCM calculation at the given geometry
                     (unchanged from the previous init1/init2).
        Step 2 (v1): single point at the step-1 geometry (geom=checkpoint).
                     The dipole moment is read from this step, so all three
                     are now at the same level as the SCEE v1 calculations.
        """
        sol_keyword = self.settings.molecule.sol_keyword
        if kind == 'Vacuum':
            scrf, eps_line, opt = '', '', ' fopt=tight'
        elif kind == 'PCM1':
            eps = self.settings.molecule.dielectric_constant
            scrf, eps_line, opt = f' scrf=(pcm,solvent={sol_keyword},read)', f'eps={eps}\n\n', ''
        elif kind == 'PCM2':
            eps = self.settings.molecule.calculated_dielectric
            scrf, eps_line, opt = f' scrf=(pcm,solvent={sol_keyword},read)', f'eps={eps}\n\n', ''
        else:
            raise ValueError(f"kind must be 'Vacuum', 'PCM1' or 'PCM2', got {kind!r}")

        # --- Step 1: v0 ---
        text  = f'%chk={chk_name}.chk\n'
        text += f'%nprocshared={self.nproc}\n'
        text += f'%mem={self.mem}\n'
        text += f'#p {self.method_v0}/{self.basis_v0} gfprint{opt}{scrf}\n'
        text += '# nosymm scf=(verytight) density=current\n'
        text += '# integral=(ultrafine,NoXCTest)\n\n'
        step1 = 'optimisation' if kind == 'Vacuum' else 'v0_calc'
        text += f'{kind}_{step1} {chk_name}\n\n'
        text += '0 1\n'
        for dat_atom, dat_x, dat_y, dat_z in zip(dat_atoms, dat_xs, dat_ys, dat_zs):
            text += f'{dat_atom:s} {dat_x:9.4f} {dat_y:8.4f} {dat_z:8.4f}\n'
        text += '\n'
        text += eps_line

        # --- Step 2: v1 single point on the step-1 geometry ---
        text += '--Link1--\n'
        text += f'%chk={chk_name}.chk\n'
        text += f'%nprocshared={self.nproc}\n'
        text += f'%mem={self.mem}\n'
        text += f'#p {self.method_v1}/{self.basis_v1} gfprint{scrf}\n'
        text += '# nosymm pop=full scf=(verytight) density=current\n'
        text += '# integral=(ultrafine,NoXCTest) geom=checkpoint\n\n'
        text += f'{kind}_dipole_calc {chk_name}\n\n'
        text += '0 1\n\n'
        text += eps_line
        text += '\n\n'
        return text
################################################################################
    def _gas_job_complete(self, outfile):
        """True if both steps (v0 and v1) of a gas-phase job terminated normally."""
        if not os.path.exists(outfile):
            return False
        with open(outfile, errors='ignore') as f:
            return f.read().count('Normal termination') >= 2
################################################################################
    def run_gas_phase(self, kind, gro_files, run=True):
        """Run Vacuum, PCM1 or PCM2 on every gas-phase configuration and
        average the v1 dipole moments.

        Writes Dipole_<kind>_configs.csv (one row per configuration) and
        Dipole_<kind>.csv (mean, SD, SE, N). Returns (mean, SE).
        Set run=False to only re-read existing .out files.
        """
        rundir = f'{kind}/'
        os.makedirs(rundir, exist_ok=True)

        names = []
        pending = []
        for i, gro in enumerate(gro_files):
            name = f'{kind}_c{i}'
            names.append(name)
            # Leave finished jobs alone, so a killed/restarted run only redoes
            # the unfinished ones (Gaussian overwrites a .out as soon as it starts)
            if self._gas_job_complete(rundir + f'{name}.out'):
                continue
            dat_atoms, dat_xs, dat_ys, dat_zs = self.gro_to_dat(gro)
            with open(rundir + f'{name}.dat', 'w') as f:
                f.write(self._gas_phase_input(kind, name, dat_atoms, dat_xs, dat_ys, dat_zs))
            pending.append(f'{name}.dat')

        if run and pending:
            print(f'{kind}: {len(names) - len(pending)} already complete, running {len(pending)}')
            self.run_gaussian(rundir, f'{kind}_c*', step=0, datfiles=pending)
        elif run:
            print(f'{kind}: all {len(names)} already complete, nothing to run')

        rows = []
        for name, gro in zip(names, gro_files):
            outfile = rundir + f'{name}.out'
            mu = self.get_single_dipole(rundir, name) if os.path.exists(outfile) else np.nan
            rows.append({'config': name, 'source_gro': gro, 'dipole': mu})
        df = pd.DataFrame(rows)
        df.to_csv(f'Dipole_{kind}_configs.csv', index=False)

        vals = df['dipole'].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        n_failed = len(df) - vals.size
        if n_failed:
            print(f'WARNING: {n_failed} of {len(df)} {kind} calculations returned no dipole')
        if vals.size == 0:
            raise RuntimeError(f'No {kind} dipole moments could be read from {rundir}')

        mean = float(vals.mean())
        sd = float(vals.std(ddof=1)) if vals.size > 1 else np.nan
        se = sd / np.sqrt(vals.size) if vals.size > 1 else np.nan
        pd.DataFrame({kind: [mean], f'SD_{kind}': [sd], f'SE_{kind}': [se],
                      f'N_{kind}': [vals.size]}).to_csv(f'Dipole_{kind}.csv', index=False)
        self.gas_stats[kind] = {'mean': mean, 'SD': sd, 'SE': se, 'N': int(vals.size)}
        print(f'{kind}: mean = {mean:.4f} D, SD = {sd:.4f}, SE = {se:.4f}, N = {vals.size}')
        return mean, se
################################################################################
    def init0(self, gro_files=('nvt_vacuum2.gro',), run=True):
        mu_Vacuum, _ = self.run_gas_phase('Vacuum', list(gro_files), run=run)
        return mu_Vacuum
################################################################################
    def init1(self, gro_files=('nvt_vacuum2.gro',), run=True):
        PCM1, _ = self.run_gas_phase('PCM1', list(gro_files), run=run)
        return PCM1
################################################################################
    def init2(self, gro_files=('nvt_vacuum2.gro',), run=True):
        PCM2, _ = self.run_gas_phase('PCM2', list(gro_files), run=run)
        return PCM2
################################################################################
    def init3(self, run=True):
    
        rundir = 'SCEE/'
        
        text  = f'%nprocshared={self.nproc}\n'
        text += f'%mem={self.mem}\n'
        text += f'#p oniom({self.method_v0}/{self.basis_v0}:amber=(print,SoftFirst,LastEquiv))=embedcharge gfprint' + '\n'
        text += '# opt=quadmacro iop(1/98=66,1/19=7)\n'
        text += '# nosymm scf=(verytight) density=current\n'
        text += '# integral=(ultrafine,NoXCTest) geom=connectivity\n\n'
        text += 'SCEE_optimisation\n\n'
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
        if run:
            self.run_gaussian(rundir,run_str,step=0)
        
################################################################################
    def init4(self, run=False):
        rundir = 'SCEE/'

        text  = f'%nprocshared={self.nproc}\n'
        text += f'%mem={self.mem}\n'
        text += f'#p {self.method_v1}/{self.basis_v1} gfprint charge' + '\n'
        text += '# nosymm pop=full scf=(verytight) density=current' + '\n'
        text += '# integral=(ultrafine,NoXCTest)\n\n'
        text += 'SCEE_dipole_calc\n\n'
        text += '0 1\n'
        
        re_str = 'Z-Matrix orientation:(?:.*\n){' + str(5+self.natom) + '}'
        
        filelist = glob.glob('*_c*_q?_chg.inp')
        print(f'{len(filelist)} SCEE charge files found')
        
        
        Model_Dipole, _ =Analysis.get_dipole_model_liquid() 
        ratio=Model_Dipole / self.qmax
        num1=self.settings.advanced.charge_scaling.qr1*ratio*self.qmax
        num2=self.settings.advanced.charge_scaling.qr2*ratio*self.qmax
        num3=self.settings.advanced.charge_scaling.qr3*ratio*self.qmax
        # The v1 inputs are only needed if the v1 jobs are going to be run.
        # When just re-reading results, skip rebuilding them (this means
        # reading all 600 large ONIOM .out files, which is slow).
        for inpfile in (filelist if run else []):

            file0_str = inpfile.replace('_chg.inp', '')
            file0_str = file0_str.split('/')[-1]
            
            text_str = f'%chk={file0_str}.chk' + '\n' + text
            
            outfile = rundir + file0_str + '.out'
            datfile = rundir + file0_str + '_v1.dat'
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
        if run:
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
