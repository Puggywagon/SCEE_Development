#!/usr/bin/python3
import gromacs
import subprocess
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import socket
import time

################################################################################
class Gro_Simulations(object):
    def __init__(self):
        pass
################################################################################
    def run_md1(self, settings, role):
        if role == 'pure':
            central_gro = 'Solvent.gro'
            topology = 'Pure.top'
            suffix = ''
        elif role == 'mixture':
            central_gro = 'Solute.gro'
            topology = 'Mixture.top'
            suffix = '_solute'
        else:
            raise ValueError(f"role must be 'pure' or 'mixture', got {role!r}")
    
            
        L = settings.advanced.configuration.box_length_nm + 1.7
        # The +1.7 nm is empirically determined: gmx insert-molecules
        # needs box headroom to place molecules without clipping. Values
        # 1.6 and 1.8 nm have been tested and only 1.7 worked reliably.
    
        print(f'performing vacuum equilibration of {central_gro}')
        gromacs.check(c=central_gro)
        gromacs.editconf(f=central_gro, box=[L, L, L], o=central_gro)
        gromacs.check(c=central_gro)

        # Energy minimisation
        gromacs.grompp(
            f='minim_vacuum.mdp', c=central_gro, p=topology,
            o=f'em{suffix}.tpr', maxwarn=2,
        )
        gromacs.mdrun(
            '-v', deffnm=f'em{suffix}',
            ntmpi=self.ntmpi, ntomp=self.ntomp,
            pin=self.pin, pinoffset=self.pinoffset,
        )
        
        # NVT vacuum equilibration
        gromacs.grompp(
            f='nvt_vacuum.mdp', c=f'em{suffix}.gro', p=topology,
            o=f'nvt_vacuum2{suffix}.tpr', maxwarn=2,
        )
        gromacs.mdrun(
            '-v', deffnm=f'nvt_vacuum2{suffix}',
            ntmpi=self.ntmpi, ntomp=self.ntomp,
            pin=self.pin, pinoffset=self.pinoffset,
        )
        
        self.insert_Molecules(settings, role)
        self.Write_to_top(settings, role)
################################################################################
    def insert_Molecules(self, settings, role):
        # Derive the suffix and central-gro from role
        if role == 'pure':
            suffix = ''
        elif role == 'mixture':
            suffix = '_solute'
        else:
            raise ValueError(f"role must be 'pure' or 'mixture', got {role!r}")
    
        # Base: the equilibrated central molecule (output of run_md1)
        base_gro = f'nvt_vacuum2{suffix}.gro'
        
        # Adding: the molecule to be replicated into the box
        if settings.molecule.aa_surround:
            # In mixture-AA, the bath is the (un-equilibrated) AA solvent
            adding_gro = 'Solvent.gro'
        else:
            # UA bath: user-supplied UA solvent gro
            adding_gro = 'Solvent_UA.gro'
        
        output_gro = f'out{suffix}.gro'
        
        print(f'inserting {settings.initial_molecules} copies of {adding_gro} into {base_gro}')
        gromacs.insert_molecules(
            f=base_gro,
            ci=adding_gro,
            nmol=str(settings.initial_molecules),
            o=output_gro,
        )
################################################################################
    def Write_to_top(self,settings, role):
        if role == 'pure':
            topology = 'Pure.top'
        elif role == 'mixture':
            topology = 'Mixture.top'
        else:
            raise ValueError(f"role must be 'pure' or 'mixture', got {role!r}")
    
        if settings.molecule.aa_surround:
            bath_name = 'Solvent'
        else:
            bath_name = 'Solvent_UA'
        with open(topology, 'a') as file:
            file.write(f'{bath_name}            {settings.initial_molecules}\n')
################################################################################
    def run_md2(self, settings, role):
        if role == 'pure':
            topology = 'Pure.top'
            suffix = ''
        elif role == 'mixture':
            topology = 'Mixture.top'
            suffix = '_solute'
        else:
            raise ValueError(f"role must be 'pure' or 'mixture', got {role!r}")
    
        # Energy minimisation of the full box
        gromacs.grompp(
            f='minim.mdp', c=f'out{suffix}.gro', p=topology,
            o=f'em1{suffix}.tpr', maxwarn=2,
        )
        gromacs.mdrun(
            '-v', deffnm=f'em1{suffix}',
            ntmpi=self.ntmpi, ntomp=self.ntomp,
            pin=self.pin, pinoffset=self.pinoffset,
        )
    
        # NVT equilibration
        gromacs.grompp(
            f='nvt.mdp', c=f'em1{suffix}.gro', p=topology,
            o=f'nvt_eq{suffix}.tpr', maxwarn=2,
        )
        gromacs.mdrun(
            '-v', deffnm=f'nvt_eq{suffix}',
            ntmpi=self.ntmpi, ntomp=self.ntomp,
            pin=self.pin, pinoffset=self.pinoffset,
        )
    
        # NPT equilibration
        gromacs.grompp(
            f='npt.mdp', c=f'nvt_eq{suffix}.gro', p=topology,
            o=f'npt{suffix}.tpr', maxwarn=2,
        )
        gromacs.mdrun(
            '-v', deffnm=f'npt{suffix}',
            ntmpi=self.ntmpi, ntomp=self.ntomp,
            pin=self.pin, pinoffset=self.pinoffset,
        )
################################################################################
    def create_mdpfile(self,HOMEDIR,mdpfile,T,p):

        T_kelvin = T  # temperature / K
        p_bar = p  # pressure / bar
        replace_dict = {'TEMPERATURE': f'{T_kelvin}',
                        'PRESSURE': f'{p_bar}'}
        search_text = 'TEMPERATURE', 'PRESSURE'
        replace_text = f'{T_kelvin}', f'{p_bar}'
        with open(HOMEDIR + '/template_2.mdp', 'r') as file: #Note that template 1 is the version with PR barostat while template 2 is the version with the c-resync barostat
            data = file.read()
            for search_text, replace_text in replace_dict.items():
                data = data.replace(search_text, replace_text)
        with open(mdpfile, 'w') as file:
            file.write(data)
################################################################################
    def process_trajectory(self,settings, role):
        print('processing trajectory')
        if role == 'pure':
            runname = 'Pure_QMMM_md3'
        elif role == 'mixture':
            runname = 'Mixture_QMMM_md3'
        else:
            raise ValueError(f"role must be 'pure' or 'mixture', got {role!r}")
    
        print('processing trajectory')
        tmp = gromacs.tools.Trjconv(
            f=f'{runname}.xtc',
            s=f'{runname}.tpr',
            o='conf_.gro',
            input='0\n',
            pbc='whole',
            b=1000, dt=20, sep=True,
        )
        chk, stdout, stderr = tmp.run()

################################################################################
    def extract_vacuum_configurations(self, settings, role='pure'):
        """Extract single-molecule configurations from the gas-phase NVT run
        (nvt_vacuum2) for the gas-phase QM calculations (Vacuum, PCM1, PCM2).

        Frames are taken every dt_ps after start_ps, and the first
        n_configurations of them are used. Returns a list of .gro paths.
        If n_configurations <= 1, the old behaviour is kept: only the final
        frame (nvt_vacuum2.gro) is used.
        """
        if role == 'pure':
            suffix = ''
        elif role == 'mixture':
            suffix = '_solute'
        else:
            raise ValueError(f"role must be 'pure' or 'mixture', got {role!r}")

        gas = settings.advanced.gas_sampling
        if gas.n_configurations <= 1:
            return [f'nvt_vacuum2{suffix}.gro']

        # Use the compressed trajectory if the vacuum run wrote one, otherwise the .trr
        traj = f'nvt_vacuum2{suffix}.xtc'
        if not os.path.exists(traj):
            traj = f'nvt_vacuum2{suffix}.trr'
        if not os.path.exists(traj):
            raise FileNotFoundError(
                f"No trajectory (nvt_vacuum2{suffix}.xtc or .trr) found for gas-phase "
                f"configuration extraction. Check the output settings in nvt_vacuum.mdp."
            )

        outdir = 'Gas_Configs'
        os.makedirs(outdir, exist_ok=True)
        # Clear old frames so a shorter re-extraction doesn't mix with older files
        for old in os.listdir(outdir):
            if old.startswith('gas_conf_') and old.endswith('.gro'):
                os.remove(os.path.join(outdir, old))

        print(f'extracting gas-phase configurations from {traj}')
        tmp = gromacs.tools.Trjconv(
            f=traj,
            s=f'nvt_vacuum2{suffix}.tpr',
            o=f'{outdir}/gas_conf_.gro',
            input='0\n',
            pbc='whole',
            b=gas.start_ps, dt=gas.dt_ps, sep=True,
        )
        chk, stdout, stderr = tmp.run()

        # trjconv -sep numbers frames from 0; sort numerically, keep the first N
        frames = [f for f in os.listdir(outdir)
                  if f.startswith('gas_conf_') and f.endswith('.gro')]
        frames.sort(key=lambda s: int(s[len('gas_conf_'):-len('.gro')]))
        if len(frames) < gas.n_configurations:
            print(f'WARNING: only {len(frames)} gas-phase frames available '
                  f'(requested {gas.n_configurations}). The vacuum run may be too '
                  f'short for start_ps={gas.start_ps}, dt_ps={gas.dt_ps}.')
        frames = frames[:gas.n_configurations]
        return [os.path.join(outdir, f) for f in frames]
################################################################################
    def run_md3(self, settings, role, HOMEDIR, T, p):
        gromacs.environment.flags['capture_output'] = False
        if role == 'pure':
            topology = HOMEDIR+f'/Pure.top'
            runname = f'Pure_QMMM_md3'
            suffix = ''
        elif role == 'mixture':
            topology = 'Mixture.top'
            suffix = '_solute'
            runname = HOMEDIR+f'/Mixture_QMMM_md3'
        else:
            raise ValueError(f"role must be 'pure' or 'mixture', got {role!r}")
    
        mdpfile='Prod.mdp'
        self.create_mdpfile(HOMEDIR,mdpfile,T,p)
        print('performing molecular dynamics')
        grofile = HOMEDIR+f'/npt{suffix}.gro'
        edrfile = f'{runname}.edr'
        groout = f'{runname}.gro'
        xtcfile = f'{runname}.xtc'
        trrfile = f'{runname}.trr'
        tprfile = f'{runname}.tpr'

        gromacs.grompp(f=mdpfile, c=grofile, p=topology, o=tprfile, maxwarn=2)
        gromacs.mdrun(
            '-v', s=tprfile, c=groout, o=trrfile, x=xtcfile, e=edrfile,
            ntmpi=self.ntmpi, ntomp=self.ntomp,
            pin=self.pin, pinoffset=self.pinoffset,
        )
    
        self.process_trajectory(settings, role)
################################################################################
