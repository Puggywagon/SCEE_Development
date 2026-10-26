#!/usr/bin/python3
import gromacs
import subprocess
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import socket
import time
import pyedr

################################################################################
class Gro_Simulations(object):
    def __init__(self):
        pass
################################################################################
    def run_md1(self, settings, role):
        files = settings.files
        if role == 'pure':
            central_gro = files.solvent_gro
            topology = files.pure_topology
            suffix = ''
        elif role == 'mixture':
            central_gro = files.solute_gro
            topology = files.mixture_topology
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
            f=settings.files.mdp_vacuum_em, c=central_gro, p=topology,
            o=f'em{suffix}.tpr', maxwarn=2,
        )
        gromacs.mdrun(
            '-v', deffnm=f'em{suffix}',
            ntmpi=self.ntmpi, ntomp=self.ntomp,
            pin=self.pin, pinoffset=self.pinoffset,
        )
        
        # NVT vacuum equilibration
        gromacs.grompp(
            f=settings.files.mdp_vacuum_nvt, c=f'em{suffix}.gro', p=topology,
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
            adding_gro = settings.files.solvent_gro
        else:
            # UA bath: user-supplied UA solvent gro
            adding_gro = settings.files.solvent_ua_gro
        
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
            topology = settings.files.pure_topology
        elif role == 'mixture':
            topology = settings.files.mixture_topology
        else:
            raise ValueError(f"role must be 'pure' or 'mixture', got {role!r}")
    
        bath_itp = (settings.files.solvent_itp if settings.molecule.aa_surround
                    else settings.files.solvent_ua_itp)
        bath_name = self._moleculetype_name(bath_itp)

        # Make sure we start on a new line, even if the topology file
        # was saved without a trailing newline
        with open(topology, 'rb') as file:
            file.seek(0, os.SEEK_END)
            needs_newline = file.tell() > 0 and (file.seek(-1, os.SEEK_END) or file.read(1) != b'\n')
        with open(topology, 'a') as file:
            if needs_newline:
                file.write('\n')
            file.write(f'{bath_name}            {settings.initial_molecules}\n')
################################################################################
    def bath_molecule_count(self, settings, role='pure'):
        """Number of bath molecules actually listed in [ molecules ] of the
        topology (the last entry), rather than assuming all requested
        molecules were inserted."""
        topology = settings.files.pure_topology if role == 'pure' else settings.files.mixture_topology
        bath_itp = (settings.files.solvent_itp if settings.molecule.aa_surround
                    else settings.files.solvent_ua_itp)
        bath_name = self._moleculetype_name(bath_itp)
        in_mol, count = False, None
        with open(topology) as f:
            for line in f:
                text = line.split(';')[0].strip()
                if not text:
                    continue
                if text.startswith('['):
                    in_mol = text.replace(' ', '').lower() == '[molecules]'
                    continue
                if in_mol:
                    parts = text.split()
                    if len(parts) >= 2 and parts[0] == bath_name:
                        count = int(parts[1])     # last matching entry = the bath
        return count
################################################################################
    def _moleculetype_name(self, itp_path):
        """Return the molecule name from the first [ moleculetype ] section of
        an .itp file, so the [ molecules ] entry matches whatever the force-field
        generator called it (e.g. 'phenol' rather than 'Solvent')."""
        in_section = False
        with open(itp_path) as f:
            for line in f:
                stripped = line.split(';')[0].strip()
                if not stripped:
                    continue
                if stripped.startswith('['):
                    in_section = stripped.replace(' ', '').lower() == '[moleculetype]'
                    continue
                if in_section:
                    return stripped.split()[0]
        raise ValueError(f'No [ moleculetype ] name found in {itp_path}')
################################################################################
    def run_md2(self, settings, role):
        if role == 'pure':
            topology = settings.files.pure_topology
            suffix = ''
        elif role == 'mixture':
            topology = settings.files.mixture_topology
            suffix = '_solute'
        else:
            raise ValueError(f"role must be 'pure' or 'mixture', got {role!r}")
    
        # Energy minimisation of the full box
        gromacs.grompp(
            f=settings.files.mdp_box_em, c=f'out{suffix}.gro', p=topology,
            o=f'em1{suffix}.tpr', maxwarn=2,
        )
        gromacs.mdrun(
            '-v', deffnm=f'em1{suffix}',
            ntmpi=self.ntmpi, ntomp=self.ntomp,
            pin=self.pin, pinoffset=self.pinoffset,
        )
    
        # NVT equilibration
        gromacs.grompp(
            f=settings.files.mdp_box_nvt, c=f'em1{suffix}.gro', p=topology,
            o=f'nvt_eq{suffix}.tpr', maxwarn=2,
        )
        gromacs.mdrun(
            '-v', deffnm=f'nvt_eq{suffix}',
            ntmpi=self.ntmpi, ntomp=self.ntomp,
            pin=self.pin, pinoffset=self.pinoffset,
        )
    
        # NPT equilibration
        gromacs.grompp(
            f=settings.files.mdp_box_npt, c=f'nvt_eq{suffix}.gro', p=topology,
            o=f'npt{suffix}.tpr', maxwarn=2,
        )
        gromacs.mdrun(
            '-v', deffnm=f'npt{suffix}',
            ntmpi=self.ntmpi, ntomp=self.ntomp,
            pin=self.pin, pinoffset=self.pinoffset,
        )

        # Make sure the density has settled before production, then that the
        # equilibrated box is still large enough for the QM/MM cluster
        self.ensure_density_converged(settings, suffix)
        self.check_box_size(settings, f'npt{suffix}.gro')
################################################################################
    @staticmethod
    def read_gro_box(gro_path):
        """Box vector lengths (nm) from the last line of a .gro file."""
        with open(gro_path) as f:
            last = [line for line in f if line.strip()][-1]
        return [float(v) for v in last.split()[:3]]
################################################################################
    def check_box_size(self, settings, gro_path):
        """Stop if the box is smaller than twice the cluster radius.

        The cluster around the central molecule must fit inside the box with
        the minimum-image convention, i.e. L >= 2 * cluster radius.
        """
        rc = settings.advanced.configuration.cluster_radius_nm
        box = self.read_gro_box(gro_path)
        L = min(box)
        print(f'Box check ({gro_path}): smallest box length {L:.3f} nm, '
              f'2 x cluster radius = {2 * rc:.3f} nm')
        if L < 2 * rc:
            raise RuntimeError(
                f"Box too small for the cluster: smallest box length {L:.3f} nm in "
                f"{gro_path} is less than 2 x cluster_radius_nm = {2 * rc:.3f} nm. "
                f"Increase box_length_nm (more molecules) or reduce cluster_radius_nm "
                f"in Settings.yml.")
################################################################################
    @staticmethod
    def density_drift(edr_path, fraction=0.5):
        """Compare the mean density of the two halves of the last `fraction`
        of an NPT run. Returns (first_half_mean, second_half_mean, % difference)."""
        data = pyedr.edr_to_dict(edr_path)
        rho = np.asarray(data['Density'], dtype=float)
        tail = rho[int(len(rho) * (1 - fraction)):]
        if tail.size < 4:
            raise RuntimeError(f'Too few density points in {edr_path} to check convergence')
        a, b = np.array_split(tail, 2)
        m1, m2 = a.mean(), b.mean()
        return m1, m2, 100.0 * abs(m2 - m1) / m2
################################################################################
    def ensure_density_converged(self, settings, suffix=''):
        """Extend the NPT run until the density has stopped drifting.

        Converged when the mean densities of the two halves of the second half
        of the run differ by less than density_tolerance_pct.
        """
        eq = settings.advanced.equilibration
        if not eq.density_check:
            return
        name = f'npt{suffix}'
        for attempt in range(eq.max_extensions + 1):
            m1, m2, pct = self.density_drift(f'{name}.edr')
            print(f'NPT density check: {m1:.2f} -> {m2:.2f} kg/m^3 '
                  f'({pct:.2f}% change, tolerance {eq.density_tolerance_pct}%)')
            if pct <= eq.density_tolerance_pct:
                print('NPT density converged.')
                return
            if attempt == eq.max_extensions:
                break
            print(f'Density not converged; extending NPT by {eq.extend_ps} ps '
                  f'(extension {attempt + 1} of {eq.max_extensions}).')
            gromacs.convert_tpr(s=f'{name}.tpr', extend=eq.extend_ps, o=f'{name}.tpr')
            gromacs.mdrun(
                '-v', deffnm=name, cpi=f'{name}.cpt',
                ntmpi=self.ntmpi, ntomp=self.ntomp,
                pin=self.pin, pinoffset=self.pinoffset,
            )
        raise RuntimeError(
            f'NPT density still drifting after {eq.max_extensions} extensions '
            f'({pct:.2f}% > {eq.density_tolerance_pct}%). Check the system, or '
            f'raise max_extensions / extend_ps in Settings.yml.')
################################################################################
    def create_mdpfile(self, HOMEDIR, mdpfile, T, p, template='template_2.mdp'):

        T_kelvin = T  # temperature / K
        p_bar = p  # pressure / bar
        replace_dict = {'TEMPERATURE': f'{T_kelvin}',
                        'PRESSURE': f'{p_bar}'}
        search_text = 'TEMPERATURE', 'PRESSURE'
        replace_text = f'{T_kelvin}', f'{p_bar}'
        with open(os.path.join(HOMEDIR, template), 'r') as file: #Note that template 1 is the version with PR barostat while template 2 is the version with the c-resync barostat
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
            topology = os.path.join(HOMEDIR, settings.files.pure_topology)
            runname = 'Pure_QMMM_md3'
            suffix = ''
        elif role == 'mixture':
            topology = os.path.join(HOMEDIR, settings.files.mixture_topology)
            runname = 'Mixture_QMMM_md3'
            suffix = '_solute'
        else:
            raise ValueError(f"role must be 'pure' or 'mixture', got {role!r}")
    
        mdpfile='Prod.mdp'
        self.create_mdpfile(HOMEDIR, mdpfile, T, p, template=settings.files.mdp_production_template)
        print('performing molecular dynamics')
        grofile = os.path.join(HOMEDIR, f'npt{suffix}.gro')
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
