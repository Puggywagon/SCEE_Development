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
class Pre_Eq_Simulations(object):
    def __init__(self):
        hostname = socket.gethostname()

        # basic GROMACS threading setup
        self.ntmpi = 1
        self.pin = 'auto'
        self.pinoffset = 0   # we'll keep it simple for now
        if 'Tower3' in hostname:
            self.ntomp = 8
            self.max_heavy_jobs = 4 # max mdrun processes on Tower3
        else:
            self.ntomp = 6
            self.max_heavy_jobs = 3   # max mdrun processes on office PC

        # if True, wait until *no* Gaussian jobs are running before starting MD
        self.respect_gaussian = True

################################################################################
    def _wait_for_md_slot(self):
        while True:
            try:
                md_count = int(subprocess.check_output(['pgrep', '-c', 'mdrun']))
            except subprocess.CalledProcessError:
                md_count = 0

            try:
                g09_count = int(subprocess.check_output(['pgrep', '-c', 'g09']))
            except subprocess.CalledProcessError:
                g09_count = 0

            # combined heavy load: MD + Gaussian
            if (md_count + g09_count) < self.max_heavy_jobs:
                break

            time.sleep(10)
################################################################################
    def Pre_Eq_Solute(self,Gro_File,Topology_File,L):
        L=L+1.7
        print('performing molecular dynamics')
        gromacs.check(c=Gro_File,)
        gromacs.editconf(f=Gro_File,, box=[L,L,L], o=Gro_File,)
        gromacs.check(c=Gro_File,)
        
        # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
        gromacs.grompp(f='minim_vacuum.mdp', c=Gro_File, p=Topology_File, o='em.tpr', maxwarn=2)
    
        # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
        self._wait_for_md_slot()
        gromacs.mdrun('-v', deffnm='em', ntmpi=self.ntmpi, ntomp=self.ntomp, pin=self.pin, pinoffset=self.pinoffset)
        
        # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
        gromacs.grompp(f='nvt_vacuum.mdp', c='em.gro', p=Topology_File, o='nvt_vacuum2.tpr', maxwarn=2)
    
        # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
        self._wait_for_md_slot()
        gromacs.mdrun('-v',deffnm='nvt_vacuum2', ntmpi=self.ntmpi, ntomp=self.ntomp, pin=self.pin, pinoffset=self.pinoffset)
        
################################################################################
    def insert_Molecules(self,Gro_File,initial_molecules):
        Solute=f'nvt_vacuum2.gro'
        Solvent=Gro_File
        gromacs.insert_molecules(f=Solute,ci=Solvent, nmol=f'{initial_molecules}', o='out.gro')
################################################################################
    def Write_to_top(self,Topology_File,initial_molecules,solresnametop):
        with open(Topology_File, 'a') as file:
            file.write(f'{solresnametop}            {initial_molecules}\n')
################################################################################
    def Pre_Eq_System(self,Topology_File):
        # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
        gromacs.grompp(f='minim.mdp', c='out.gro',p=Topology_File, o='em1.tpr', maxwarn=2)
    
        # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
        self._wait_for_md_slot()
        gromacs.mdrun('-v', deffnm='em1', ntmpi=self.ntmpi, ntomp=self.ntomp, pin=self.pin, pinoffset=self.pinoffset)
        
        # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
        gromacs.grompp(f='nvt.mdp', c='em1.gro', p=Topology_File, o=f'nvt_eq.tpr', maxwarn=2)
    
        # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
        self._wait_for_md_slot()
        gromacs.mdrun('-v', deffnm=f'nvt_eq', ntmpi=self.ntmpi, ntomp=self.ntomp, pin=self.pin, pinoffset=self.pinoffset)     
        
        # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
        gromacs.grompp(f='npt.mdp', c='nvt_eq.gro', p=Topology_File, o=f'system.tpr', maxwarn=2)
    
        # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
        self._wait_for_md_slot()
        gromacs.mdrun('-v', deffnm=f'system', ntmpi=self.ntmpi, ntomp=self.ntomp, pin=self.pin, pinoffset=self.pinoffset)

################################################################################
    def get_dipole_model(self): # Think we move this one and the next one into a different section in the script. Maybe we combine the MD and SCEE loops so that it is all happening in one replica folder?
        gromacs.environment.flags['capture_output'] = True
        input_str = '0\n'
        tmp = gromacs.tools.Dipoles(f=f'nvt_vacuum2',
                                s=f'nvt_vacuum2',
                                input=input_str)
        chk, dipole_output, stderr = tmp.run()

        f = open(f'Dipole.txt', 'w')
        f.write(dipole_output)
        f.close()

        dipole_lines = dipole_output.splitlines()  # Split text into lines
        dipole_line = dipole_lines[8]  # Assuming epsilon is on the last line
        dipole_model = float(dipole_line.split()[2])  # Extract the numerical value
        return dipole_model
######################################################
    def get_dipole_model_liquid(self): # Add in 
        gromacs.environment.flags['capture_output'] = True
        input_str = '0\n'
        tmp = gromacs.tools.Dipoles(f=f'Pure_QMMM_md3.trr',
                                s=f'Pure_QMMM_md3.tpr',
                                input=input_str)
        chk, dipole_output, stderr = tmp.run()

        f = open(f'Dipole2.txt', 'w')
        f.write(dipole_output)
        f.close()

        dipole_lines = dipole_output.splitlines()  # Split text into lines
        dipole_line = dipole_lines[8]  # Assuming epsilon is on the last line
        Model_Dipole = float(dipole_line.split()[2])  # Extract the numerical value
        epsilon_line = dipole_lines[-1]  # Assuming epsilon is on the last line
        epsilon = float(epsilon_line.split()[2])  # Extract the numerical value
        
        return Model_Dipole, Epsilon
######################################################
