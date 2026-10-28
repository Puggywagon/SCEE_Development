#!/usr/bin/python3
import gromacs
import subprocess
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
################################################################################
class Pre_Eq_Simulations(object):
    def __init__(self):
        pass
################################################################################
    def Pre_Eq_Solute(self,central,L):
        L=L+1.7
        print('performing molecular dynamics')
        gromacs.check(c=f'{central}_AA.gro')
        gromacs.editconf(f=f'{central}_AA.gro', box=[L,L,L], o=f'{central}_AA.gro')
        gromacs.check(c=f'{central}_AA.gro')
        
        # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
        gromacs.grompp(f='minim_vacuum.mdp', c=f'{central}_AA.gro', p=f'{central}_AA.top', o='em.tpr', maxwarn=2)
    
        # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
        gromacs.mdrun('-v', deffnm='em',ntmpi=1, ntomp=6, pin='on', pinoffset=12)
        
        # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
        gromacs.grompp(f='nvt_vacuum.mdp', c='em.gro', p=f'{central}_AA.top', o='nvt_vacuum2.tpr', maxwarn=2)
    
        # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
        gromacs.mdrun('-v',deffnm='nvt_vacuum2', ntmpi=1, ntomp=6, pin='on', pinoffset=12)
################################################################################
    def insert_Molecules(self,central,solvent,system_title,initial_molecules):
        Solute=f'nvt_vacuum2.gro'
        Solvent=f'{solvent}_UA.gro'
        gromacs.insert_molecules(f=Solute,ci=Solvent, nmol=f'{initial_molecules}', o='out.gro')
################################################################################
    def Write_to_top(self,system_title,initial_molecules,solresnametop):
        with open(f'{system_title}.top', 'a') as file:
            file.write(f'{solresnametop}            {initial_molecules}\n')
################################################################################
    def Pre_Eq_System(self,system_title):
        # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
        gromacs.grompp(f='minim.mdp', c='out.gro',p=f'{system_title}.top', o='em1.tpr', maxwarn=2)
    
        # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
        gromacs.mdrun('-v', deffnm='em1', ntmpi=1, ntomp=6, pin='on', pinoffset=12)
        
        # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
        gromacs.grompp(f='nvt.mdp', c='em1.gro', p=f'{system_title}.top', o=f'nvt_eq.tpr', maxwarn=2)
    
        # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
        gromacs.mdrun('-v', deffnm=f'nvt_eq',ntmpi=1, ntomp=6, pin='on', pinoffset=12)        
        
        # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
        gromacs.grompp(f='npt.mdp', c='nvt_eq.gro', p=f'{system_title}.top', o=f'{system_title}.tpr', maxwarn=2)
    
        # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
        gromacs.mdrun('-v', deffnm=f'{system_title}', ntmpi=1, ntomp=6, pin='on', pinoffset=12)
################################################################################
    def get_dipole_model(self):
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
    def get_dipole_model_liquid(self,system_title):
        gromacs.environment.flags['capture_output'] = True
        input_str = '0\n'
        tmp = gromacs.tools.Dipoles(f=f'Simulations/replica_1/298.0K/1.0Bar/{system_title}_QMMM_md3',
                                s=f'Simulations/replica_1/298.0K/1.0Bar/{system_title}_QMMM_md3',
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
        
        return Model_Dipole, epsilon
######################################################
