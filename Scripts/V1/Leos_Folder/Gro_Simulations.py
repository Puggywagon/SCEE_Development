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
    def insert_Molecules(self,initial_molecules,Mixture):
        if Mixture==No:
            Solute=f'nvt_vacuum2.gro'
            Adding=f'nvt_vacuum2.gro'          
            gromacs.insert_molecules(f=Solute,ci=Adding, nmol=f'{initial_molecules}', o='out.gro')
        elif Mixture==Yes:
            Solute=f'system.gro'
            Adding=f'nvt_vacuum3.gro'
            gromacs.insert_molecules(f=Solute,ci=Adding, nmol=f'{initial_molecules}', o='out2.gro')
################################################################################
    def Write_to_top(self,Topology_File,initial_molecules,solresnametop):
        with open(Topology_File, 'a') as file:
            file.write(f'{solresnametop}            {initial_molecules}\n')
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
    def process_trajectory(self,system_title,Mixture):
        print('processing trajectory')
        input_str = '0\n'
        if Mixture == 'Yes':
            tmp = gromacs.tools.Trjconv(f=f'{system_title}_QMMM_md3_2.xtc',
                                s=f'{system_title}_QMMM_md3_2.tpr',
                                o=f'conf2_.gro',
                                input=input_str,
                                pbc='whole',
                                b=1000, dt=20, sep=True)
            chk, stdout, stderr = tmp.run()
            
        elif Mixture == 'No':
            tmp = gromacs.tools.Trjconv(f=f'{system_title}_QMMM_md3.xtc',
                                s=f'{system_title}_QMMM_md3.tpr',
                                o=f'conf_.gro',
                                input=input_str,
                                pbc='whole',
                                b=1000, dt=20, sep=True)
            chk, stdout, stderr = tmp.run()
################################################################################
    def process_gro(self):
        print('processing conf_*.gro files')
        logfile = open('junk.log', 'w')
        cmd5 = ["gfortran -o Shell_Oniom Shell_Oniom.f90 && ./Shell_Oniom"]
        subprocess.call(cmd5, shell=True, stdout=logfile, stderr=logfile)
################################################################################
    def run_md(self,Gro_File,Topology_File,L,initial_molecules,solresnametop,Mixture, MD='Vacuum'):
        L=L+1.7
        print('performing molecular dynamics')
        gromacs.check(c=Gro_File,)
        gromacs.editconf(f=Gro_File, box=[L,L,L], o=Gro_File,)
        gromacs.check(c=Gro_File,)
        if Mixture == 'Yes':
            # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
            gromacs.grompp(f='minim_vacuum.mdp', c=Gro_File, p=Topology_File, o='em.2tpr', maxwarn=2)
    
            # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
            self._wait_for_md_slot()
            gromacs.mdrun('-v', deffnm='em2', ntmpi=self.ntmpi, ntomp=self.ntomp, pin=self.pin, pinoffset=self.pinoffset)
        
            # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
            gromacs.grompp(f='nvt_vacuum.mdp', c='em2.gro', p=Topology_File, o='nvt_vacuum3.tpr', maxwarn=2)
    
            # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
            self._wait_for_md_slot()
            gromacs.mdrun('-v',deffnm='nvt_vacuum3', ntmpi=self.ntmpi, ntomp=self.ntomp, pin=self.pin, pinoffset=self.pinoffset)
        
        elif Mixture == 'No':
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
        
        self.insert_Molecules(Gro_File,initial_molecules,Mixture)
        self.Write_to_top(Topology_File,initial_molecules,solresnametop)
################################################################################
    def run_md(self,Topology_File,Mixture, MD='Box'):
        if Mixture== 'Yes':
            # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
            gromacs.grompp(f='minim.mdp', c='out2.gro',p=Topology_File, o='em3.tpr', maxwarn=2)
    
            # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
            self._wait_for_md_slot()
            gromacs.mdrun('-v', deffnm='em3', ntmpi=self.ntmpi, ntomp=self.ntomp, pin=self.pin, pinoffset=self.pinoffset)
        
            # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
            gromacs.grompp(f='nvt.mdp', c='em3.gro', p=Topology_File, o=f'nvt_eq2.tpr', maxwarn=2)
    
           # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
            self._wait_for_md_slot()
            gromacs.mdrun('-v', deffnm=f'nvt_eq2', ntmpi=self.ntmpi, ntomp=self.ntomp, pin=self.pin, pinoffset=self.pinoffset)     
        
            # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
            gromacs.grompp(f='npt.mdp', c='nvt_eq2.gro', p=Topology_File, o=f'system2.tpr', maxwarn=2)
    
            # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
            self._wait_for_md_slot()
            gromacs.mdrun('-v', deffnm=f'system2', ntmpi=self.ntmpi, ntomp=self.ntomp, pin=self.pin, pinoffset=self.pinoffset)
       
        elif Mixture == 'No':
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
    def run_md(self,mdpfile,HOMEDIR,system_title,T,p,Mixture, MD='Production'):
        self.create_mdpfile(HOMEDIR,mdpfile,T,p)
        print('performing molecular dynamics')

        mdpfile = f'{mdpfile}'
        if Mixture=='Yes':
            grofile = HOMEDIR+f'/system2.gro'
        elif Mixture == 'No':
            grofile = HOMEDIR+f'/system.gro'
        topol = HOMEDIR+f'/{system_title}.top'
        print(f'mdpfile:{mdpfile} \n grofile: {grofile}\n topol: {topol}')

        runname = f'{system_title}_QMMM_md3'
        edrfile = f'{runname}.edr'
        groout = f'{runname}.gro'
        xtcfile = f'{runname}.xtc'
        trrfile = f'{runname}.trr'
        tprfile = f'{runname}.tpr'
        # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
        
        if Mixture == 'Yes':
            gromacs.grompp(f=mdpfile, c=grofile, p=topol, o=tprfile+'_2', maxwarn=2)
            # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
            self._wait_for_md_slot()
            gromacs.mdrun('-v', s=tprfile+'_2', c=groout+'_2', o=trrfile+'_2', x=xtcfile+'_2', e=edrfile+'_2', ntmpi=self.ntmpi, ntomp=self.ntomp, pin=self.pin, pinoffset=self.pinoffset)
        
        elif Mixture == 'No':
            gromacs.grompp(f=mdpfile, c=grofile, p=topol, o=tprfile, maxwarn=2)
            # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
            self._wait_for_md_slot()
            gromacs.mdrun('-v', s=tprfile, c=groout, o=trrfile, x=xtcfile, e=edrfile, ntmpi=self.ntmpi, ntomp=self.ntomp, pin=self.pin, pinoffset=self.pinoffset)
        
        self.process_trajectory(system_title,Mixture)
        self.process_gro(Mixture)
################################################################################
