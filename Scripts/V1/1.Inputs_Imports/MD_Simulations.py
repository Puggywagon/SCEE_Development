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
class MD_Simulations(object):
    def __init__(self):
        hostname = socket.gethostname()

        self.ntmpi = 1
        self.pin = 'auto'
        self.pinoffset = 0

        if 'Tower3' in hostname:
            self.ntomp = 8
            self.max_heavy_jobs = 4
        else:
            self.ntomp = 6
            self.max_heavy_jobs = 3

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
    def run_md(self,mdpfile,HOMEDIR,System_Gro,Topology_File):
        print('performing molecular dynamics')
        
        if pure_solvent == 'Yes':
            system_title='Pure'
        else:
            system_title='Mixture'
  
        mdpfile = f'{mdpfile}'
        grofile = HOMEDIR+f'/{system_title}.gro'
        topol = HOMEDIR+f'/{system_title}.top'
        print(f'mdpfile:{mdpfile} \n grofile: {grofile}\n topol: {topol}')

        runname = f'{system_title}_QMMM_md3'
        edrfile = f'{runname}.edr'
        groout = f'{runname}.gro'
        xtcfile = f'{runname}.xtc'
        trrfile = f'{runname}.trr'
        tprfile = f'{runname}.tpr'
        # gmx grompp -f md.mdp -c argon_start.pdb -p argon.top
        gromacs.grompp(f=mdpfile, c=grofile, p=topol, o=tprfile, maxwarn=2)
        # gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0
        self._wait_for_md_slot()
        gromacs.mdrun('-v', s=tprfile, c=groout, o=trrfile, x=xtcfile, e=edrfile, ntmpi=self.ntmpi, ntomp=self.ntomp, pin=self.pin, pinoffset=self.pinoffset)
        return System_title
################################################################################
    def create_mdpfile(self,HOMEDIR,mdpfile,T,p):
        T_kelvin = T  # temperature / K
        p_bar = p  # pressure / bar
        replace_dict = {'TEMPERATURE': f'{T_kelvin}',
                        'PRESSURE': f'{p_bar}'}
        search_text = 'TEMPERATURE', 'PRESSURE'
        replace_text = f'{T_kelvin}', f'{p_bar}'
        with open(HOMEDIR + '/template_2.mdp', 'r') as file:
            data = file.read()
            for search_text, replace_text in replace_dict.items():
                data = data.replace(search_text, replace_text)
        with open(mdpfile, 'w') as file:
            file.write(data)
################################################################################
    def process_trajectory(self,system_title):
        print('processing trajectory')
        input_str = '0\n'
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
##############################################################################################

    
