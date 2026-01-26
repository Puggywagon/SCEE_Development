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
class Simulations_Analysis(object):
    def __init__(self):
        pass
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
